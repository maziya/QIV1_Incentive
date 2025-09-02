library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tibble)
library(readxl)
library(edgeR)
library(limma)
library(ggplot2)


#mapping to ensembl genes and target ids
load_gene_annotations <- function() {
  Gn <- genes(EnsDb.Hsapiens.v86, columns = c("gene_id", "gene_name"))
  Gn <- as_tibble(Gn)
  Gn <- dplyr::rename(Gn, target_id = gene_id)
  Gn <- dplyr::select(Gn, "target_id", "gene_name")
  return(Gn)
}


#Load and filter responder data
# @param file_path Path to responder CSV file
# @param strain_column Column name for the strain (e.g., "RES4.0.B.Phuket")
# @return filtered responder data
load_responder_data = function(file_path, strain_column) {
  responder_data = read.table(file_path, sep = ',', header = TRUE)
  filtered_data = responder_data %>%
    dplyr::filter(!!sym(strain_column) == 'TRUE')
  return(filtered_data)
}

#' Load and prepare metadata
#' @param file_path Path to Excel metadata file
#' @param responder_ids Vector of responder subject IDs
#' @param subjects_to_exclude Vector of subject IDs to exclude
#' @return processed metadata
prepare_metadata <- function(file_path, responder_ids, subjects_to_exclude = NULL) {
  meta <- read_excel(file_path, sheet = 1)
  meta_filt <- meta[, c(2, 3, 6, 7, 30)]
  colnames(meta_filt)[1] <- "SubjectID"
  
  meta_filt <- meta_filt %>%
    mutate(SubjectIDNew = paste0(SubjectID, Visit))
  
  #Remove excluded subjects
  if (!is.null(subjects_to_exclude)) {
    meta_filt <- meta_filt %>%
      dplyr::filter(!SubjectID %in% subjects_to_exclude)
  }
  
  #Assign responder status
  meta_filt$Status <- "NonResponder"
  meta_filt$Status[meta_filt$SubjectID %in% responder_ids] <- "Responder"
  
  return(meta_filt)
}

# ============================================================================
# COUNT MATRIX PROCESSING 
# ============================================================================

#' Process count matrix: remove version IDs, join with gene names, filter duplicates
#' @param count_matrix Count matrix with gene IDs as rownames
#' @param gene_annotations Gene annotation tibble
#' @param gf_id Optional gene filter IDs
#' @return processed count matrix
process_count_matrix <- function(count_matrix,Gn, gf_id = NULL) {
  mtx <- rownames_to_column(count_matrix, var = 'target_id')
  
  #Remove version numbers from Ensembl IDs
  mtx$target_id <- sub("\\..*", "", mtx$target_id)
  
  #Join with HGNC gene symbols
  mtx <- inner_join(mtx,Gn, by = 'target_id')
  
  #Calculate row sums for duplicate filtering
  numeric_data <- mtx[, 2:(ncol(mtx) - 1)] %>% dplyr::select(where(is.numeric))
  row_sums <- rowSums(numeric_data, na.rm = TRUE)
  mtx$rowsums <- row_sums
  
  #Filter duplicates by keeping max rowsum per gene
  mtx_filtered <- mtx %>%
    dplyr::group_by(gene_name) %>%
    dplyr::filter(rowsums == max(rowsums)) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  mtx_filtered <- column_to_rownames(mtx_filtered, var = 'target_id')
  mtx_filtered <- mtx_filtered[, 1:(ncol(mtx_filtered) - 2)]
  result_matrix <- as.matrix(mtx_filtered)
  #only protein coding genes gf_id
  if (!is.null(gf_id)) {
    result_matrix <- result_matrix[rownames(result_matrix) %in% gf_id, , drop = FALSE]
  }
  
  #Remove zero-sum rows
  row_sums <- rowSums(result_matrix)
  result_matrix <- result_matrix[which(row_sums != 0), ]

  return(result_matrix)
}

#' Filter and order count matrix by metadata samples
#' @param count_matrix Processed count matrix
#' @param metadata Metadata tibble
#' @param id_column Column name for sample IDs in metadata
#' @return filtered and ordered count matrix
filter_count_matrix_by_metadata <- function(count_matrix, metadata, id_column = "SubjectID") {
  #Filter columns that exist in metadata
  valid_samples <- colnames(count_matrix)[colnames(count_matrix) %in% metadata[[id_column]]]
  filtered_matrix <- count_matrix[, valid_samples, drop = FALSE]
  
  #Reorder columns to match metadata order
  ordered_matrix <- filtered_matrix[, metadata[[id_column]], drop = FALSE]
  
  return(ordered_matrix)
}

#============================================================================
#DIFFERENTIAL EXPRESSION ANALYSIS 
#============================================================================

#' Create design matrix for DEG analysis
#' @param metadata Metadata tibble
#' @param group_column Column name for main comparison group
#' @param covariates Vector of covariate column names
#' @param group_levels Optional vector specifying factor levels for main group
#' @return design matrix
create_design_matrix <- function(metadata, group_column, covariates = c("SEX", "AGE", "Library_Batch"), group_levels = NULL) {
  # Create main group factor
  if (!is.null(group_levels)) {
    group <- factor(metadata[[group_column]], levels = group_levels)
  } else {
    group <- factor(metadata[[group_column]])
  }
  
  #Create covariate factors/numerics
  design_data <- metadata
  design_data$group <- group
  
  # Convert covariates to appropriate types
  if ("AGE" %in% covariates) {
    design_data$AGE <- as.numeric(design_data$AGE)
  }
  if ("Library_Batch" %in% covariates) {
    design_data$Library_Batch <- factor(design_data$Library_Batch)
  }
  if ("SEX" %in% covariates) {
    design_data$SEX <- factor(design_data$SEX)
  }
  
  # Create formula and design matrix
  covariate_formula <- paste(covariates, collapse = "+")
  design_formula <- as.formula(paste("~ group +", covariate_formula))
  design <- model.matrix(design_formula, data = design_data)
  
  return(design)
}

#' Perform differential expression analysis using limma-voom
#' @param count_matrix Filtered count matrix
#' @param design Design matrix
#' @param coef_name Coefficient name for topTable extraction
#' @param gene_annotations Gene annotation tibble
#' @param block_factor Optional blocking factor for paired analysis (vector of subject IDs)
#' @return DEG results tibble
perform_deg_analysis <- function(count_matrix, design, coef_name, Gn, block_factor = NULL) {
  DGEList <- DGEList(count_matrix)
  DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
  
  if (!is.null(block_factor)) {
    v.DGEList.norm <- voom(DGEList.norm, design, plot = FALSE, normalize.method = "quantile")
    corfit <- duplicateCorrelation(v.DGEList.norm, design, block = block_factor)
    v.DGEList.norm <- voom(DGEList.norm, design, normalize.method = "quantile",
                           block = block_factor,
                           correlation = corfit$consensus.correlation)
    
    fit <- lmFit(v.DGEList.norm, design, block = block_factor,
                 correlation = corfit$consensus.correlation)
  } else {
    v.DGEList.norm <- voom(DGEList.norm, design, plot = FALSE)
    fit <- lmFit(v.DGEList.norm, design)
  }
  
  fit_ebayes <- eBayes(fit, trend = TRUE)
  
  deg_results <- topTable(fit_ebayes, adjust = "BH", coef = coef_name, 
                          number = Inf, sort.by = "logFC")
  deg_results$target_id <- rownames(deg_results)
  
  all_genes <- inner_join(deg_results, Gn, by = "target_id")
  
  return(all_genes)
}

# ============================================================================
# VISUALIZATION 
# ============================================================================

#' Create volcano plot from DEG results
#' @param deg_results DEG results tibble
#' @param title Plot title
#' @param pval_threshold P-value threshold for significance
#' @param fc_threshold Log fold change threshold for significance
#' @return ggplot object
create_volcano_plot <- function(deg_results, title, pval_threshold = 0.05, fc_threshold = 1) {
  deg_results$significance <- ifelse(
    deg_results$P.Value < pval_threshold & abs(deg_results$logFC) > fc_threshold,
    ifelse(deg_results$logFC > fc_threshold, "Upregulated", "Downregulated"),
    "Not significant"
  )
  
  vplot <- ggplot(deg_results) +
    aes(y = -log10(P.Value), x = logFC, text = paste("Symbol:", gene_name), color = significance) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "Upregulated" = "#A020F0",
      "Downregulated" = "#A0522D",
      "Not significant" = "grey60"
    )) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "longdash", colour = "black", size = 1) +
    geom_vline(xintercept = fc_threshold, linetype = "longdash", colour = "black", size = 1) +
    geom_vline(xintercept = -fc_threshold, linetype = "longdash", colour = "black", size = 1) +
    labs(title = title) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  return(vplot)
}

#' Save DEG results and volcano plot
#' @param deg_results DEG results tibble with significance column
#' @param volcano_plot ggplot volcano plot object
#' @param output_prefix Prefix for output files
#' @param width Plot width
#' @param height Plot height
save_deg_outputs <- function(deg_results, volcano_plot, output_prefix, width = 15, height = 15) {
  csv_file <- paste0(output_prefix, ".csv")
  write.table(deg_results, csv_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  png_file <- paste0(output_prefix, ".png")
  ggsave(png_file, plot = volcano_plot, width = width, height = height)
  }

# ============================================================================
# WRAPPER FUNCTIONS
# ============================================================================

#' Complete DEG analysis pipeline for single timepoint comparison
#' @param count_matrix Raw count matrix
#' @param metadata Metadata for samples
#' @param gene_annotations Gene annotation tibble
#' @param comparison_name Name for the comparison (used in outputs)
#' @param group_column Column name for comparison groups
#' @param group_levels Factor levels for main group
#' @param coef_name Coefficient name for extraction
#' @param gf_id Optional gene filter
#' @param covariates Covariates to include
#' @param id_column Sample ID column in metadata
run_single_timepoint_deg <- function(count_matrix, metadata, Gn, 
                                     comparison_name, group_column = "Status", 
                                     group_levels = c("NonResponder", "Responder"),
                                     coef_name = "groupResponder", gf_id = NULL,
                                     covariates = c("SEX", "AGE", "Library_Batch"),
                                     id_column = "SubjectID") {
  
  processed_matrix <- process_count_matrix(count_matrix, Gn, gf_id)
  filtered_matrix <- filter_count_matrix_by_metadata(processed_matrix, metadata, id_column)
  design <- create_design_matrix(metadata, group_column, covariates, group_levels)
  deg_results <- perform_deg_analysis(filtered_matrix, design, coef_name, Gn)
  volcano_plot <- create_volcano_plot(deg_results, comparison_name)
  save_deg_outputs(deg_results, volcano_plot, comparison_name)
  
  return(deg_results)
}

#' Complete DEG analysis pipeline for paired timepoint comparison
#' @param count_matrix Raw count matrix (should contain both timepoints)
#' @param metadata Metadata for samples (should contain both timepoints)
#' @param gene_annotations Gene annotation tibble
#' @param comparison_name Name for the comparison
#' @param visit_column Column name for visit/timepoint
#' @param visit_levels Factor levels for visits
#' @param coef_name Coefficient name for extraction
#' @param gf_id Optional gene filter
#' @param covariates Covariates to include
#' @param id_column Sample ID column in metadata
#' @param block_column Column name for blocking factor (subject ID for paired analysis)
run_paired_timepoint_deg <- function(count_matrix, metadata,Gn,
                                     comparison_name, visit_column = "Visit",
                                     visit_levels = c("V1", "V2"), coef_name = "VisitV2",
                                     gf_id = NULL, covariates = c("SEX", "AGE", "Library_Batch"),
                                     id_column = "SubjectIDNew", block_column = "SubjectID") {
  
  processed_matrix <- process_count_matrix(count_matrix,Gn,gf_id)
  filtered_matrix <- filter_count_matrix_by_metadata(processed_matrix, metadata, id_column)
  design <- create_design_matrix(metadata,visit_column,covariates,visit_levels)
  
  block_factor <- metadata[[block_column]]
  
  deg_results <- perform_deg_analysis(filtered_matrix,design,coef_name,Gn,block_factor)
  volcano_plot <- create_volcano_plot(deg_results,comparison_name)
  save_deg_outputs(deg_results, volcano_plot, comparison_name)
  
  return(deg_results)
}
# ============================================================================
# MAIN ANALYSIS
# ============================================================================

main_analysis <- function() {
  #Gene annotations
  Gn = load_gene_annotations()
  
  #Responder data 
  # responder_file = "/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv"
  responder_file = "C:\\Users\\maziy\\OneDrive - IIT-Madras(IC&SR)\\INCENTIVE\\incentiveR\\2024-05-10_Responder-Non-Responder_v6.1\\QIV1-Responder.csv"
  
  Phuket_res = load_responder_data(responder_file, "RES4.0.B.Phuket")
  
  #Prepare metadata
  # meta_file = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/QIV1_India_MetaData_Transcriptomics.xlsx"
  meta_file = "C:\\Users\\maziy\\OneDrive - IIT-Madras(IC&SR)\\INCENTIVE\\transcriptomics\\QIV1_DEG_Analysis\\data\\QIV1_India_MetaData_Transcriptomics.xlsx"
  excluded_subjects = c('QIV1KEM085', 'QIV1KEM033', 'QIV1KEM037', 
                         'QIV1KEM019', 'QIV1KEM055', 'QIV1KEM064', 'QIV1KEM077')
  
  QIV1_meta_filt <- prepare_metadata(meta_file, Phuket_res$SubjectID, excluded_subjects)
  
  #Split metadata by visit
  QIV1_meta_filt_V1 <- QIV1_meta_filt %>% dplyr::filter(Visit == 'V1')
  QIV1_meta_filt_V2 <- QIV1_meta_filt %>% dplyr::filter(Visit == 'V2')
  QIV1_meta_Res <- QIV1_meta_filt %>% dplyr::filter(Status == 'Responder')
  
  #Count data
  # QIV1count = read.table("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BA088-INCENTIVE_eukaryota_count_table.txt")
  QIV1count = read.table("C:\\Users\\maziy\\OneDrive - IIT-Madras(IC&SR)\\INCENTIVE\\transcriptomics\\QIV1_DEG_Analysis\\data\\BA088-INCENTIVE_eukaryota_count_table.txt")
  #Split count matrix by visit and clean column names
  QIV1_V1 <- QIV1count[, endsWith(colnames(QIV1count), "V2")]
  QIV1_V2 <- QIV1count[, endsWith(colnames(QIV1count), "V3")]
  QIV1_all <- cbind(QIV1_V1, QIV1_V2)
  
  colnames(QIV1_V1) <- gsub("V2$", "", colnames(QIV1_V1))
  colnames(QIV1_V2) <- gsub("V3$", "", colnames(QIV1_V2))
  
  #Fix column names for combined matrix
  colnames_old <- colnames(QIV1_all)
  colnames_new <- colnames_old %>%
    sub("V2$", "V1", .) %>%
    sub("V3$", "V2", .)
  colnames(QIV1_all) <- colnames_new
  
  tmp_env <- new.env()
  load("C:\\Users\\maziy\\OneDrive - IIT-Madras(IC&SR)\\INCENTIVE\\transcriptomics\\QIV1_DEG_Analysis\\data\\protein_coding.RData", envir = tmp_env)
  gf_id = tmp_env$gf_id
  rm(tmp_env)
  
  #V1 Responder vs Non-Responder
  deg_v1 <- run_single_timepoint_deg(
    QIV1_V1, QIV1_meta_filt_V1, Gn,
    "QIV1_DEGs_Phuket_ResvsNR_day0_V1",
    gf_id = gf_id
  )
  
  #V2 Responder vs Non-Responder  
  deg_v2 <- run_single_timepoint_deg(
    QIV1_V2, QIV1_meta_filt_V2, Gn,
    "QIV1_DEGs_Phuket_ResvsNR_day3_V2",
    gf_id = gf_id
  )
  
  #Responders V1 vs V2 
  deg_paired <- run_paired_timepoint_deg(
    QIV1_all, QIV1_meta_Res, Gn,
    "QIV1_DEGs_Phuket_Res_V1vsV2",
    gf_id = gf_id
  )
    return(list(deg_v1 = deg_v1, deg_v2 = deg_v2, deg_paired = deg_paired))
}

#Run
results <- main_analysis()
