#====================================
#FSGEA with BTM ----
#====================================
library(dplyr)
library(fgsea)
library(data.table)
run_GSEA_BTM <- function(deg_csv_path,
                         output_csv_path,
                         gmt_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data_24082026/BloodGen3_modules_v2.gmt",
                         rank_column = "logFC",
                         gene_column = "HGNC_symbol",
                         fdr_cutoff = 0.01,
                         minSize = 10,
                         maxSize = 500,
                         nperm = 10000) {
  
  BTM <- fgsea::gmtPathways(gmt_path)
  deg_df <- data.table::fread(deg_csv_path, data.table = FALSE, fill= TRUE)
  # Validate required columns
  if (!all(c(rank_column, gene_column) %in% colnames(deg_df))) {
    stop(paste("Missing required columns in input file:", 
               paste(setdiff(c(rank_column, gene_column), colnames(deg_df)), collapse = ", ")))
  }
  
  #named vector for rankings
  #if abs logFC is to be used
  # deg_df$abs_rank_value <- abs(deg_df[[rank_column]])
  # rankings <- setNames(deg_df$abs_rank_value, deg_df[[gene_column]])
  
  rankings <- setNames(deg_df[[rank_column]], deg_df[[gene_column]])
  rankings <- sort(rankings, decreasing = TRUE)
  
  #Run GSEA
  gsea_res <- fgsea::fgsea(
    pathways = BTM,
    stats = rankings,
    scoreType = "std", #change to pos for only positive enrichment 
    minSize = minSize,
    maxSize = maxSize,
    nproc = 1,
    nPermSimple = nperm
  )
  
  #change the leadingEdge column which is a list of 
  #characters into character vector with semicolon
  gsea_res$leadingEdge <- vapply(
    gsea_res$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )
  
  #Order and filter by adjusted p-value
  gsea_res <- gsea_res[order(gsea_res$padj), ]
  gsea_fdr <- dplyr::filter(gsea_res, padj < fdr_cutoff)
  data.table::fwrite(gsea_fdr, output_csv_path)
  return(gsea_fdr)
}

# fgsea_output <- run_GSEA_BTM(
#   deg_csv_path = "/home/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/FGSEA_logFC/Washington_spearcorr.csv",
#   output_csv_path = paste0("FGSEA",basename(deg_csv_path),".csv"))

deg_files <- list.files(
  path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/aggregate_DEG_kallisto_filterbyExpr_groupspecific",
  pattern = "\\.csv$",
  full.names = TRUE
)

gsea_outputs <- lapply(deg_files, function(deg_csv) {
  output_file <- paste0("FGSEA_", basename(deg_csv))
  run_GSEA_BTM(
    deg_csv_path = deg_csv,
    output_csv_path = output_file
  )
})

#Make a consolidated table with the pathways and scores for each gene-ranked list
library(purrr)
library(tidyr)
library(stringr)
library(tibble)

FGSEA_results = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/FGSEA_aggregate_kallisto_filterbyExpr_groupspecific"

csv_files = list.files(FGSEA_results,
  pattern = "\\.csv$",
  full.names = TRUE)

read_fgsea_nes = function(file_path) {
  df = fread(file_path, fill = TRUE,data.table = FALSE)
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # ensure required columns exist
  if (!all(c("pathway", "NES") %in% colnames(df))) {
    return(NULL)
  }
  
  # remove blank pathways
  df = df %>%
    filter(!is.na(pathway) & trimws(pathway) != "")
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  df_out = df %>%
    dplyr::select(pathway, NES) %>%
    dplyr::rename(!!basename(file_path) := NES)
  return(df_out)
}

nes_list_raw = lapply(csv_files,read_fgsea_nes)

nes_list = nes_list_raw[!sapply(nes_list_raw, is.null)]

nes_matrix = reduce( nes_list, full_join,by = "pathway")
# clean column names
colnames(nes_matrix)[-1] = str_remove(  colnames(nes_matrix)[-1], "\\.csv$")
colnames(nes_matrix) = stringr::str_remove(colnames(nes_matrix),"^FGSEA_")
colnames(nes_matrix) =stringr::str_remove(colnames(nes_matrix), "_DEG_results")
nes_matrix <- nes_matrix[!grepl("^TBD", nes_matrix$pathway), ]
fwrite(nes_matrix,file = "NES_matrix_all_fgsea.csv")

rownames(nes_matrix) = NULL
#plotting heatmap for the commmon pathways with NES scores for each condition
nes_matrix = column_to_rownames(nes_matrix, var = "pathway")
nes_matrix = as.matrix(nes_matrix)
nes_matrix[is.na(nes_matrix)] = 0


library(ComplexHeatmap)
library(circlize)
heatmap_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

ht <- Heatmap(
  nes_matrix,
  name = "NES",
  col = heatmap_colors,
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  border = TRUE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "left",      
  row_dend_side = "right",      
  row_names_max_width = unit(15, "cm"),
  row_names_gp = gpar(fontsize = 13),
  column_names_gp = gpar(fontsize = 13),
  column_names_rot = 45,
  
  heatmap_legend_param = list(
    title = "Normalized Enrichment Score",
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 12)
  )
)

png("NES_heatmap_BTM_strainwise_bloodgen3_HRLR.png", width = 18, height = 18, units = "in", res = 300)
draw(
  ht,
  heatmap_legend_side = "right",
  padding = unit(c(5, 10, 5, 50), "mm")
)
dev.off()

# create a fontsize vector — default 14 for all rows
fontsizes <- rep(14, nrow(nes_matrix))
fontcolors <- rep("black", nrow(nes_matrix))
# set specific rows to bigger font by name
big_rows <- c("Neutrophil activation(M10.4)", "B cells(M12.8)",
              "B cells(M16.57)","T cells(M12.6)","Lymphocytes(M13.27)",
              "Erythroid cells(M9.2)","Cytotoxic lymphocytes(M9.1)",
              "Type 1 Interferon(M8.3)","Interferon(M10.1)", 
              "Inflammation(M13.12)", "Interferon(M13.17)",
              "Cytokines/chemokines(M13.16)")  
fontsizes[rownames(nes_matrix) %in% big_rows] <- 18
fontcolors[rownames(nes_matrix) %in% big_rows] <- "blue"
ht <- Heatmap(
  nes_matrix,
  name = "NES",
  col = heatmap_colors,
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  border = TRUE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "left",
  row_dend_side = "right",
  row_names_max_width = unit(15, "cm"),
  
  row_names_gp = gpar(fontsize = fontsizes, col = fontcolors),
  column_names_gp = gpar(fontsize = 14),
  column_names_rot = 45,
  
  heatmap_legend_param = list(
    title = "Normalized Enrichment Score",
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 10)
  )
)

png("NES_agg_v1.png", width = 18, height = 18, units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", padding = unit(c(5, 10, 5, 5), "mm"))
dev.off()
