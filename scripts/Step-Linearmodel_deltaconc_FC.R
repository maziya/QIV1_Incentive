library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

foldchange  =read.table('/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv', sep=',',
                            header = TRUE)

foldchange = foldchange[,c(1,14,15,16,17)]

strains <- c("FC.A.HongKong", "FC.A.Victoria", "FC.B.Phuket", "FC.B.Washington")

run_FCcorr <- function(matrix_V1, matrix_V2, metadata, strains, feature_type = "feature", include_batch = TRUE) {
  
  colnames(matrix_V1) <- gsub("V1$", "", colnames(matrix_V1))
  colnames(matrix_V2) <- gsub("V2$", "", colnames(matrix_V2))
  
  valid_subjects <- Reduce(intersect, list(colnames(matrix_V1), colnames(matrix_V2), metadata$SubjectID))
  
  mat1 <- matrix_V1[, valid_subjects]
  mat2 <- matrix_V2[, valid_subjects]
  meta <- metadata[match(valid_subjects, metadata$SubjectID), ]
  
  # Calculate Delta 
  delta_mat <- as.matrix(mat2 - mat1)
  
  # Run models across all strains
  all_results <- lapply(strains, function(strain) {
    
    # Pre-extract the specific response vector to speed up the loop
    response_vec <- meta[[strain]]
    
    res_list <- lapply(1:nrow(delta_mat), function(i) {
      x <- delta_mat[i, ]
      
      # Base dataframe (without batch)
      df <- data.frame(response = response_vec, delta = x, sex = meta$SEX, age = meta$AGE)
      
      # Conditional Linear Model based on include_batch
      if (include_batch) {
        df$batch <- meta$Batch
        fit <- lm(response ~ delta + sex + batch + age, data = df)
      } else {
        fit <- lm(response ~ delta + sex + age, data = df)
      }
      
      coef <- summary(fit)$coefficients
      
      # Correlation
      test <- cor.test(x, response_vec, method = "spearman", exact = FALSE, use = "complete.obs")
      
      c(beta = as.numeric(coef["delta", "Estimate"]),
        se   = as.numeric(coef["delta", "Std. Error"]),
        t    = as.numeric(coef["delta", "t value"]),
        pval = as.numeric(coef["delta", "Pr(>|t|)"]),
        cor  = test$estimate,
        cor_pval = test$p.value)
    })
    
    # Compile results for this strain
    res_df <- as.data.frame(do.call(rbind, res_list))
    res_df[[feature_type]] <- rownames(delta_mat)
    res_df$strain <- strain
    res_df$padj <- p.adjust(res_df$pval, method = "BH")
    res_df$cor_padj <- p.adjust(res_df$cor_pval, method = "BH")
    
    return(res_df)
  })
  
  return(do.call(rbind, all_results))
}


# Run for Proteomics 
metadata = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/untargprot_metadata_FC.csv")

untargeted_matrix_log = read.csv('/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/proteomics_log.csv', row.names = 1)
prot_V1 <- untargeted_matrix_log[, grep("V1$", colnames(untargeted_matrix_log))]
prot_V2 <- untargeted_matrix_log[, grep("V2$", colnames(untargeted_matrix_log))]


prot_results <- run_FCcorr(
  matrix_V1 = prot_V1, 
  matrix_V2 = prot_V2, 
  metadata = metadata,
  strains = strains,
  feature_type = "protein",
  include_batch = TRUE
)

# Run for Metabolomics

metadata = read.csv('/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/mets_metadata_FC.csv')
metab = read.csv('/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/log_batchadj_scaled_metabolomics.csv', row.names = 1)
metabolomics_V1 <-metab[, grep("V1$", colnames(metab))]
metabolomics_V2 <- metab[, grep("V2$", colnames(metab))]

met_results <- run_FCcorr(
  matrix_V1 = metabolomics_V1, 
  matrix_V2 = metabolomics_V2, 
  metadata = metadata,
  strains = strains,
  feature_type = "metabolite",
  include_batch = FALSE
)

# Plots
plot_volcano <- function(results_df, feature_col, mapping_df = NULL) {
  for (s in unique(results_df$strain)) {
    df <- results_df %>% 
      filter(strain == s) %>%
      mutate(
        log10p = -log10(pval),
        significant = pval < 0.05
      )
    
    top_hits <- df %>% arrange(pval) %>% head(10)
    
    # If a mapping dataframe is provided (like for proteins to HGNC), join it
    if (!is.null(mapping_df)) {
      top_hits <- inner_join(top_hits, mapping_df, by = setNames(feature_col, feature_col))
      label_col <- "HGNC_symbol"
    } else {
      label_col <- feature_col
    }
    
    p <- ggplot(df, aes(x = beta, y = log10p, color = significant)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      scale_color_manual(values = c("grey", "blue")) +
      geom_text_repel(data = top_hits, aes(label = .data[[label_col]]), size = 3) +
      theme_bw() +
      labs(title = s, x = "Effect Size (Beta)", y = "-log10(p-value)")
    
    ggsave(paste0("volcano_", feature_col, "_", s, ".png"), plot = p, height = 8, width = 8)
  }
}
untargeteduniprot_ENSG_HGNC <- read_csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/untargeteduniprot_ENSG_HGNC.csv")
colnames(untargeteduniprot_ENSG_HGNC)[1] = "protein"
plot_volcano(prot_results, "protein", untargeteduniprot_ENSG_HGNC)
plot_volcano(met_results, "metabolite")

# ==============================================================================
# Significant proteins and mets
# ==============================================================================
prot_results_annotated <- prot_results %>% 
  left_join(untargeteduniprot_ENSG_HGNC, by = "protein")

lapply(strains, function(s) {
  prot_results_annotated %>% 
    filter(strain == s) %>% 
    write.csv(paste0(s, "_prot_lm_FC.csv"), row.names = FALSE, quote = FALSE)
})

lapply(strains, function(s) {
  met_results %>% filter(strain == s) %>% 
    write.csv(paste0(s, "_mets_lm_FC.csv"), row.names = FALSE, quote = FALSE)
})


library(dplyr)
library(readr)

# sig proteins and mets from either the linear model or the corr coeff analysis
get_significant_hits <- function(file_paths) {
  significant_data <- lapply(file_paths, read_csv, show_col_types = FALSE) %>%
    bind_rows() %>%
    filter(pval < 0.05 & cor_pval < 0.05)
  
  return(significant_data)
}


prot_files <- list.files(pattern = "_prot_lm_FC.csv", full.names = TRUE)
mets_files <- list.files(pattern = "_mets_lm_FC.csv", full.names = TRUE)

sig_proteins <- get_significant_hits(prot_files) #took "&" condition of pval- cor_pval
sig_metabolites <- get_significant_hits(mets_files) #took "|" condition of pval - cor_pval

