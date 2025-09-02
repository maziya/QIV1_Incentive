#====================================
#FSGEA with BTM ----
#====================================
library(dplyr)
library(fgsea)
run_GSEA_BTM <- function(deg_csv_path,
                         output_csv_path,
                         gmt_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BTM_for_GSEA_20131008.gmt",
                         rank_column = "logFC",
                         gene_column = "HGNC_symbol",
                         fdr_cutoff = 0.01,
                         minSize = 10,
                         maxSize = 500,
                         nperm = 10000) {
  
  BTM <- fgsea::gmtPathways(gmt_path)
  deg_df <- read.csv(deg_csv_path)
  # Validate required columns
  if (!all(c(rank_column, gene_column) %in% colnames(deg_df))) {
    stop(paste("Missing required columns in input file:", 
               paste(setdiff(c(rank_column, gene_column), colnames(deg_df)), collapse = ", ")))
  }
  
  #named vector for rankings
  rankings <- setNames(deg_df[[rank_column]], deg_df[[gene_column]])
  rankings <- sort(rankings, decreasing = TRUE)
  
  #Run GSEA
  gsea_res <- fgsea::fgsea(
    pathways = BTM,
    stats = rankings,
    scoreType = "pos", #only positive enrichment
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
  
  write.table(gsea_fdr, file = output_csv_path, row.names = FALSE, sep = ",", quote = FALSE)
  return(gsea_fdr)
}

# fgsea_output <- run_GSEA_BTM(
#   deg_csv_path = "/home/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/FGSEA_logFC/Washington_spearcorr.csv",
#   output_csv_path = paste0("FGSEA",basename(deg_csv_path),".csv"))

deg_files <- list.files(
  path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/FGSEA_spearman_correlation_genes",
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

FGSEA_results =  "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/FGSEA_spearman_correlation_genes"

csv_files = list.files(FGSEA_results, pattern = "\\.csv$", full.names = TRUE)
read_fgsea_nes = function(file_path) {
  df = read.csv(file_path)
  #if some fgsea results in no enriched pathways
  if (all(trimws(df$pathway) == "")) {
    return(NULL)
  }
  df_out = df %>%
    dplyr::select(pathway, NES) %>%
    dplyr::rename(!!basename(file_path) := NES) 
  return(df_out)
}

nes_list_raw = lapply(csv_files, read_fgsea_nes)
nes_list = nes_list_raw[!sapply(nes_list_raw, is.null)]
nes_matrix = reduce(nes_list, full_join, by = "pathway")

#format colnames to save 
colnames(nes_matrix)[-1] = str_remove(colnames(nes_matrix)[-1], "\\.csv$")
colnames(nes_matrix) = stringr::str_remove(colnames(nes_matrix), "^FGSEA_")
colnames(nes_matrix) = stringr::str_remove(colnames(nes_matrix), "_DEG_results")
write.csv(nes_matrix, file = "NES_matrix_all_fgsea.csv", quote = FALSE,row.names = FALSE)

#plotting heatmap for the commmon pathways with NES scores for each condition
nes_matrix = column_to_rownames(nes_matrix, var = "pathway")
nes_matrix = as.matrix(nes_matrix)
nes_matrix[is.na(nes_matrix)] = 0

library(pheatmap)
pdf("NES_heatmap.pdf", width = 30,height = 30)
pheatmap(nes_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_col = 10, 
         cellwidth = 25,       
         cellheight = 10,
         color = colorRampPalette(c("grey", "lightblue", "darkblue"))(100), 
         main = "Normalized Enrichment Score with BTM GSEA",height =40,width=10,
         legend = TRUE)
dev.off()
