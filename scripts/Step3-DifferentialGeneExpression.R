library(limma)
library(sva)
library(edgeR)
library(variancePartition)
library(BiocParallel)
library(dplyr)
library(ggplot2)

# Enable parallel processing for dream
register(MulticoreParam(workers = 8))
protein_coding = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/protein_coding_ensemble_hgnclist.csv")

# =========================================================================
# FILTERING & NORMALIZATION 
# =========================================================================
dge_global <- DGEList(QIV1_all.num)

#  Filter based strictly on Visit to ensure a unified gene universe for ALL models
keep_global <- filterByExpr(dge_global, group = QIV1_meta_filt_cond$Visit)
QIV1_filter <- dge_global[keep_global, , keep.lib.sizes = FALSE]
DGEList.norm <- calcNormFactors(QIV1_filter, method = "TMM")
results_list <- list()

# =========================================================================
# STRAIN-WISE DEG ANALYSIS----
# =========================================================================

strains <- c("HongKong", "Victoria", "Phuket", "Washington")


for (strain in strains) {
  response_col <- strain
  interaction_col <- paste0("Status_Visit_", strain)
  
  QIV1_meta_filt_cond[[response_col]] <- factor(QIV1_meta_filt_cond[[response_col]], levels = c("LR", "HR", "MR"))
  QIV1_meta_filt_cond$Visit <- factor(QIV1_meta_filt_cond$Visit, levels = c("V1", "V2"))
  QIV1_meta_filt_cond[[interaction_col]] <- interaction(
    QIV1_meta_filt_cond[[response_col]], 
    QIV1_meta_filt_cond$Visit, 
    drop = TRUE, sep = "_")
  #=====================================
  #Responder vs NonResponder both visits
  #=====================================
  design <- as.formula(paste0("~ 0 + ", response_col, 
                              " + Visit + SEX + AGE + Library_Batch + (1 | SubjectID)"))
  dge          <- DGEList(QIV1_all.num)
  keep         <- edgeR::filterByExpr(dge, group = QIV1_meta_filt_cond[[response_col]])
  QIV1_filter  <- dge[keep, , keep.lib.sizes = FALSE]
  DGEList.norm <- calcNormFactors(QIV1_filter, method = "TMM")
  vobj         <- voomWithDreamWeights(DGEList.norm, design, QIV1_meta_filt_cond, plot = TRUE)
  
  contrast <- makeContrastsDream(
    design,
    QIV1_meta_filt_cond,
    contrasts = c(
      HRvsLR = paste0(response_col, "HR - ", response_col, "LR"),
      MRvsLR = paste0(response_col, "MR - ", response_col, "LR"),
      HRvsMR = paste0(response_col, "HR - ", response_col, "MR")
    )
  )
  
  fit_final <- dream(vobj, design, QIV1_meta_filt_cond, L = contrast)
  
  
  results_list[[paste0(strain, "_HRvsLR")]] <- topTable(fit_final, coef = "HRvsLR", number = Inf)
  results_list[[paste0(strain, "_MRvsLR")]] <- topTable(fit_final, coef = "MRvsLR", number = Inf)
  results_list[[paste0(strain, "_HRvsMR")]] <- topTable(fit_final, coef = "HRvsMR", number = Inf)
  
  
  
  #===========================================
  #V2 vs V1 in Res & NonRes and Res separately
  #============================================
  design2 <- as.formula(paste0("~ 0 + ", interaction_col, 
                               " + SEX + AGE + Library_Batch + (1 | SubjectID)"))
  dge          <- DGEList(QIV1_all.num)
  keep         <- edgeR::filterByExpr(dge, group = QIV1_meta_filt_cond[[interaction_col]])
  QIV1_filter  <- dge[keep, , keep.lib.sizes = FALSE]
  DGEList.norm <- calcNormFactors(QIV1_filter, method = "TMM")
  vobj         <- voomWithDreamWeights(DGEList.norm, design2, QIV1_meta_filt_cond, plot = TRUE)
  
  contrast2 <- makeContrastsDream(
    design2,
    QIV1_meta_filt_cond,
    contrasts = c(
      V2vsV1_HR  = paste0(interaction_col, "HR_V2 - ", interaction_col, "HR_V1"),
      V2vsV1_LR  = paste0(interaction_col, "LR_V2 - ", interaction_col, "LR_V1"),
      V2vsV1_MR  = paste0(interaction_col, "MR_V2 - ", interaction_col, "MR_V1"),
      HRV2vsLRV2 = paste0(interaction_col, "HR_V2 - ", interaction_col, "LR_V2"),
      HRV1vsLRV1 = paste0(interaction_col, "HR_V1 - ", interaction_col, "LR_V1"),
      HRV2vsMRV2 = paste0(interaction_col, "HR_V2 - ", interaction_col, "MR_V2"),
      HRV1vsMRV1 = paste0(interaction_col, "HR_V1 - ", interaction_col, "MR_V1"),
      MRV2vsLRV2 = paste0(interaction_col, "MR_V2 - ", interaction_col, "LR_V2"),
      MRV1vsLRV1 = paste0(interaction_col, "MR_V1 - ", interaction_col, "LR_V1")
    )
  )
  
  fit_final2 <- dream(vobj, design2, QIV1_meta_filt_cond, L = contrast2)
  
  results_list[[paste0(strain, "_V2vsV1_HR")]]   <- topTable(fit_final2, coef = "V2vsV1_HR",  number = Inf)
  results_list[[paste0(strain, "_V2vsV1_LR")]]   <- topTable(fit_final2, coef = "V2vsV1_LR",  number = Inf)
  results_list[[paste0(strain, "_V2vsV1_MR")]]   <- topTable(fit_final2, coef = "V2vsV1_MR",  number = Inf)
  results_list[[paste0(strain, "_HRV2vsLRV2")]]  <- topTable(fit_final2, coef = "HRV2vsLRV2", number = Inf)
  results_list[[paste0(strain, "_HRV1vsLRV1")]]  <- topTable(fit_final2, coef = "HRV1vsLRV1", number = Inf)
  results_list[[paste0(strain, "_HRV2vsMRV2")]]  <- topTable(fit_final2, coef = "HRV2vsMRV2", number = Inf)
  results_list[[paste0(strain, "_HRV1vsMRV1")]]  <- topTable(fit_final2, coef = "HRV1vsMRV1", number = Inf)
  results_list[[paste0(strain, "_MRV2vsLRV2")]]  <- topTable(fit_final2, coef = "MRV2vsLRV2", number = Inf)
  results_list[[paste0(strain, "_MRV1vsLRV1")]]  <- topTable(fit_final2, coef = "MRV1vsLRV1", number = Inf)
}


# =========================================================================
# OVERALL AGGREGATED RESPONDER ANALYSIS----
# =========================================================================

interaction_col <- paste0("Status_Visit_Agg")

QIV1_meta_filt_cond$responder_group <- factor(QIV1_meta_filt_cond$responder_group, levels = c("LR", "HR", "MR"))
QIV1_meta_filt_cond$Visit <- factor(QIV1_meta_filt_cond$Visit, levels = c("V1", "V2"))
QIV1_meta_filt_cond[[interaction_col]] <- interaction(
  QIV1_meta_filt_cond$responder_group, 
  QIV1_meta_filt_cond$Visit, 
  drop = TRUE, sep = "_")


#===========================================
#V2 vs V1 in Res & NonRes and Res separately
#============================================
design_agg <- as.formula(paste0("~ 0 + ", interaction_col, 
                             " + SEX + AGE + Library_Batch + (1 | SubjectID)"))
#dge          <- DGEList(QIV1_all.num)
#keep         <- edgeR::filterByExpr(dge, group =  QIV1_meta_filt_cond$responder_group)
#QIV1_filter  <- dge[keep, , keep.lib.sizes = FALSE]
#DGEList.norm <- calcNormFactors(QIV1_filter, method = "TMM")
vobj         <- voomWithDreamWeights(DGEList.norm, design_agg, QIV1_meta_filt_cond, plot = TRUE)

contrast_agg <- makeContrastsDream(
  design_agg,
  QIV1_meta_filt_cond,
  contrasts = c(
    V2vsV1_HR  = paste0(interaction_col, "HR_V2 - ", interaction_col, "HR_V1"),
    V2vsV1_LR  = paste0(interaction_col, "LR_V2 - ", interaction_col, "LR_V1"),
    V2vsV1_MR  = paste0(interaction_col, "MR_V2 - ", interaction_col, "MR_V1"),
    HRV2vsLRV2 = paste0(interaction_col, "HR_V2 - ", interaction_col, "LR_V2"),
    HRV1vsLRV1 = paste0(interaction_col, "HR_V1 - ", interaction_col, "LR_V1"),
    HRV2vsMRV2 = paste0(interaction_col, "HR_V2 - ", interaction_col, "MR_V2"),
    HRV1vsMRV1 = paste0(interaction_col, "HR_V1 - ", interaction_col, "MR_V1"),
    MRV2vsLRV2 = paste0(interaction_col, "MR_V2 - ", interaction_col, "LR_V2"),
    MRV1vsLRV1 = paste0(interaction_col, "MR_V1 - ", interaction_col, "LR_V1")
  )
)

fit_final_agg <- dream(vobj, design_agg, QIV1_meta_filt_cond, L = contrast_agg)

results_list[[paste0("V2vsV1_HR")]]   <- topTable(fit_final_agg, coef = "V2vsV1_HR",  number = Inf)
results_list[[paste0("V2vsV1_LR")]]   <- topTable(fit_final_agg, coef = "V2vsV1_LR",  number = Inf)
results_list[[paste0("V2vsV1_MR")]]   <- topTable(fit_final_agg, coef = "V2vsV1_MR",  number = Inf)
results_list[[paste0("HRV2vsLRV2")]]  <- topTable(fit_final_agg, coef = "HRV2vsLRV2", number = Inf)
results_list[[paste0("HRV1vsLRV1")]]  <- topTable(fit_final_agg, coef = "HRV1vsLRV1", number = Inf)
results_list[[paste0("HRV2vsMRV2")]]  <- topTable(fit_final_agg, coef = "HRV2vsMRV2", number = Inf)
results_list[[paste0("HRV1vsMRV1")]]  <- topTable(fit_final_agg, coef = "HRV1vsMRV1", number = Inf)
results_list[[paste0("MRV2vsLRV2")]]  <- topTable(fit_final_agg, coef = "MRV2vsLRV2", number = Inf)
results_list[[paste0("MRV1vsLRV1")]]  <- topTable(fit_final_agg, coef = "MRV1vsLRV1", number = Inf)



# =========================================================================
# SAVE DEG RESULTS ----
# =========================================================================
for (name in names(results_list)) {
  deg_results <- results_list[[name]]
  
  deg_results$target_id <- rownames(deg_results)
  deg_results = inner_join(deg_results,protein_coding,by = "target_id")
  deg_results = deg_results %>%
    arrange(desc(abs(logFC))) %>%   
    distinct(HGNC_symbol, .keep_all = TRUE)
  
  deg_results$significance <- ifelse(
    deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1,
    ifelse(deg_results$logFC > 1, "Upregulated", "Downregulated"),
    "Not significant")
  write.csv(deg_results, file = paste0(name, "_DEG_results.csv"), row.names = FALSE, quote = FALSE)
}


#======================================
# VOLCANO PLOTS ----
#======================================

for (name in names(results_list)) {
  deg_results <- results_list[[name]]
  deg_results$HGNC_symbol <- rownames(deg_results)
  
  deg_results$significance <- ifelse(
    deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1,
    ifelse(deg_results$logFC > 1, "Upregulated", "Downregulated"),
    "Not significant")
  
  vplot <- ggplot(deg_results) +
    aes(y = -log10(P.Value), x = logFC, color = significance) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "Upregulated" = "#A020F0",
      "Downregulated" = "#A0522D",
      "Not significant" = "grey60")) +
    geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour = "black", size = 1) +
    geom_vline(xintercept = 1, linetype = "longdash", colour = "black", size = 1) +
    geom_vline(xintercept = -1, linetype = "longdash", colour = "black", size = 1) +
    labs(title = name) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12))
  
  ggsave(paste0(name, ".png"), plot = vplot, width = 15, height = 15)
}

