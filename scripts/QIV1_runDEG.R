library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(readxl)
library(tibble)
library(sva)

#mapping to ensembl genes and target ids
protein_coding = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/protein_coding_ensemble_hgnclist.csv")

#========================================================================
# PREPROCESSING metadata, count matrix & serology data----
# =======================================================================
#from HI assay identify total responders and filter count matrix accordingly
#info available for 95 samples
QIV1_responder = read.table("/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv", sep=',',
                            header = TRUE)
QIV1_total_resp = QIV1_responder %>% dplyr::filter(TotResp4.0 == 4) #32
QIV1_double_resp = QIV1_responder %>% dplyr::filter(TotResp4.0 == 2) #20
QIV1_triple_resp = QIV1_responder %>% dplyr::filter(TotResp4.0 == 3)#38
QIV1_single_resp = QIV1_responder %>% dplyr::filter(TotResp4.0 == 1)#4
QIV1_non_resp = QIV1_responder %>% dplyr::filter(TotResp4.0 == 0)#1

#transcriptomics metadata info available for 100/99, V1/V2 samples

QIV1_meta =  read_excel("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/QIV1_India_MetaData_Transcriptomics.xlsx", sheet = 1)
QIV1_meta_filt = QIV1_meta[,c(2,3,6,7,30)]
colnames(QIV1_meta_filt)[1] = "SubjectID"

#create a new column combining visitID and subjectID
QIV1_meta_filt <- QIV1_meta_filt %>%
  mutate(SubjectIDNew = paste0(SubjectID, Visit))

#remove certain samples
#'QIV1KEM085' no data available for V2
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM085')
# 33 is a non-responder to all strains
#QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM033')
# 37 no serology data
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM037')
# no consent
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM019')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM055')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM064')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM077')

#samples with very little read counts for genes outliers for at least one visit
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM001')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM079')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM063')

#strain responders
HongKong_res = QIV1_responder %>%
  dplyr::filter(RES4.0.A.HongKong == 'TRUE')
Victoria_res = QIV1_responder %>%
  dplyr::filter(RES4.0.A.Victoria == 'TRUE')
Phuket_res = QIV1_responder %>%
  dplyr::filter(RES4.0.B.Phuket == 'TRUE')
Washington_res = QIV1_responder %>%
  dplyr::filter(RES4.0.B.Washington == 'TRUE')


#changing labels to HR and LR high and low responders
QIV1_meta_filt$Status = "LR"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_total_resp$SubjectID] = "HR"
# QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_triple_resp$SubjectID] = "HR"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_single_resp$SubjectID] = "NR"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_non_resp$SubjectID] = "NR"


#adding responder status for each strain in the metadata
QIV1_meta_filt$Responder_HongKong = "NonResponder"
QIV1_meta_filt$Responder_HongKong[QIV1_meta_filt$SubjectID %in% HongKong_res$SubjectID] = "Responder"
QIV1_meta_filt$Responder_HongKong = factor(QIV1_meta_filt$Responder_HongKong, levels = c("NonResponder", "Responder"))

QIV1_meta_filt$Responder_Victoria = "NonResponder"
QIV1_meta_filt$Responder_Victoria[QIV1_meta_filt$SubjectID %in% Victoria_res$SubjectID] = "Responder"
QIV1_meta_filt$Responder_Victoria = factor(QIV1_meta_filt$Responder_Victoria, levels = c("NonResponder", "Responder"))

QIV1_meta_filt$Responder_Phuket = "NonResponder"
QIV1_meta_filt$Responder_Phuket[QIV1_meta_filt$SubjectID %in% Phuket_res$SubjectID] = "Responder"
QIV1_meta_filt$Responder_Phuket = factor(QIV1_meta_filt$Responder_Phuket, levels = c("NonResponder", "Responder"))

QIV1_meta_filt$Responder_Washington = "NonResponder"
QIV1_meta_filt$Responder_Washington[QIV1_meta_filt$SubjectID %in% Washington_res$SubjectID] = "Responder"
QIV1_meta_filt$Responder_Washington = factor(QIV1_meta_filt$Responder_Washington, levels = c("NonResponder", "Responder"))

#=============================
#COUNT matrix preprocessing
#=============================
QIV1count = read.delim("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BA088-INCENTIVE_eukaryota_count_table.txt",
                       header = TRUE,check.names = FALSE, row.names = 1)

#split count matrix based on visit V2 or V3
#rename dfs to V1 and V2 from V2 and V3 to avoid confusion
QIV1_V1 <- QIV1count[, endsWith(colnames(QIV1count), "V2")]
QIV1_V2 <- QIV1count[, endsWith(colnames(QIV1count), "V3")]
#V1 and V2 merged count matrix
QIV1_all = cbind(QIV1_V1,QIV1_V2)
#remove V2 and V3 suffix from the colnames
colnames(QIV1_V1) <- gsub("V2$", "", colnames(QIV1_V1))
colnames(QIV1_V2) <- gsub("V3$", "", colnames(QIV1_V2))

#remove the Ensemble version details 
QIV1_all_mtx <-as.matrix(QIV1_all)
QIV1_all_mtx<- rownames_to_column(QIV1_all,var = 'target_id')
QIV1_all_mtx$target_id <- sub("\\..*", "", QIV1_all_mtx$target_id)
QIV1_all_mtx <- inner_join(QIV1_all_mtx,protein_coding,by='target_id')

numeric_data <- QIV1_all_mtx[, 2:(ncol(QIV1_all_mtx) - 1)] |> dplyr::select(where(is.numeric))
row_sums <- rowSums(numeric_data, na.rm = TRUE)
QIV1_all_mtx$rowsums <- row_sums

#for duplicate HGNC gene names keep only those with max counts
QIV1_all_mtx.filtered <- QIV1_all_mtx %>%
  dplyr::group_by(HGNC_symbol) %>%
  dplyr::filter(rowsums == max(rowsums)) %>%
  dplyr::slice(1) %>%
  ungroup()

QIV1_all_mtx.filtered  = column_to_rownames(QIV1_all_mtx.filtered, var = 'target_id')
QIV1_all_mtx.filtered.new = QIV1_all_mtx.filtered[,1:(ncol(QIV1_all_mtx.filtered) - 2)]

#change V2 to V1 and then V3 to V2
colnames_old <- colnames(QIV1_all_mtx.filtered.new)
colnames_new <- colnames_old |>
  sub("V2$", "V1", x = _) |>
  sub("V3$", "V2", x = _)

colnames(QIV1_all_mtx.filtered.new) <- colnames_new

#subset count matrix for only samples in metadata df
QIV1_all.num <-as.matrix(QIV1_all_mtx.filtered.new)
QIV1_all.num = QIV1_all.num[, (colnames(QIV1_all.num) %in% QIV1_meta_filt$SubjectIDNew)]
#onlyprotein-coding genes
QIV1_all.num = QIV1_all.num[rownames(QIV1_all.num) %in% protein_coding$target_id, , drop = FALSE]
#remove zero sum rows
row_sums<- rowSums(QIV1_all.num)
QIV1_all.num<- QIV1_all.num[which(row_sums != 0),]

#check columnorder of countmatrix and rows of metadata
identical(QIV1_meta_filt$SubjectIDNew,colnames(QIV1_all.num))
QIV1_all.num = QIV1_all.num[,QIV1_meta_filt$SubjectIDNew]                                              

#make factors for status and visit
QIV1_meta_filt$Status = factor(QIV1_meta_filt$Status, levels = c("HR","LR", "NR"))
QIV1_meta_filt$Visit = factor(QIV1_meta_filt$Visit, levels = c("V1", "V2"))
QIV1_meta_filt$Status_Visit <- paste0(QIV1_meta_filt$Status, "_", QIV1_meta_filt$Visit)

#covariates
AGE <-as.numeric(QIV1_meta_filt$AGE)
Library_Batch <- factor(QIV1_meta_filt$Library_Batch)
SEX <-factor(QIV1_meta_filt$SEX)

#again split QIV1_all.num based on V1 and V2

QIV1_V1.num = QIV1_all.num[, endsWith(colnames(QIV1_all.num), "V1")]
QIV1_V2.num = QIV1_all.num[, endsWith(colnames(QIV1_all.num), "V2")]


#=====================================================================
# DEG analysis for overall effect HR vs LR----
# ====================================================================
#design matrix for comparing HR vs LR overall
design1 <- model.matrix(~0 + Status + Visit + SEX + AGE + Library_Batch, data = QIV1_meta_filt)
contrast1 <- makeContrasts(HRvsLR = StatusHR-StatusLR, levels = design1)
contrast1 <- makeContrasts(HRvsNR = StatusHR-StatusNR, levels = design1)


DGEList <- DGEList(QIV1_all.num)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
v.DGEList.norm <- voom(DGEList.norm, design1, plot = TRUE, normalize.method = "quantile")

corfit = duplicateCorrelation(v.DGEList.norm ,design1,block = QIV1_meta_filt$SubjectID)
corfit$consensus.correlation
v.DGEList.norm = voom(DGEList.norm,design1,normalize.method = "quantile",
                      block = QIV1_meta_filt$SubjectID,correlation = corfit$consensus.correlation)
fit <- lmFit(v.DGEList.norm, design1,block = QIV1_meta_filt$SubjectID,correlation = corfit$consensus.correlation)
fit.contrasts = contrasts.fit(fit, contrasts = contrast1)
fit_ebayes <- eBayes(fit.contrasts, trend=TRUE)

HRvsLR = topTable(fit_ebayes,coef = 1,number = Inf)

#==========================================================
#design matrix for comparing between timepoints & status 
#==========================================================

design2 <- model.matrix(~0 + Status_Visit + SEX + AGE + Library_Batch, data = QIV1_meta_filt)
contrast2 <- makeContrasts(HRvsLR_V1 = Status_VisitHR_V1 - Status_VisitLR_V1,
  HRvsLR_V2 = Status_VisitHR_V2 - Status_VisitLR_V2,
  HR_V2vsV1 = Status_VisitHR_V2 - Status_VisitHR_V1,
  LR_V2vsV1 = Status_VisitLR_V2 - Status_VisitLR_V1,
  levels = design2)
print(identical(QIV1_meta_filt$SubjectIDNew,colnames(QIV1_all.num)))
v.DGEList.norm <- voom(DGEList.norm, design2, plot = FALSE, normalize.method = "quantile")

corfit = duplicateCorrelation(v.DGEList.norm ,design2,block = QIV1_meta_filt$SubjectID)
corfit$consensus.correlation
v.DGEList.norm = voom(DGEList.norm,design2,normalize.method = "quantile",
                      block = QIV1_meta_filt$SubjectID,
                      correlation = corfit$consensus.correlation)
fit2 <- lmFit(v.DGEList.norm, design2,block = QIV1_meta_filt$SubjectID,
             correlation = corfit$consensus.correlation)
fit2.contrasts = contrasts.fit(fit2, contrasts = contrast2)
fit2_ebayes <- eBayes(fit2.contrasts, trend=TRUE)

HRvsLR_V1 = topTable(fit2_ebayes,coef = 1,number = Inf)
HRvsLR_V2 = topTable(fit2_ebayes,coef = 2,number = Inf)
HR_V2vsV1 = topTable(fit2_ebayes,coef = 3,number = Inf)
LR_V2vsV1 = topTable(fit2_ebayes,coef = 4,number = Inf)


overallresults = list(HRvsLR=HRvsLR,HRvsLR_V1=HRvsLR_V1,HRvsLR_V2=HRvsLR_V2,HR_V2vsV1=HR_V2vsV1,LR_V2vsV1=LR_V2vsV1)

#=============================================================================
# DEG analysis strain wise (for each of the 4 strains in the vaccine)----
#=============================================================================


strains <- c("HongKong", "Victoria", "Phuket", "Washington")
results_list <- list()

for (strain in strains) {
  response_col <- paste0("Responder_", strain)
  interaction_col <- paste0("Status_Visit_", strain)
  
  QIV1_meta_filt[[response_col]] <- factor(QIV1_meta_filt[[response_col]], levels = c("NonResponder", "Responder"))
  QIV1_meta_filt$Visit <- factor(QIV1_meta_filt$Visit, levels = c("V1", "V2"))
  QIV1_meta_filt[[interaction_col]] <- interaction(QIV1_meta_filt[[response_col]], QIV1_meta_filt$Visit, drop = TRUE)
  
  #=====================================
  #Responder vs NonResponder both visits
  #=====================================
  #remove batch effect using Combat-seq
  QIV1_adjusted <- ComBat_seq(QIV1_all.num, batch=QIV1_meta_filt$Library_Batch, group=NULL)
  
  design1 <- model.matrix(~0 + QIV1_meta_filt[[response_col]] + SEX + AGE + Library_Batch + Visit, data = QIV1_meta_filt)
  colnames(design1)[1:2] <- c("NonResponder", "Responder")

  contrast1 <- makeContrasts(ResponderVsNonResponder = Responder - NonResponder,levels = design1)
  
  DGE <- DGEList(QIV1_adjusted)
  DGE <- calcNormFactors(DGE, method = "TMM")
  print(identical(QIV1_meta_filt$SubjectIDNew,colnames(QIV1_all.num)))
  v1 <- voom(DGE, design1, normalize.method = "quantile", plot = FALSE)
  corfit1 <- duplicateCorrelation(v1, design1, block = QIV1_meta_filt$SubjectID)
  v1 <- voom(DGE, design1, normalize.method = "quantile", block = QIV1_meta_filt$SubjectID, correlation = corfit1$consensus.correlation)
  
  fit1 <- lmFit(v1, design1, block = QIV1_meta_filt$SubjectID, correlation = corfit1$consensus.correlation)
  fit1 <- contrasts.fit(fit1, contrast1)
  fit1 <- eBayes(fit1, trend = TRUE)
  
  results_list[[paste0(strain, "_ResponderVsNonResponder")]] <- topTable(fit1, coef = 1, number = Inf)
  
  #===========================================
  #V2 vs V1 in Res & NonRes and Res separately
  #============================================
  design2 <- model.matrix(~0 + QIV1_meta_filt[[interaction_col]] + SEX + AGE + Library_Batch, data = QIV1_meta_filt)
  colnames(design2)[1:4] = c("NonResponder.V1", "Responder.V1", "NonResponder.V2","Responder.V2")

  contrast2 <- makeContrasts(
    V2vsV1_Responder = Responder.V2 - Responder.V1,
    V2vsV1_NonResponder = NonResponder.V2 - NonResponder.V1,
    ResponderV2vsNonResponderV2 = Responder.V2 - NonResponder.V2,
    ResponderV1vsNonResponderV1 = Responder.V1 - NonResponder.V1,
    levels = design2)
  
  v2 <- voom(DGE, design2, normalize.method = "quantile", plot = FALSE)
  corfit2 <- duplicateCorrelation(v2, design2, block = QIV1_meta_filt$SubjectID)
  v2 <- voom(DGE, design2, normalize.method = "quantile", block = QIV1_meta_filt$SubjectID, correlation = corfit2$consensus.correlation)
  
  fit2 <- lmFit(v2, design2, block = QIV1_meta_filt$SubjectID, correlation = corfit2$consensus.correlation)
  fit2 <- contrasts.fit(fit2, contrast2)
  fit2 <- eBayes(fit2, trend = TRUE)
  
  results_list[[paste0(strain, "_Responder_V2vsV1")]] <- topTable(fit2, coef = 1, number = Inf)
  results_list[[paste0(strain, "_NonResponder_V2vsV1")]] <- topTable(fit2, coef = 2, number = Inf)
  results_list[[paste0(strain, "_ResponderV2vsNonResponderV2")]] <- topTable(fit2, coef = 3, number = Inf)
  results_list[[paste0(strain, "_ResponderV1vsNonResponderV1")]] <- topTable(fit2, coef = 3, number = Inf)
}

#Save DEG results as csv files
for (name in names(results_list)) {
  deg_results <- results_list[[name]]
  deg_results$target_id <- rownames(deg_results)
  
  all_genes <- left_join(deg_results, protein_coding, by = "target_id")
  
  all_genes$significance <- ifelse(
    all_genes$P.Value < 0.05 & abs(all_genes$logFC) > 1,
    ifelse(all_genes$logFC > 1, "Upregulated", "Downregulated"),
    "Not significant")
  write.csv(all_genes, file = paste0(name, "_DEG_results.csv"), row.names = FALSE, quote = FALSE)
}


#======================================
#Volcano plots for the DEGs----
#======================================

for (name in names(results_list)) {
  deg_results <- results_list[[name]]
  deg_results$target_id <- rownames(deg_results)
  
  all_genes <- left_join(deg_results, protein_coding, by = "target_id")
  
  all_genes$significance <- ifelse(
    all_genes$P.Value < 0.05 & abs(all_genes$logFC) > 1,
    ifelse(all_genes$logFC > 1, "Upregulated", "Downregulated"),
    "Not significant")
  
  vplot <- ggplot(all_genes) +
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
