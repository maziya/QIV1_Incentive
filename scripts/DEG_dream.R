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
QIV1_meta_filt = QIV1_meta_filt %>%
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
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM030')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM040')

#high percentage of rRNA
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM023')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM014')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM043')


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
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_triple_resp$SubjectID] = "HR"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_single_resp$SubjectID] = "NR"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_non_resp$SubjectID] = "NR"


#adding responder status for each strain in the metadata
QIV1_meta_filt$Responder_HongKong = "LR"
QIV1_meta_filt$Responder_HongKong[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$HongKong == "HR"]] = "HR"
QIV1_meta_filt$Responder_HongKong[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$HongKong == "MR"]] = "MR"
QIV1_meta_filt$Responder_HongKong = factor(QIV1_meta_filt$Responder_HongKong, levels = c("LR", "HR", "MR"))

QIV1_meta_filt$Responder_Victoria = "LR"
QIV1_meta_filt$Responder_Victoria[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Victoria == "HR"]] = "HR"
QIV1_meta_filt$Responder_Victoria[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Victoria == "MR"]] = "MR"
QIV1_meta_filt$Responder_Victoria = factor(QIV1_meta_filt$Responder_Victoria, levels = c("LR", "HR", "MR"))

QIV1_meta_filt$Responder_Phuket = "LR"
QIV1_meta_filt$Responder_Phuket[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Phuket == "HR"]] = "HR"
QIV1_meta_filt$Responder_Phuket[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Phuket == "MR"]] = "MR"
QIV1_meta_filt$Responder_Phuket = factor(QIV1_meta_filt$Responder_Phuket, levels = c("LR", "HR", "MR"))

QIV1_meta_filt$Responder_Washington = "LR"
QIV1_meta_filt$Responder_Washington[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Washington == "HR"]] = "HR"
QIV1_meta_filt$Responder_Washington[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Washington == "MR"]] = "MR"
QIV1_meta_filt$Responder_Washington = factor(QIV1_meta_filt$Responder_Washington, levels = c("LR", "HR", "MR"))

#=============================
#COUNT matrix preprocessing
#=============================
QIV1count = read.delim("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BA088-INCENTIVE_eukaryota_count_table.txt",
                       header = TRUE,check.names = FALSE, row.names = 1)

#split count matrix based on visit V2 or V3
#rename dfs to V1 and V2 from V2 and V3 to avoid confusion
QIV1_V1 = QIV1count[, endsWith(colnames(QIV1count), "V2")]
QIV1_V2 = QIV1count[, endsWith(colnames(QIV1count), "V3")]
#V1 and V2 merged count matrix
QIV1_all = cbind(QIV1_V1,QIV1_V2)
#remove V2 and V3 suffix from the colnames
colnames(QIV1_V1) = gsub("V2$", "", colnames(QIV1_V1))
colnames(QIV1_V2) = gsub("V3$", "", colnames(QIV1_V2))

#remove the Ensemble version details 
QIV1_all_mtx =as.matrix(QIV1_all)
QIV1_all_mtx= rownames_to_column(QIV1_all,var = 'target_id')
QIV1_all_mtx$target_id = sub("\\..*", "", QIV1_all_mtx$target_id)
QIV1_all_mtx = inner_join(QIV1_all_mtx,protein_coding,by='target_id')

numeric_data = QIV1_all_mtx[, 2:(ncol(QIV1_all_mtx) - 1)] %>% dplyr::select(where(is.numeric))
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
colnames_new <- colnames_old %>%
  sub("V2$", "V1", x = _) %>%
  sub("V3$", "V2", x = _)

colnames(QIV1_all_mtx.filtered.new) <- colnames_new

#subset count matrix for only samples in metadata df
QIV1_all.num <-as.matrix(QIV1_all_mtx.filtered.new)
QIV1_all.num = QIV1_all.num[, (colnames(QIV1_all.num) %in% QIV1_meta_filt_cond$SubjectIDNew)]
#onlyprotein-coding genes
QIV1_all.num = QIV1_all.num[rownames(QIV1_all.num) %in% protein_coding$target_id, , drop = FALSE]
#remove zero sum rows
row_sums<- rowSums(QIV1_all.num)
QIV1_all.num<- QIV1_all.num[which(row_sums != 0),]

#check columnorder of countmatrix and rows of metadata
identical(QIV1_meta_filt_cond$SubjectIDNew,colnames(QIV1_all.num))
QIV1_all.num = QIV1_all.num[,QIV1_meta_filt_cond$SubjectIDNew]                                              

#split QIV1_all.num based on V1 and V2

QIV1_V1.num = QIV1_all.num[, endsWith(colnames(QIV1_all.num), "V1")]
QIV1_V2.num = QIV1_all.num[, endsWith(colnames(QIV1_all.num), "V2")]

#make factors for status and visit
QIV1_meta_filt$Status = factor(QIV1_meta_filt$Status, levels = c("HR","MR", "LR"))
QIV1_meta_filt$Visit = factor(QIV1_meta_filt$Visit, levels = c("V1", "V2"))
QIV1_meta_filt$Status_Visit <- paste0(QIV1_meta_filt$Status, "_", QIV1_meta_filt$Visit)

#covariates
AGE <-as.numeric(QIV1_meta_filt_cond$AGE)
Library_Batch <- factor(QIV1_meta_filt_cond$Library_Batch)
SEX <-factor(QIV1_meta_filt_cond$SEX)



DGEList = DGEList(QIV1_all.num)
keep = edgeR::filterByExpr(DGEList, design= design)
QIV1_filter = DGEList[keep,,keep.lib.sizes=FALSE]
DGEList.norm = calcNormFactors(QIV1_filter, method = "TMM")

QIV1_meta_filt_cond$Visit =relevel(QIV1_meta_filt_cond$Visit, ref = "V1")

#Dream pipeline for DEGs
# form = ~ 0 + Visit + SEX + AGE + Library_Batch + (1 | SubjectID)
# vobj = voomWithDreamWeights(DGEList.norm,form,QIV1_meta_filt_cond, plot = TRUE)
# 
# fit_temp = dream(vobj, form, QIV1_meta_filt_cond)
# coef_names = colnames(coef(fit_temp))
# contrast = makeContrasts(
#   Visit2vsVisit1 = VisitV2 - VisitV1,
#   levels = coef_names)
#  
# fit_final = dream(vobj,form,QIV1_meta_filt_cond,L = contrast)
# fit_final = eBayes(fit_final)
# 
# res_V2vsV1 = topTable(fit_final, coef = "Visit2vsVisit1",number = Inf,sort.by = "P", adjust.method = "BH")

strains <- c("HongKong", "Victoria", "Phuket", "Washington")
results_list <- list()

for (strain in strains) {
  response_col <- paste0("Responder_", strain)
  interaction_col <- paste0("Status_Visit_", strain)
  
  QIV1_meta_filt[[response_col]] <- factor(QIV1_meta_filt[[response_col]], levels = c("LR", "HR", "MR"))
  QIV1_meta_filt$Visit <- factor(QIV1_meta_filt$Visit, levels = c("V1", "V2"))
  QIV1_meta_filt[[interaction_col]] <- interaction(QIV1_meta_filt[[response_col]], QIV1_meta_filt$Visit, drop = TRUE)
  
  #=====================================
  #Responder vs NonResponder both visits
  #=====================================
  
  design =  ~ 0 + Visit + SEX + AGE + Library_Batch + (1 | SubjectID)
  vobj = voomWithDreamWeights(DGEList.norm,design1,QIV1_meta_filt_cond, plot = TRUE)
  DGEList = DGEList(QIV1_all.num)
  keep = edgeR::filterByExpr(DGEList, design= design)
  QIV1_filter = DGEList[keep,,keep.lib.sizes=FALSE]
  DGEList.norm = calcNormFactors(QIV1_filter, method = "TMM")

  fit_temp = dream(vobj, design1, QIV1_meta_filt_cond)
  coef_names = colnames(coef(fit_temp))
  contrast <- makeContrasts(HRvsLR = HR-LR, MRvsLR = MR-LR, HRvsMR = HR-MR, levels = coef_names)
  
  fit_final = dream(vobj,design,QIV1_meta_filt_cond,L = contrast)
  fit_final = eBayes(fit_final)
  
  results_list[[paste0(strain, "_HRvsLR")]] <- topTable(fit_final, coef = 1, number = Inf)
  results_list[[paste0(strain, "_MRvsLR")]] <- topTable(fit_final, coef = 2, number = Inf)
  results_list[[paste0(strain, "_HRvsMR")]] <- topTable(fit_final, coef = 3, number = Inf)
  
  
  #===========================================
  #V2 vs V1 in Res & NonRes and Res separately
  #============================================
  design2 = 0 + Visit + QIV1_meta_filt[[interaction_col]] +SEX + AGE + Library_Batch + (1 | SubjectID)
  colnames(design2)[1:6] = c("LR.V1","HR.V1", "MR.V1", "LR.V2","HR.V2","MR.V2")
  
  contrast2 <- makeContrasts(
    V2vsV1_HR = HR.V2 - HR.V1,
    V2vsV1_LR = LR.V2 - LR.V1,
    V2vsV1_MR = MR.V2 - MR.V1,
    HRV2vsLRV2 = HR.V2 - LR.V2,
    HRV1vsLRV1 = HR.V1 - LR.V1,
    HRV2vsMRV2 = HR.V2 - MR.V2,
    HRV1vsMRV1 = HR.V1 - MR.V1,
    MRV2vsLRV2 = MR.V2 - LR.V2,
    MRV1vsLRV1 = MR.V1 - LR.V1,
    levels = design2)
  
  DGE = DGEList(QIV1_all.num)
  keep = edgeR::filterByExpr(DGE, design= design2)
  DGE = DGE[keep,,keep.lib.sizes=FALSE]
  DGE = calcNormFactors(DGE, method = "TMM")
  print(identical(QIV1_meta_filt$SubjectIDNew,colnames(QIV1_all.num)))
 
  fit_temp = dream(vobj, design2, QIV1_meta_filt_cond)
  coef_names = colnames(coef(fit_temp))
  
  fit_final2 = dream(vobj,design2,QIV1_meta_filt_cond,L = contrast2)
  fit_final2 = eBayes(fit_final2)
  
  results_list[[paste0(strain, "_V2vsV1_HR")]] <- topTable(fit_final2, coef = 1, number = Inf)
  results_list[[paste0(strain, "_V2vsV1_LR")]] <- topTable(fit_final2, coef = 2, number = Inf)
  results_list[[paste0(strain, "_V2vsV1_MR")]] <- topTable(fit_final2, coef = 3, number = Inf)
  results_list[[paste0(strain, "_HRV2vsLRV2")]] <- topTable(fit_final2, coef = 4, number = Inf)
  results_list[[paste0(strain, "_HRV1vsLRV1")]] <- topTable(fit_final2, coef = 5, number = Inf)
  results_list[[paste0(strain, "_HRV2vsMRV2")]] <- topTable(fit_final2, coef = 6, number = Inf)
  results_list[[paste0(strain, "_HRV1vsMRV1")]] <- topTable(fit_final2, coef = 7, number = Inf)
  results_list[[paste0(strain, "_MRV2vsLRV2")]] <- topTable(fit_final2, coef = 8, number = Inf)
  results_list[[paste0(strain, "_MRV1vsLRV1")]] <- topTable(fit_final2, coef = 9, number = Inf)
}
