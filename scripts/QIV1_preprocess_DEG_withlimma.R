library(ggplot2)
library(readxl)
library(dplyr)
library(tibble)
library(limma)
library(edgeR)
library(readxl)


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

#changing labels to HR and LR high and low responders
QIV1_meta_filt$Status = "LowResponder"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_total_resp$SubjectID] = "HighResponder"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_triple_resp$SubjectID] = "HighResponder"


#create a new column combining visitID and subjectID
QIV1_meta_filt <- QIV1_meta_filt %>%
  mutate(SubjectIDNew = paste0(SubjectID, Visit))

#remove certain samples
#'QIV1KEM085' no data available for V2
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM085')
# 33 is a non-responder to all strains
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM033')
# 37 no serology data
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM037')
# no consent
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM019')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM055')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM064')
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(!SubjectID == 'QIV1KEM077')

#get separate metadata for HR and LR
QIV1_meta_HR = QIV1_meta_filt %>% dplyr::filter(Status == 'HighResponder')
QIV1_meta_LR = QIV1_meta_filt %>% dplyr::filter(Status == 'LowResponder')

###analysis for day0 HRvsLR, day3 HRvsLR----

QIV1count = read.table("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BA088-INCENTIVE_eukaryota_count_table.txt")
QIV1tpm = read.table("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BA088-INCENTIVE_eukaryota_TPM_table.txt")
#split count matrix based on visit V2 or V3
#rename dfs to V1 and V2 from V2 and V3
QIV1_V1 <- QIV1count[, endsWith(colnames(QIV1count), "V2")]
QIV1_V2 <- QIV1count[, endsWith(colnames(QIV1count), "V3")]

#V1 and V2 merged count matrix
QIV1_all = cbind(QIV1_V1,QIV1_V2)

#remove V2 and V3 suffix from the colnames
colnames(QIV1_V1) <- gsub("V2$", "", colnames(QIV1_V1))
colnames(QIV1_V2) <- gsub("V3$", "", colnames(QIV1_V2))

QIV1_meta_filt_V2 = QIV1_meta_filt %>% dplyr::filter(Visit == 'V2')
QIV1_meta_filt_V1 = QIV1_meta_filt %>% dplyr::filter(Visit == 'V1')

#removing version ID for each gene and keeping rows with
#max rowsums and later removing rows with zero sums

##############for QIV1_V1##################
QIV1_V1_mtx <-as.matrix(QIV1_V1)
QIV1_V1_mtx<- rownames_to_column(QIV1_V1,var = 'target_id')
QIV1_V1_mtx$target_id <- sub("\\..*", "", QIV1_V1_mtx$target_id)
QIV1_V1_mtx <- inner_join(QIV1_V1_mtx,Gn,by='target_id')

numeric_data <- QIV1_V1_mtx[, 2:(ncol(QIV1_V1_mtx) - 1)] |> dplyr::select(where(is.numeric))
row_sums <- rowSums(numeric_data, na.rm = TRUE)
QIV1_V1_mtx$rowsums <- row_sums

QIV1_V1_mtx.filtered <- QIV1_V1_mtx %>%
  dplyr::group_by(gene_name) %>%
  dplyr::filter(rowsums == max(rowsums)) %>%
  dplyr::slice(1) %>%
  ungroup()

QIV1_V1_mtx.filtered  = column_to_rownames(QIV1_V1_mtx.filtered, var = 'target_id')
QIV1_V1_mtx.filtered.new = QIV1_V1_mtx.filtered[,1:(ncol(QIV1_V1_mtx.filtered) - 2)]

QIV1_V1.num <-as.matrix(QIV1_V1_mtx.filtered.new)

QIV1_V1.num = QIV1_V1.num[, (colnames(QIV1_V1.num) %in% QIV1_meta_filt_V1$SubjectID)]
QIV1_V1.num = QIV1_V1.num[rownames(QIV1_V1.num) %in% gf_id, , drop = FALSE]
row_sums<- rowSums(QIV1_V1.num)
QIV1_V1.num<- QIV1_V1.num[which(row_sums != 0),]

identical(QIV1_meta_filt_V1$SubjectID,colnames(QIV1_V1.num))
QIV1_V1.num = QIV1_V1.num[,QIV1_meta_filt_V1$SubjectID]

#covariates for DEG
group <- factor(QIV1_meta_filt_V1$Status, levels = c("LowResponder", "HighResponder"))
AGE <-as.numeric(QIV1_meta_filt_V1$AGE)
Library_Batch <- factor(QIV1_meta_filt_V1$Library_Batch)
SEX <-factor(QIV1_meta_filt_V1$SEX)
covariates <- c("SEX", "AGE", "Library_Batch")
covariate_formula <- paste(covariates, collapse = "+")
design_formula <- as.formula(paste("~ group +", covariate_formula))
design <- model.matrix(design_formula, data = QIV1_meta_filt_V1)


#DEG with limma###
DGEList <- DGEList(QIV1_V1.num)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
v.DEGList.norm <- voom(DGEList.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.norm, design)
fit_ebayes <- eBayes(fit, trend=TRUE)
deg_results<- topTable(fit_ebayes,adjust ="BH", coef="groupHighResponder", number = Inf, sort.by = "logFC")
deg_results$target_id = rownames(deg_results)
all_genes = inner_join(deg_results, Gn, by = "target_id")


#volcano plots ----
all_genes$significance <- ifelse(all_genes$P.Value < 0.05 & abs(all_genes$logFC) > 1,
                                 ifelse(all_genes$logFC > 1, "Upregulated", "Downregulated"),
                                 "Not significant")
output_file <- paste0("QIV1_DEGs_HRvsLR_Day0proteincoding_limma", ".csv")
write.table(all_genes, output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


vplot <- ggplot(all_genes) +
  aes(y=-log10(P.Value), x=logFC, text = paste("Symbol:", gene_name), color = significance) +
  geom_point(size=3) +
  scale_color_manual(values = c("Upregulated" = "#A020F0",
                                "Downregulated" = "#A0522D",
                                "Not significant" = "grey60")) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="black", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="black", size=1) +
  labs(title="QIV1_DEGs_HRvsLR_Day0proteincoding_limma ",
       subtitle = "QIV1-HRvsLR-V1-BeforeVaccination") +
  theme_bw()
ggsave('QIV1_DEGs_HRvsLR_Day0proteincoding_limma.png',plot = vplot,width = 15, height=15)


##############for QIV1_V2##################
QIV1_V2_mtx <-as.matrix(QIV1_V2)
QIV1_V2_mtx<- rownames_to_column(QIV1_V2,var = 'target_id')
QIV1_V2_mtx$target_id <- sub("\\..*", "", QIV1_V2_mtx$target_id)
QIV1_V2_mtx <- inner_join(QIV1_V2_mtx,Gn,by='target_id')

numeric_data <- QIV1_V2_mtx[, 2:(ncol(QIV1_V2_mtx) - 1)] |> dplyr::select(where(is.numeric))
row_sums <- rowSums(numeric_data, na.rm = TRUE)
QIV1_V2_mtx$rowsums <- row_sums

QIV1_V2_mtx.filtered <- QIV1_V2_mtx %>%
  dplyr::group_by(gene_name) %>%
  dplyr::filter(rowsums == max(rowsums)) %>%
  dplyr::slice(1) %>%
  ungroup()

QIV1_V2_mtx.filtered  = column_to_rownames(QIV1_V2_mtx.filtered, var = 'target_id')
QIV1_V2_mtx.filtered.new = QIV1_V2_mtx.filtered[,1:(ncol(QIV1_V2_mtx.filtered) - 2)]

QIV1_V2.num <-as.matrix(QIV1_V2_mtx.filtered.new)

QIV1_V2.num = QIV1_V2.num[, (colnames(QIV1_V2.num) %in% QIV1_meta_filt_V2$SubjectID)]
QIV1_V2.num = QIV1_V2.num[rownames(QIV1_V2.num) %in% gf_id, , drop = FALSE]
row_sums<- rowSums(QIV1_V2.num)
QIV1_V2.num<- QIV1_V2.num[which(row_sums != 0),]

identical(QIV1_meta_filt_V2$SubjectID,colnames(QIV1_V2.num))
QIV1_V2.num = QIV1_V2.num[,QIV1_meta_filt_V2$SubjectID]

#covariates for DEG
group <- factor(QIV1_meta_filt_V2$Status, levels = c("LowResponder", "HighResponder"))
AGE <-as.numeric(QIV1_meta_filt_V2$AGE)
Library_Batch <- factor(QIV1_meta_filt_V2$Library_Batch)
SEX <-factor(QIV1_meta_filt_V2$SEX)
covariates <- c("SEX", "AGE", "Library_Batch")
covariate_formula <- paste(covariates, collapse = "+")
design_formula <- as.formula(paste("~ group +", covariate_formula))
design <- model.matrix(design_formula, data = QIV1_meta_filt_V2)


#DEG with limma###
DGEList <- DGEList(QIV1_V2.num)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
v.DEGList.norm <- voom(DGEList.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.norm, design)
fit_ebayes <- eBayes(fit, trend=TRUE)
deg_results<- topTable(fit_ebayes,adjust ="BH", coef="groupHighResponder", number = Inf, sort.by = "logFC")
deg_results$target_id = rownames(deg_results)
all_genes = inner_join(deg_results, Gn, by = "target_id")


#volcano plots ----
all_genes$significance <- ifelse(all_genes$P.Value < 0.05 & abs(all_genes$logFC) > 1,
                                 ifelse(all_genes$logFC > 1, "Upregulated", "Downregulated"),
                                 "Not significant")
output_file <- paste0("QIV1_DEGs_HRvsLR_Day3proteincoding_limma", ".csv")
write.table(all_genes, output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


vplot <- ggplot(all_genes) +
  aes(y=-log10(P.Value), x=logFC, text = paste("Symbol:", gene_name), color = significance) +
  geom_point(size=3) +
  scale_color_manual(values = c("Upregulated" = "#A020F0",
                                "Downregulated" = "#A0522D",
                                "Not significant" = "grey60")) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="black", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="black", size=1) +
  labs(title="QIV1_DEGs_HRvsLR_Day3proteincoding_limma ",
       subtitle = "QIV1-HRvsLR-V2-AfterVaccination") +
  theme_bw()
ggsave('QIV1_DEGs_HRvsLR_Day3proteincoding_limma.png',plot = vplot,width = 15, height=15)



###analysis for HR (day0 vs day3) or LR(day0 vs day3)----
QIV1_all_mtx <-as.matrix(QIV1_all)
QIV1_all_mtx<- rownames_to_column(QIV1_all,var = 'target_id')
QIV1_all_mtx$target_id <- sub("\\..*", "", QIV1_all_mtx$target_id)
QIV1_all_mtx <- inner_join(QIV1_all_mtx,Gn,by='target_id')

numeric_data <- QIV1_all_mtx[, 2:(ncol(QIV1_all_mtx) - 1)] |> dplyr::select(where(is.numeric))
row_sums <- rowSums(numeric_data, na.rm = TRUE)
QIV1_all_mtx$rowsums <- row_sums

QIV1_all_mtx.filtered <- QIV1_all_mtx %>%
  dplyr::group_by(gene_name) %>%
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


QIV1_all.num <-as.matrix(QIV1_all_mtx.filtered.new)

QIV1_HR.num = QIV1_all.num[, (colnames(QIV1_all.num) %in% QIV1_meta_HR$SubjectIDNew)]
QIV1_LR.num = QIV1_all.num[, (colnames(QIV1_all.num) %in% QIV1_meta_LR$SubjectIDNew)]

identical(QIV1_meta_LR$SubjectIDNew,colnames(QIV1_LR.num))
QIV1_LR.num = QIV1_LR.num[,QIV1_meta_LR$SubjectIDNew]                                              



QIV1_meta_LR$Visit <- factor(QIV1_meta_LR$Visit, levels = c("V1", "V2"))
AGE <-as.numeric(QIV1_meta_LR$AGE)
Library_Batch <- factor(QIV1_meta_LR$Library_Batch)
SEX <-factor(QIV1_meta_LR$SEX)
covariates <- c("SEX", "AGE", "Library_Batch")
covariate_formula <- paste(covariates, collapse = "+")
design_formula <- as.formula(paste("~ Visit +", covariate_formula))
design <- model.matrix(design_formula, data = QIV1_meta_LR)

QIV1_LR.num = QIV1_LR.num[rownames(QIV1_LR.num) %in% gf_id, , drop = FALSE]
myDGEList <- DGEList(QIV1_LR.num)
myDGEList.norm <- calcNormFactors(myDGEList, method = "TMM")
v.DEGList.norm <- voom(myDGEList.norm, design, plot = FALSE, normalize.method = "quantile")

corfit = duplicateCorrelation(v.DGEList.norm ,design,block = QIV1_meta_filt$SubjectID)
corfit$consensus.correlation
v.DGEList.norm = voom(DGEList.norm,design,normalize.method = "quantile",
                      block = QIV1_meta_filt$SubjectID,
                      correlation = corfit$consensus.correlation)
fit <- lmFit(v.DGEList.norm, design,block = QIV1_meta_filt$SubjectID,
             correlation = corfit$consensus.correlation)

fit_ebayes <- eBayes(fit, trend=TRUE)
deg_results<- topTable(fit_ebayes,adjust ="BH", coef="VisitV2", number = Inf, sort.by = "logFC")
deg_results$target_id = rownames(deg_results)
all_genes = inner_join(deg_results, Gn, by = "target_id")
