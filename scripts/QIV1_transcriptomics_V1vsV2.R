library(readxl)
library(dplyr)
library(tibble)
library(limma)
library(edgeR)

QIV1 = read.table("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/BA088-INCENTIVE_eukaryota_count_table.txt")

#from HI assay identify total responders and filter count matrix accordingly
QIV1_responder = read.table("/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv", sep=',',
                            header = TRUE)
QIV1_total_resp = QIV1_responder %>% filter(TotResp4.0 == 4) #32
QIV1_double_resp = QIV1_responder %>% filter(TotResp4.0 == 2) #20
QIV1_triple_resp = QIV1_responder %>% filter(TotResp4.0 == 3)#38
QIV1_single_resp = QIV1_responder %>% filter(TotResp4.0 == 1)#4

#transcriptomics metadata

QIV1_meta =  read_excel("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/QIV1_India_MetaData_Transcriptomics.xlsx", sheet = 1)
QIV1_meta_filt = QIV1_meta[,c(2,3,6,7,30)]
colnames(QIV1_meta_filt)[1] = "SubjectID"

#changing labels to HR and LR high and low responders
QIV1_meta_filt$Status <- "LowResponder"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_total_resp$SubjectID] <- "HighResponder"
QIV1_meta_filt$Status[QIV1_meta_filt$SubjectID %in% QIV1_triple_resp$SubjectID] <- "HighResponder"

QIV1_meta_filt <- QIV1_meta_filt[QIV1_meta_filt$SubjectID != "QIV1KEM085", ]

library(dplyr)
QIV1_meta_filt <- QIV1_meta_filt %>% 
  mutate(Visit = case_when(
    Visit == 'V2' ~ 'V3',
    Visit == 'V1' ~ 'V2',  
    TRUE ~ Visit           
  ))
QIV1_meta_filt <- QIV1_meta_filt %>%
  mutate(SampleID = paste0(SubjectID, Visit))

#Reorder metadata to match columns of the count matrix
QIV1_meta_filt <- QIV1_meta_filt[match(colnames(QIV1.num), QIV1_meta_filt$SampleID), ]



#mapping to ensembl genes and target ids
library(ensembldb)
library(EnsDb.Hsapiens.v86)
Gn <- genes(EnsDb.Hsapiens.v86, columns=c("gene_id", "gene_name"))
Gn <- as_tibble(Gn)
Gn <- dplyr::rename(Gn, target_id = gene_id)
Gn <- dplyr::select(Gn, "target_id", "gene_name")


QIV1_mtx <-as.matrix(QIV1)
QIV1_mtx<- rownames_to_column(QIV1,var = 'target_id')
QIV1_mtx$target_id <- sub("\\..*", "", QIV1_mtx$target_id)
QIV1_mtx <- inner_join(QIV1_mtx,Gn,by='target_id')

numeric_data <- QIV1_mtx[, 2:(ncol(QIV1_mtx) - 1)] |> dplyr::select(where(is.numeric))
row_sums <- rowSums(numeric_data, na.rm = TRUE)
QIV1_mtx$rowsums <- row_sums

QIV1_mtx.filtered <- QIV1_mtx %>%
  dplyr::group_by(gene_name) %>%
  dplyr::filter(rowsums == max(rowsums)) %>%
  dplyr::slice(1) %>%
  ungroup()

QIV1_mtx.filtered  = column_to_rownames(QIV1_mtx.filtered, var = 'target_id')
QIV1_mtx.filtered.new = QIV1_mtx.filtered[,1:(ncol(QIV1_mtx.filtered) - 2)]

QIV1.num <-as.matrix(QIV1_mtx.filtered.new)
row_sums<- rowSums(QIV1.num)
QIV1.num<- QIV1.num[which(row_sums != 0),]
QIV1.num = QIV1.num[, !(colnames(QIV1.num) %in% 'QIV1KEM085V2')]

QIV1.num = QIV1.num[rownames(QIV1.num) %in% gf_id, , drop = FALSE]
identical(QIV1_meta_filt$SampleID, colnames(QIV1.num))

responder <- factor(QIV1_meta_filt$Status, levels = c("LowResponder", "HighResponder"))
age <-as.numeric(QIV1_meta_filt$AGE)
Library_Batch <- factor(QIV1_meta_filt$Library_Batch)
sex <-factor(QIV1_meta_filt$SEX)

QIV1_meta_filt$Visit <- factor(QIV1_meta_filt$Visit, levels = c("V2", "V3"))
visit <- QIV1_meta_filt$Visit

covariates <- c("sex", "age", "Library_Batch", "visit")
covariate_formula <- paste(covariates, collapse = "+")
design_formula <- as.formula(paste("~ responder +", covariate_formula))
design <- model.matrix(design_formula, data = QIV1_meta_filt)

design_inter <- model.matrix(~ responder * visit + SEX + AGE + Library_Batch, data = QIV1_meta_filt)


DGEList <- DGEList(QIV1.num)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")

v.DGEList.norm <- voom(DGEList.norm, design, plot = FALSE,normalize.method = "quantile")
corfit = duplicateCorrelation(v.DGEList.norm ,design,block = QIV1_meta_filt$SubjectID)
corfit$consensus.correlation #0.07379232
v.DGEList.norm = voom(DGEList.norm,design,normalize.method = "quantile",
                      block = QIV1_meta_filt$SubjectID,
                      correlation = corfit$consensus.correlation)
fit <- lmFit(v.DGEList.norm, design,block = QIV1_meta_filt$SubjectID,
             correlation = corfit$consensus.correlation)
fit_ebayes <- eBayes(fit, trend=TRUE)
deg_visit<- topTable(fit_ebayes,adjust ="BH", coef="visitV3", number = Inf, sort.by = "logFC")
deg_responder <- topTable(fit_ebayes, coef = "responderHighResponder", adjust = "BH", number = Inf, sort.by = "logFC")

# contrast.matrix <- makeContrasts(
#   V3_vs_V2 = groupV3,
#   High_vs_Low = responderHighResponder,
#   levels = design
# )
# 
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2, trend = TRUE)
# 
# deg_visit <- topTable(fit2, coef = "V3_vs_V2", adjust = "BH", number = Inf)
# deg_resp  <- topTable(fit2, coef = "High_vs_Low", adjust = "BH", number = Inf)

deg_responder <- topTable(fit_ebayes, coef = "responderHighResponder:visitV3", number = Inf, adjust = "BH")
top_interaction$target_id = rownames(top_interaction)
all_genes = inner_join(top_interaction, Gn, by = "target_id")
output_file <- paste0("QIV1_LR_DEGs_V1vsV2proteincoding_limma", ".csv")
write.table(all_genes, output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# dge <- estimateDisp(DGEList.norm, design_inter)
# fit_2 <- glmQLFit(dge, design_inter)
# qlf <- glmQLFTest(fit_2, coef = 7)
# pvals <- qlf$table
# pvals$adj_pval <- stats::p.adjust(pvals$PValue, method = 'BH')
# pvals$target_id <- rownames(pvals)
# all_genes <- inner_join(pvals, Gn, by = "target_id")


