
library(edgeR)
library(GSVA)
library(limma)
library(fgsea)

#===============================
# files for ssGSEA
# ==============================
protein_coding = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/protein_coding_ensemble_hgnclist.csv")
gmt_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/BloodGen3_modules_v2.gmt"
BTM <- fgsea::gmtPathways(gmt_path)

#get hgnc gene symbols since BTM has gene symbols
expr_ensg = as.data.frame(rownames(QIV1_all.num))
colnames(expr_ensg) = "target_id"
expr_ensg = inner_join(expr_ensg, protein_coding, by = "target_id")
rownames(QIV1_all.num) = expr_ensg$HGNC_symbol

#===========================================================================
# filtering genes which were only tested for in DEG analysis 15750 (kallisto)
# 14374(salmon) correcting for library batch effect
# ==========================================================================

QIV1_filtered = QIV1_all.num[Washington_MRV1vsLRV1_DEG_results$HGNC_symbol,]
QIV1_adjusted <- sva::ComBat_seq(
  counts = as.matrix(QIV1_filtered),
  batch = QIV1_meta_filt_cond$Library_Batch,
  group =QIV1_meta_filt_cond$responder_group)
QIV1_adjusted <- sva::ComBat_seq(
  counts = as.matrix(QIV1_filtered), 
  batch = QIV1_meta_filt_cond$Library_Batch, 
  group =NULL)

DGEList <- DGEList(QIV1_adjusted) 
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
cpm.norm <- edgeR::cpm(DGEList.norm, log = TRUE)

#==================
# GSVA----
# ==================
param <- GSVA::ssgseaParam(
  exprData = as.matrix(cpm.norm),
  geneSets = BTM,
  minSize = 10,
  maxSize = Inf,
  alpha = 0.25,       
  normalize = TRUE    
)
ssgsea_scores = GSVA::gsva(param)
ssgsea_V1 <- ssgsea_scores[, endsWith(colnames(ssgsea_scores), "V1")]
ssgsea_V2 <- ssgsea_scores[, endsWith(colnames(ssgsea_scores), "V2")]
colnames(ssgsea_V1) <- gsub("_V1$", "", colnames(ssgsea_V1))
colnames(ssgsea_V2) <- gsub("_V2$", "", colnames(ssgsea_V2))

# remove BTMs which do not have a defined name (TBD)
ssgsea_scores_V1 = ssgsea_V1[!grepl("^TBD", rownames(ssgsea_V1)), ]
ssgsea_scores_V2 = ssgsea_V2[!grepl("^TBD", rownames(ssgsea_V2)), ]

write.csv(ssgsea_scores_V1, "ssgsea_scores_V1_null.csv", quote = FALSE)
write.csv(ssgsea_scores_V2, "ssgsea_scores_V2_null.csv", quote = FALSE)


#====================
# Plotting----
# ===================

annotation_col <- QIV1_meta_filt_cond[match(colnames(ssgsea_scores), QIV1_meta_filt_cond$SubjectIDNew), ]
annotation_col <- as.data.frame(annotation_col)
annotation_col_vis <- annotation_col[, "HongKong", drop = FALSE]
rownames(annotation_col_vis) <- annotation_col$SubjectIDNew

annotation_colors <- list(
  HongKong = c("LR" = "#FFA500", "MR" = "#CC9767", "HR" = "#0077B6"))
library(pheatmap)
pdf("ssgsea_V1_null_HK.pdf", width = 14, height = 25)  
pheatmap(
  ssgsea_scores_V1,
  scale = "row",
  annotation_col = annotation_col_vis,
  annotation_colors = annotation_colors,
  fontsize_row = 4,              
  fontsize_col = 8,
  show_rownames = TRUE,
  show_colnames = TRUE,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("blue", "white", "red"))(100),border_color = NA
)
dev.off()






