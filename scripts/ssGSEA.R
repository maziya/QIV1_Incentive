gmt_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BTM_for_GSEA_20131008.gmt"
gmt_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/ReactomePathways.gmt"
BTM <- fgsea::gmtPathways(gmt_path)

#get hgnc gene symbols since BTM has gene symbols
expr_ensg = as.data.frame(rownames(QIV1_adjusted))
colnames(expr_ensg) = "target_id"
expr_ensg = inner_join(expr_ensg, protein_coding, by = "target_id")
rownames(QIV1_adjusted) = expr_ensg$HGNC_symbol


QIV1_V1 <- QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V1")]
QIV1_V2 <- QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V2")]

DGEList <- DGEList(QIV1_V1)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
cpm.norm <- edgeR::cpm(DGEList.norm, log = FALSE)
cpm.log <- log2(cpm.norm + 1)

library(GSVA)

param <- GSVA::ssgseaParam(
  exprData = as.matrix(cpm.log),
  geneSets = BTM,
  minSize = 10,
  maxSize = Inf,
  alpha = 0.25,       
  normalize = TRUE    
)


ssgsea_scores = GSVA::gsva(param)

ssgsea_scores_fil = ssgsea_scores[!grepl("^TBA", rownames(ssgsea_scores)), ]

annotation_col <- QIV1_meta_filt_cond[match(colnames(ssgsea_scores), QIV1_meta_filt_cond$SubjectIDNew), ]
annotation_col <- annotation_col[, c("Visit", "SubjectIDNew", "Responder_Victoria")]
rownames(annotation_col) <- annotation_col$SubjectIDNew
annotation_col <- as.data.frame(annotation_col)
annotation_col_vis <- annotation_col[, "ResponderGroup", drop = FALSE]
rownames(annotation_col_vis) <- annotation_col$SubjectIDNew

annotation_colors <- list(
   # Visit = c("V1" = "#FF00FF", "V2" = "#00FFFF"),
   ResponderGroup = c("LR" = "#FFA500", "MR" = "#CC9767", "HR" = "#0077B6"))

pdf("ssgsea_heatmap.pdf", width = 14, height = 25)  
pheatmap(
  ssgsea_scores_fil,
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

library(corrplot)




spearman_corr <- cor(ssgsea_scores, method = "spearman")
spearman_dist <- as.dist(1 - spearman_corr)
pdf("ssgsea_corrplot2.pdf", width = 14, height = 20) 

pheatmap(
  spearman_corr,
  clustering_distance_rows = as.dist(1 - spearman_corr),
  clustering_distance_cols = as.dist(1 - spearman_corr),
  annotation_col = annotation_col_vis,
  annotation_colors = annotation_colors,
  clustering_method = "average",
  show_rownames = FALSE, show_colnames = FALSE,
  color = colorRampPalette(c("blue","white","red"))(100),
  legend = TRUE,
  border_color = NA
)
dev.off()


