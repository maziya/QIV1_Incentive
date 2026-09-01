# k can also be set to NA, in which case nemo chooses its value.
# nemo.affinity.graph is the integrated affinity graph.
# 
if (!requireNamespace("SNFtool", quietly = TRUE)) {
  install.packages("SNFtool")
}
library(SNFtool)
devtools::install_github('Shamir-Lab/NEMO/NEMO')
library(dplyr)
library(tibble)
library(NEMO)
library(mclust)
library(ComplexHeatmap)
library(circlize)
library(limma)

QIV1_meta_filt_cond_V2 = QIV1_metadata_processed %>%  filter(Visit == "V2")
proteomics_V2 = column_to_rownames(proteomics_V2,var='feature')
Olink_V2 = column_to_rownames(Olink_V2,var='feature')

valid_subjects <- Reduce(intersect, list(
  colnames(ssgsea_scores_V2), 
  colnames(proteomics_V2), 
  colnames(Olink_V2),
  QIV1_meta_filt_cond_V2$SubjectIDNew
))



# Subset all matrices to the exact same subjects in the exact same order
BTM_sel   <- ssgsea_scores_V2[, valid_subjects]
prot_sel  <- proteomics_V2[, valid_subjects]
olink_sel <- Olink_V2[, valid_subjects]

# Single variance filtering pass for BTM
btm_variance <- apply(BTM_sel, 1, var)
BTM_sel <- BTM_sel[btm_variance > quantile(btm_variance, 0.2), ]

# Scale BTM for clustering and heatmaps
BTM_scaled <- t(scale(t(BTM_sel)))

# ==============================================================================
# 2. NEMO Clustering
# ==============================================================================
omics.list <- list(BTM = BTM_scaled, Proteomics = prot_sel, Olink = olink_sel)

set.seed(42)
affinity.graph <- nemo.affinity.graph(omics.list, k = 10)
num.clusters   <- nemo.num.clusters(affinity.graph)
clustering     <- spectralClustering(affinity.graph, 5)

# Format clustering results
clust.df <- data.frame(
  SubjectID = gsub("V2$", "", colnames(affinity.graph)),
  NEMOCluster = as.factor(clustering)
)

# Merge metadata (ensure responder_groups is uniquely merged)
responder_groups = QIV1_meta_filt_cond_V2[,c(1,7:11)]
clust.df <- inner_join(clust.df, responder_groups, by = "SubjectID")

# ==============================================================================
# 3. Cluster Evaluation (Adjusted Rand Index)
# ==============================================================================
# Helper function for permutation testing
test_ari <- function(cluster_labels, response_labels, n_perm = 1000) {
  actual_ari <- adjustedRandIndex(cluster_labels, response_labels)
  set.seed(42)
  perm_aris <- replicate(n_perm, adjustedRandIndex(cluster_labels, sample(response_labels)))
  pval <- mean(perm_aris >= actual_ari)
  return(list(ARI = actual_ari, P_Value = pval))
}

ari_PH <- test_ari(clust.df$NEMOCluster, clust.df$Phuket)
print(paste("Phuket ARI:", round(ari_PH$ARI, 4), "| p-val:", ari_PH$P_Value))

# ==============================================================================
# 4. Visualization (ComplexHeatmap)
# ==============================================================================
ordered_metadata <- clust.df[order(clust.df$NEMOCluster), ]
rownames(ordered_metadata) <- ordered_metadata$SubjectID

# Reorder matrices to match metadata

colnames(BTM_scaled) = gsub("V2$","", colnames(BTM_scaled))
BTM_plot_data <- BTM_scaled[, rownames(ordered_metadata)]

colnames(olink_sel) = gsub("V2$","", colnames(olink_sel))
olink_plot_data <- olink_sel[, rownames(ordered_metadata)]
olink_plot_data_scaled <- t(scale(t(olink_plot_data)))

colnames(prot_sel) = gsub("V2$","", colnames(prot_sel))
prot_plot_data <- prot_sel[, rownames(ordered_metadata)]
prot_plot_data_scaled <- t(scale(t(prot_plot_data)))


responder_ha <- HeatmapAnnotation(
  Victoria = ordered_metadata$Victoria,
  NEMOCluster = ordered_metadata$NEMOCluster,
  col = list(
    NEMOCluster = c("1"="#009E73", "2"="#CC79A7", "3"="#0072B2", "4"="#F0E442", "5"="#332288"),
    Victoria = c("LR" = "#E69F00", "MR" = "#56B4E9", "HR" = "#D55E00")
  )
)

Olinkheatmap <- Heatmap(
  olink_plot_data_scaled,
  col = colorRamp2(c(-4, 0, 4), c("blue", "grey", "red")),              
  top_annotation = responder_ha,
  column_split = ordered_metadata$NEMOCluster,
  column_gap = unit(2, "mm"),
  cluster_column_slices = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  border = FALSE,
  heatmap_legend_param = list(title = "Olink Protein score scaled Day3",
                              title_gp = gpar(fontsize = 18),   
                              labels_gp = gpar(fontsize = 16))
)

pdf("Olinkheatmap_NEMO_5clustersVC_v2.pdf", width = 16, height = 16)
draw(Olinkheatmap, column_title = " Olink Proteins driving NEMO clusters ",
     column_title_gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()

# ==============================================================================
# 5. Differential Analysis (limma)
# ==============================================================================
design <- model.matrix(~0 + NEMOCluster, data = clust.df)
colnames(design) <- paste0("C", levels(clust.df$NEMOCluster))

fit <- lmFit(BTM_sel, design)

num.clusters = 3 # set as 3 clusters

contrast.matrix <- makeContrasts(C1vC2 = C1 - C2, C1vC3 = C1 - C3, C2vC3 = C2 - C3, levels = design)


fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
C1vC2 <- topTable(fit2, coef = "C1vC2", number = Inf)
C2vC3 <- topTable(fit2, coef = "C2vC3", number = Inf)
C1vC3 <- topTable(fit2, coef = "C1vC3", number = Inf)
