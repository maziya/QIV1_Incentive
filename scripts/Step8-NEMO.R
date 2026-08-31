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

valid_subjects <- Reduce(intersect, list(
  colnames(ssgsea_scores_V1), 
  colnames(proteomics_V1), 
  colnames(Olink_V1),
  QIV1_meta_filt_cond_V1$SubjectIDNew
))

# Subset all matrices to the exact same subjects in the exact same order
BTM_sel   <- ssgsea_scores_V1[, valid_subjects]
prot_sel  <- proteomics_V1[, valid_subjects]
olink_sel <- Olink_V1[, valid_subjects]

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
clustering     <- spectralClustering(affinity.graph, num.clusters)

# Format clustering results
clust.df <- data.frame(
  SubjectID = gsub("V2$", "", colnames(affinity.graph)),
  NEMOCluster = as.factor(clustering)
)

# Merge metadata (ensure responder_groups is uniquely merged)
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

ari_HK <- test_ari(clust.df$NEMOCluster, clust.df$HongKong)
print(paste("HongKong ARI:", round(ari_HK$ARI, 4), "| p-val:", ari_HK$P_Value))

# ==============================================================================
# 4. Visualization (ComplexHeatmap)
# ==============================================================================
ordered_metadata <- clust.df[order(clust.df$NEMOCluster), ]
rownames(ordered_metadata) <- ordered_metadata$SubjectID

# Reorder matrices to match metadata
BTM_plot_data <- BTM_scaled[, rownames(ordered_metadata)]

responder_ha <- HeatmapAnnotation(
  HongKong = ordered_metadata$HongKong,
  NEMOCluster = ordered_metadata$NEMOCluster,
  col = list(
    NEMOCluster = c("1"="#009E73", "2"="#CC79A7", "3"="#0072B2", "4"="#F0E442"),
    HongKong = c("LR" = "#E69F00", "MR" = "#56B4E9", "HR" = "#D55E00")
  )
)

BTMheatmap <- Heatmap(
  BTM_plot_data,
  col = colorRamp2(c(-4, 0, 4), c("blue", "grey", "red")),              
  top_annotation = responder_ha,
  column_split = ordered_metadata$NEMOCluster,
  cluster_column_slices = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2"
)

pdf("BTM_NEMO_clusters_v1.pdf", width = 16, height = 14)
draw(BTMheatmap, column_title = "BTMs driving NEMO clusters")
dev.off()

# ==============================================================================
# 5. Differential Analysis (limma)
# ==============================================================================
design <- model.matrix(~0 + NEMOCluster, data = clust.df)
colnames(design) <- paste0("C", levels(clust.df$NEMOCluster))

fit <- lmFit(BTM_sel, design)

# Dynamically handle 2-cluster vs 3-cluster outputs
if (num.clusters == 2) {
  contrast.matrix <- makeContrasts(C1vC2 = C1 - C2, levels = design)
} else if (num.clusters == 3) {
  contrast.matrix <- makeContrasts(C1vC2 = C1 - C2, C1vC3 = C1 - C3, C2vC3 = C2 - C3, levels = design)
}

fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
C1vC2 <- topTable(fit2, coef = "C1vC2", number = Inf)
