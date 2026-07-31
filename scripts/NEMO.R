# k can also be set to NA, in which case nemo chooses its value.
# nemo.affinity.graph is the integrated affinity graph.
# 
if (!requireNamespace("SNFtool", quietly = TRUE)) {
install.packages("SNFtool")
}
library(SNFtool)
devtools::install_github('Shamir-Lab/NEMO/NEMO')
library(NEMO)

BTM_wide = pivot_wider(BTM_long, names_from = sample, values_from = value)
BTM_wide = BTM_wide[,-c(2)]
BTM_wide = column_to_rownames(BTM_wide, var="Feature")
btm_variance <- apply(BTM_wide, 1, var)
variance_threshold <- quantile(btm_variance, 0.2)
BTM_sel = BTM_wide[btm_variance > variance_threshold, ]

proteomics_sel = proteomics_V1[,which(colnames(proteomics_V1)%in%colnames(BTM_sel))]
rownames(proteomics_sel)= rownames(proteomics_V1)
Olink_sel = Olink_V1[,which(colnames(Olink_V1)%in%colnames(BTM_sel))]
rownames(Olink_sel) = Olink_V1$Protein
metabolomics_sel = metabolomics_V1[,which(colnames(metabolomics_V1)%in%colnames(BTM_sel))]
rownames(metabolomics_sel) = metabolomics_V1$Metabolite

omics.list = list(BTM_sel,proteomics_sel,Olink_sel)


set.seed(42)
affinity.graph = nemo.affinity.graph(omics.list, k = 10)

  # nemo to estimate the number of clusters
num.clusters = nemo.num.clusters(affinity.graph)

# clustering is the cluster assignment vector, ordered by the columns of affinity.graph.
clustering = spectralClustering(affinity.graph, num.clusters)
names(clustering) = colnames(affinity.graph)

clust.df = as.data.frame(clustering)
clust.df = rownames_to_column(clust.df,var="SubjectID")
clust.df$SubjectID = gsub("V1$","", clust.df$SubjectID)
clust.df = inner_join(clust.df,ResponderGroups, by = "SubjectID")
clust.df = inner_join(clust.df,responder_groups, by = "SubjectID")

clustering2 = spectralClustering(affinity.graph, 3)
names(clustering2) = colnames(affinity.graph)
clust.df2 = as.data.frame(clustering2)
clust.df2 = rownames_to_column(clust.df2,var="SubjectID")
clust.df2$SubjectID = gsub("V1$","", clust.df2$SubjectID)
clust.df2 = inner_join(clust.df2,ResponderGroups, by = "SubjectID")
clust.df2 = inner_join(clust.df2,responder_groups, by = "SubjectID")


library(mclust) 

# clustering: named vector from spectralClustering, e.g. clustering["P1"] = 1
# labels:  HR/MR/LR labels, same patient order/names

ari_HK <- adjustedRandIndex(clust.df$clustering, clust.df$HongKong)
print(ari_HK)
set.seed(42)
perm_aris <- replicate(1000, adjustedRandIndex(clust.df$clustering, sample(clust.df$HongKong)))
p_value <- mean(perm_aris >= ari_HK)

adjustedRandIndex(clust.df$clustering,clust.df$HongKong)
adjustedRandIndex(clust.df$clustering, clust.df$Victoria)
adjustedRandIndex(clust.df$clustering, clust.df$Phuket)
adjustedRandIndex(clust.df$clustering, clust.df$Washington)



library(pheatmap)
cont_table <- table(fused.clusters$clustering_k2, fused.clusters$HongKong)
pheatmap(cont_table, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, main = "NEMO clusters vs Responder status")


clustering_k2 <- spectralClustering(affinity.graph, 2)
names(clustering_k2) <- colnames(affinity.graph)
table(clustering_k2)
names(clustering_k2) = colnames(affinity.graph)
fused.clusters = as.data.frame(clustering_k2)
fused.clusters = rownames_to_column(fused.clusters,var="SubjectID")
fused.clusters$SubjectID = gsub("V1$","", fused.clusters$SubjectID)
fused.clusters = inner_join(fused.clusters,ResponderGroups, by = "SubjectID")
fused.clusters = inner_join(fused.clusters,responder_groups, by = "SubjectID")


#clustering per omics layer
per.omic.graphs <- lapply(names(omics.list), function(name) {
  nemo.affinity.graph(list(omics.list[[name]]), k = 5)
})
names(per.omic.graphs) <- names(omics.list)

num.clusters <- 3 
set.seed(42) 
per.omic.clusters <- lapply(per.omic.graphs, function(aff) {
  spectralClustering(aff, K = num.clusters)
})
names(per.omic.clusters) <- names(omics.list)
set.seed(42) 
fused.aff <- nemo.affinity.graph(omics.list, k = 10)
fused.clusters <- spectralClustering(fused.aff, K = num.clusters)



names(fused.clusters) <- colnames(fused.aff)
table(fused.clusters)
names(fused.clusters) = colnames(fused.aff)
fused.clusters.df = as.data.frame(fused.clusters)
fused.clusters.df = rownames_to_column(fused.clusters.df,var="SubjectID")
fused.clusters.df$SubjectID = gsub("V1$","", fused.clusters.df$SubjectID)
fused.clusters.df = inner_join(fused.clusters.df,ResponderGroups, by = "SubjectID")
fused.clusters.df = inner_join(fused.clusters.df,responder_groups, by = "SubjectID")

adjustedRandIndex(fused.clusters.df$fused.clusters,fused.clusters.df$HongKong)
adjustedRandIndex(fused.clusters.df$fused.clusters, fused.clusters.df$Victoria)
adjustedRandIndex(fused.clusters.df$fused.clusters, fused.clusters.df$Phuket)
adjustedRandIndex(fused.clusters.df$fused.clusters, fused.clusters.df$Washington)



library(mclust)
comparisons <- expand.grid(names(per.omic.clusters), names(per.omic.clusters))
sapply(seq_len(nrow(comparisons)), function(i) {
  a <- comparisons[i, 1]; b <- comparisons[i, 2]
  adjustedRandIndex(per.omic.clusters[[a]], per.omic.clusters[[b]])
})

sapply(per.omic.clusters, function(cl) adjustedRandIndex(cl, fused.clusters))


library(pheatmap)

max.aff <- max(sapply(per.omic.graphs, max))
breaks.common <- seq(0, max.aff, length.out = 101)
for (name in names(omics.list)) {
  aff <- per.omic.graphs[[name]]
  cl <- per.omic.clusters[[name]]
  ord <- order(cl)
  
  ann <- data.frame(cluster = factor(cl[ord]))
  rownames(ann) <- rownames(aff)[ord]
  
  pheatmap(aff[ord, ord],
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = ann, annotation_col = ann,
           show_rownames = FALSE, show_colnames = FALSE,
           main = paste("Affinity graph -", name),
           breaks = breaks.common,
          filename = paste0("affinity_graph_", name, ".png"),
          width = 6, height = 5.5, dpi = 300)
  
}
# order patients by fused cluster assignment
ord.fused <- order(fused.clusters)

ann.fused <- data.frame(cluster = factor(fused.clusters[ord.fused]))
rownames(ann.fused) <- rownames(aff)[ord.fused]

pheatmap(fused.aff[ord.fused, ord.fused],
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = ann.fused, annotation_col = ann.fused,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Affinity graph - Fused",
         breaks = breaks.common,
         filename = "affinity_graph_Fused.png",
         width = 6, height = 5.5, dpi = 300)


heatmap_colors <- colorRamp2(c(-4, 0, 4), c("blue", "grey", "red"))
library(circlize)
library(ComplexHeatmap)
group_colors = c("1"="blue", "2"="orange")
group_colors = c("1"="#009E73", "2"="#CC79A7","3" = "#0072B2", "4"= "#F0E442")
responder_colors <- c("LR" = "#E69F00", "MR" = "#56B4E9", "HR" = "#D55E00")
annotation_col$ResponderGroup = factor(annotation_col$ResponderGroup)
annotation_col <- clust.df
rownames(annotation_col) <-  clust.df$SubjectID

BTM_sel_scaled <- t(scale(t(BTM_sel)))

ordered_metadata <- annotation_col[order(annotation_col$clustering), ]
BTM_sel_scaled <- BTM_sel_scaled[, rownames(ordered_metadata)]

responder_ha <- HeatmapAnnotation(
  ResponderGroup = ordered_metadata[, c(
    #"Victoria",
    #"HongKong",
    #"Phuket",
    #"Washington",
    "ResponderGroup"
  )],
  NEMOCluster = ordered_metadata$clustering,
  col = list(NEMOCluster = group_colors,
             #Victoria   = responder_colors,
             #HongKong   = responder_colors,
             #Phuket     = responder_colors,
             #Washington = responder_colors,
             ResponderGroup = responder_colors
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    NEMOCluster = list(
      title_gp = gpar(fontsize = 18),
      labels_gp = gpar(fontsize = 16)
    )
  )
)

BTMheatmap = Heatmap(
  BTM_sel_scaled,
  col = heatmap_colors,              
  top_annotation = responder_ha,
  column_split = ordered_metadata$clustering,
  column_gap = unit(2, "mm"),
  cluster_column_slices = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  #clustering_distance_rows = "pearson",
  #clustering_distance_columns = "pearson",
  
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  border = FALSE,
  heatmap_legend_param = list(title = "BTM scores Day0",
                              title_gp = gpar(fontsize = 18),   
                              labels_gp = gpar(fontsize = 16))
)
pdf("BTM_NEMO_clusters_v4.pdf", width = 16, height = 14)
draw(BTMheatmap, column_title = "BTMs driving NEMO clusters ",
     column_title_gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()
