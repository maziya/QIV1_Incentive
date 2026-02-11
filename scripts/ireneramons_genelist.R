#using the gene list from IreneRamons review paper to check 
#the expression profile in V1 and V2 with heatmaps
QIV1_V1 <- QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V1")]
QIV1_V2 <- QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V2")]

DGEList <- DGEList(QIV1_V1)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
cpm.norm.log <- edgeR::cpm(DGEList.norm, log = TRUE)
QIV1_V1.df = as.data.frame(cpm.norm.log)
QIV1_V1.df = rownames_to_column(QIV1_V1.df, var="target_id")
QIV1_V1.df = left_join(QIV1_V1.df,protein_coding, by ="target_id")
QIV1_V1.df = column_to_rownames(QIV1_V1.df, var="HGNC_symbol")
QIV1_V1.df = QIV1_V1.df[,-c(1)]

#extract gene exp for all signatures
baseline_signatures_various <- read_csv("baseline_signatures_various.csv")

QIV1_V1_immune = as.list()
QIV1_V1_immune = lapply(colnames(baseline_signatures_various), function(sig) {
  
  genes = baseline_signatures_various[[sig]]
  genes = genes[genes %in% rownames(QIV1_V1.df)]
  QIV1_V1.df[genes, , drop = FALSE]
})

names(QIV1_V1_immune) <- colnames(baseline_signatures_various)
QIV1_sig_mean_list <- lapply(names(QIV1_V1_immune), function(sig) {
  
  mat = QIV1_V1_immune[[sig]]
  mat_scaled <- t(scale(t(mat)))
  
  data.frame(
    SubjectIDNew = colnames(mat_scaled),
    mean_expression = colMeans(mat_scaled),
    Signature = sig
  )
})

QIV1_sig_mean = bind_rows(QIV1_sig_mean_list)

QIV1_sig_mean <- left_join(
  QIV1_sig_mean,
  QIV1_meta_filt_cond,
  by = "SubjectIDNew"
)
comparisons <- list(
  c("LR","MR"),
  c("MR","HR"),
  c("LR","HR")
)

library(ggpubr)
sig_mean_box <- ggplot(QIV1_sig_mean,
                         aes(x = ResponderGroup, y = mean_expression, fill = ResponderGroup)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.8, color = "black") +
  theme_classic(base_size = 16) + facet_wrap(~ Signature) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14)
  ) +labs(title = "Baseline Signature Expression in QIV1 Day 0",
          x = "Responder Status",
          y = "Mean Z-scored gene expression")+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    step.increase = 0.12
  )
ggsave("Baseline_signatures_QIV1_Day0.png",plot = sig_mean_box, dpi = 300, width =20, height = 15)


##################################################################################
immunegenes <- immunegenes_curatedlist %>%
  pivot_longer(cols = everything(), names_to = "column", values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene)
protein_coding = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/protein_coding_ensemble_hgnclist.csv")

immunegenes = as.data.frame(immunegenes)


genes_to_keep = immunegenes$gene[immunegenes$gene %in% rownames(QIV1_V1.df)]
others = immunegenes$gene[!immunegenes$gene %in% rownames(QIV1_V1.df)]
QIV1_V1_immune = QIV1_V1.df[genes_to_keep,]

QIV1_V1_immune = as.matrix(QIV1_V1_immune)

QIV1_meta_filt_cond = as.data.frame(QIV1_meta_filt_cond)
annotation_col <- QIV1_meta_filt_cond[match(colnames(QIV1_V1_immune), QIV1_meta_filt_cond$SubjectIDNew), ]
rownames(annotation_col) <- annotation_col$SubjectIDNew
annotation_col <- annotation_col[, c("Responder_Victoria","Responder_HongKong","Responder_Phuket","Responder_Washington","diabetic","SEX","covid_vaccinated",
                                     )]  # add "ResponderGroup" if using adjMFC category


responder_colors <- c("LR" = "#E69F00", "MR" = "#56B4E9", "HR" = "#D55E00")
diabetic_colors = c("Diabetic" = "pink", "Non_diabetic" = "green")
sex_colors = c("M"="blue", "F"="orange")
covid_colors = c("Covishield"="cyan", "Covaxin"="red","No_covid_vaccine"="gray")
# group_colors = c("LR" = "#E69F00", "MR" = "#56B4E9", "HR" = "#D55E00")
annotation_colors <- list(
  Responder_Victoria = responder_colors,
  Responder_HongKong = responder_colors,
  Responder_Phuket = responder_colors,
  Responder_Washington = responder_colors,
  
)

pdf("immunegene_list_V1.pdf", width = 14, height = 25)  
pheatmap(
  QIV1_V1_immune,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 4,              
  fontsize_col = 8,
  show_rownames = TRUE,
  show_colnames = FALSE,
  legend = TRUE,main = "Gene expression(log-transformed TMM)",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100),border_color = NA
)
dev.off()
responder_ha <- HeatmapAnnotation(
  df = annotation_col[, c(
    "Responder_Victoria",
    "Responder_HongKong",
    "Responder_Phuket",
    "Responder_Washington",
    "diabetic",
    "SEX","covid_vaccinated"
  )], #add "ResponderGroup" if comparing adjMFC groups
  
  col = list(
    Responder_Victoria   = responder_colors,
    Responder_HongKong   = responder_colors,
    Responder_Phuket     = responder_colors,
    Responder_Washington = responder_colors,
    diabetic             = diabetic_colors,
    SEX                  = sex_colors,
    covid_vaccinated     = covid_colors
    # ResponderGroup       = group_colors
  ),
  
  annotation_name_side = "left",
  which = "column"
)


library(circlize)
heatmap_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
QIV1_V1_immune_scaled <- t(scale(t(QIV1_V1.df)))

immuneheatmap = Heatmap(
  QIV1_V1_immune_scaled,
  name = "Expression",               
  col = heatmap_colors,              
  top_annotation = responder_ha,           
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  border = FALSE,
  heatmap_legend_param = list(title = "Gene expression (log-TMM)")
)
pdf("immunegene_list_V1_clustered2.pdf", width = 14, height = 25)
draw(immuneheatmap)
dev.off()

#extract dendogram from heatmap
ht <- draw(immuneheatmap)
col_dend <- column_dend(ht)
col_hclust <- as.hclust(col_dend)
col_clusters <- cutree(col_hclust, k = 3)
cluster1_ids = names(col_clusters[col_clusters == 1])
cluster2_ids = names(col_clusters[col_clusters == 2])
cluster3_ids = names(col_clusters[col_clusters == 3])


#make df from the luster info
cluster_df <- data.frame(
  SampleID = names(col_clusters),
  Cluster  = col_clusters
)
annot_with_cluster <- annotation_col %>%
  mutate(SampleID = rownames(annotation_col)) %>%
  inner_join(cluster_df, by = "SampleID")
long_df <- annot_with_cluster %>%
  pivot_longer(
    cols = c(Responder_Victoria, Responder_HongKong,
             Responder_Phuket, Responder_Washington),
    names_to = "Strain",
    values_to = "Responder"
  )
#calculate percentage of responder categories in each cluster per strain
percent_df <- long_df %>%
  group_by(Cluster, Strain, Responder) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster, Strain) %>%
  mutate(Percent = n / sum(n) * 100)

#plot bar chart
clusterbar = ggplot(percent_df, aes(x = Strain, y = Percent, fill = Responder)) +
  geom_col(position = "stack") +
  facet_wrap(~ Cluster, nrow = 1) +
  labs(
    title = "Responder distribution within each cluster",
    x = "Strain",
    y = "Percentage"
  ) +
  scale_fill_manual(
    values = c(LR = "#E69F00", MR = "#56B4E9", HR = "#D55E00")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("clusters_baseline_immunegenes1.png", plot = clusterbar, width = 15, height = 10, dpi = 300)

percent_df <- long_df %>%
  mutate(Cluster = as.factor(Cluster)) %>%   
  group_by(Cluster, Responder) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Percent = n / sum(n) * 100)


clusterbar <- ggplot(percent_df,
                     aes(x = Cluster, y = Percent, fill = Responder)) +
  geom_col(position = "fill") +   
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Responder distribution aggregated across all strains",
    x = "Cluster",
    y = "Percentage"
  ) +
  scale_fill_manual(
    values = c(LR = "#E69F00", MR = "#56B4E9", HR = "#D55E00")
  ) +
  theme_classic(base_size = 14)
ggsave("clusters_across_strains1.png", 
       plot = clusterbar, width = 15, height = 10, dpi = 300)



collapsed_df <- long_df %>%
distinct(Cluster, SampleID, Responder)
counts_df <- collapsed_df %>%
group_by(Cluster, Responder) %>%
summarise(n = n(), .groups = "drop") %>%
tidyr::pivot_wider(
names_from = Responder,
values_from = n,
values_fill = 0
) %>%
mutate(Total = LR + MR + HR)

run_fisher_test <- function(df, c1, c2, responder) {
df1 <- df %>% filter(Cluster == c1)
df2 <- df %>% filter(Cluster == c2)
# Build the 2x2 table
mat <- matrix(
c(df1[[responder]],
df1$Total - df1[[responder]],
df2[[responder]],
df2$Total - df2[[responder]]),
nrow = 2,
byrow = TRUE
)
test <- fisher.test(mat)
tibble(
Cluster1  = c1,
Cluster2  = c2,
Responder = responder,
p_value   = test$p.value
)
}
clusters <- sort(unique(counts_df$Cluster))
responders <- c("LR", "MR", "HR")
results <- purrr::map_dfr(
combn(clusters, 2, simplify = FALSE),
\(pair) purrr::map_dfr(
responders,
\(resp) run_fisher_test(counts_df, pair[1], pair[2], resp)
)
)
results_fdr <- results %>%
mutate(FDR = p.adjust(p_value, method = "fdr"))




deg_results <- results_list[["HongKong_HRV1vsLRV1"]]
 deg_results$target_id <- rownames(deg_results)
 all_genes <- left_join(deg_results, protein_coding, by = "target_id")
 all_genes$significance <- ifelse(
      all_genes$P.Value < 0.01 & abs(all_genes$logFC) > 1,
      ifelse(all_genes$logFC > 1, "Upregulated", "Downregulated"),
       "Not significant")
 all_genes_predictor = all_genes %>% filter(HGNC_symbol %in% immunegenes$HGNC_symbol)
HK_immunegeneDEG_V1 = all_genes_predictor

 

 
df_all <- HK_immunegeneDEG_V1 %>%
   inner_join(VC_immunegeneDEG_V1, by = "HGNC_symbol", suffix = c("_HK", "_VC")) %>%
   inner_join(PH_immunegeneDEG_V1, by = "HGNC_symbol") %>%
   rename(logFC_PH = logFC, P.Value_PH = P.Value) %>%
   inner_join(WH_immunegeneDEG_V1, by = "HGNC_symbol") %>%
   rename(logFC_WH = logFC, P.Value_WH = P.Value)
 
pos_all <- df_all %>% filter(if_all(starts_with("logFC_"), \(x) x > 0))
neg_all <- df_all %>% filter(if_all(starts_with("logFC_"), \(x) x < 0))


#TGsig expression in QIV1
TGsig = immunegenes %>% filter(gene %in% c("C15orf57","LONP2","PAPSS2","EPHB1","ADAM12","SMC1A","RETN","ENPP1","CD101","C2orf63"))

QIV1_tgsig = as.data.frame(QIV1_V1_immune_scaled) %>% filter(rownames(QIV1_V1_immune_scaled) %in% TGsig$gene)

QIV1_tgsig_long <- QIV1_tgsig %>%
  mutate(gene = rownames(QIV1_tgsig)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "SubjectIDNew",
    values_to = "expression"
  ) %>%
  left_join(QIV1_meta_filt_cond, by = "SubjectIDNew")
QIV1_tgsig_plot <- QIV1_tgsig_long %>%
  pivot_longer(
    cols = starts_with("Responder_"),
    names_to = "Virus",
    values_to = "Response"
  ) %>%
  mutate(
    Virus = gsub("Responder_", "", Virus)
  )


QIV1_tgsig_mean <- QIV1_tgsig %>%
  mutate(gene = rownames(QIV1_tgsig)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "SubjectIDNew",
    values_to = "expression"
  ) %>%
  group_by(SubjectIDNew) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
  left_join(QIV1_meta_filt_cond, by = "SubjectIDNew") %>%
  select(SubjectIDNew,
         mean_expression,
         Responder_Victoria,
         Responder_HongKong,
         Responder_Phuket,
         Responder_Washington,
         ResponderGroup)
#add ResponderGroup if using adjMFC labels

QIV1_tgsig_mean_facet <- QIV1_tgsig_mean %>%
  pivot_longer(
    cols = starts_with("Responder_"),
    names_to = "Virus",
    values_to = "Response"
  ) %>%
  mutate(
    Virus = sub("Responder_", "", Virus),
    Response = factor(Response)
  )


comparisons <- list(
  c("LR","MR"),
  c("MR","HR"),
  c("LR","HR")
)

#boxplot of genes individually
tgsig_box = ggplot(QIV1_tgsig_plot,   
  aes(x = gene, y = expression, fill = Response)) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~ Virus, ncol = 2) +
  theme_bw(base_size = 16) +
  theme(strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 14)
  )
ggsave("tgsig_exp_QIV1_day0_v2.png",plot = tgsig_box, dpi = 300, width =20, height = 15)

#boxplot of mean expression of all genes
tgsig_mean_box <- ggplot(QIV1_tgsig_mean,
aes(x = ResponderGroup, y = mean_expression, fill = ResponderGroup)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.8, color = "black") +
  # facet_wrap(~ Virus, ncol = 2) +
  theme_classic(base_size = 16) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14)
  ) +labs(title = "QIV1 signature(from 19kgene list) average expression in QIV1 Day0",
          x = "Responder Status",
          y = "Average gene expression")+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    step.increase = 0.12
  )
ggsave("QIV1sig_baseline_bestpredictorlist_from19klist.png",plot = tgsig_mean_box, dpi = 300, width =20, height = 15)

