#only protein coding genes from count matrix using gencode data

QIV1_V2.num = QIV1_V2.num[rownames(QIV1_V2.num) %in% gf_id, , drop = FALSE]


#plots on the logcpm values with and without filtering
#and normalization with TMM----
myDGEList <- DGEList(QIV1.num.V1)
log2.cpm <- cpm(myDGEList, log=TRUE)
sampleLabels = colnames(log2.cpm)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = 2:ncol(log2.cpm.df),
                                  names_to = "samples",
                                  values_to = "expression")

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 124,
               size = 3,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM) HR & NR",
       subtitle="only protein coding unfiltered, non-normalized") +
  theme_bw() + coord_flip()
ggsave('log2cpm_HR_NR.png',plot = p1,width = 20, height = 15,dpi = 300)

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=32 #genes which have greater counts per million in at least 32 samples
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = 2:ncol(log2.cpm.filtered.df),
                                           names_to = "samples",
                                           values_to = "expression")

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 124,
               size = 3,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM) V2",
       subtitle="filtered, non-normalized") +
  theme_bw()+coord_flip()
ggsave('log2cpm_filtered_cpm1_V2.png',plot = p2,width = 20, height = 15,dpi = 300)

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = 2:ncol(log2.cpm.filtered.norm.df),
                                                names_to = "samples",
                                                values_to = "expression")


p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 124,
               size = 3,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM) V2",
       subtitle="filtered,TMM normalized") +
  theme_bw()+coord_flip()

ggsave('log2cpm_filtered_norm_V2.png',plot = p3,width = 20, height = 15,dpi = 300)

#without any filtering with cpm
DGEList <- DGEList(QIV1_V1.num)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
log2.cpm.norm <- cpm(DGEList.norm, log=TRUE)
log2.cpm.norm.df <- as_tibble(log2.cpm.norm, rownames = "geneID")
colnames(log2.cpm.norm.df) <- c("geneID", sampleLabels)
log2.cpm.norm.df.pivot <- pivot_longer(log2.cpm.norm.df,
                                                cols = 2:ncol(log2.cpm.norm.df),
                                                names_to = "samples",
                                                values_to = "expression")


p4 <- ggplot(log2.cpm.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 124,
               size = 3,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM) V1",
       subtitle="protein coding,TMM normalized") +
  theme_bw()+coord_flip()

ggsave('log2cpm_proteincoding_norm_V1.png',plot = p4,width = 20, height = 15,dpi = 300)

#hierarchical clustering ----
#function to color labels in dendogram
set_labels_colors <- function(dend, labels_to_color1, color1, labels_to_color2, color2) {
  all_labels <- labels(dend)
  new_labels <- ifelse(all_labels %in% labels_to_color1 | all_labels %in% labels_to_color2, all_labels, "")
  label_colors <- ifelse(all_labels %in% labels_to_color1, color1,
                         ifelse(all_labels %in% labels_to_color2, color2, "black"))
  dend <- set(dend, "labels", new_labels)
  dend <- set(dend, "labels_col", label_colors)
  return(dend)
}

clusters <- as.dendrogram(hclust(dist(t(log2.cpm.norm), method = "maximum"),method = 'complete'))
library(dendextend)
responder_ids <- QIV1_meta_filt_V1 %>%
  dplyr::filter(Status == "Responder") %>%
  dplyr::select(SubjectID)%>%
  pull()
non_responder_ids <- QIV1_meta_filt_V1 %>%
  dplyr::filter(Status == "NonResponder") %>%
  dplyr::select(SubjectID)%>%
  pull()

dend <- set_labels_colors(clusters, responder_ids, "purple", non_responder_ids, "orange")
png("samples_dendogram_proteincoding_norm_V1.png",units="in", width=15,height=10,res=600)
plot(dend)
dev.off()


#PCA plots ----

sampleLabels <- colnames(log2.cpm.norm)
group <- QIV1_meta_filt_V2$Status
names(group) <- QIV1_meta_filt_V2$SubjectID
group <- group[sampleLabels]
group <- factor(group)


pca.res <- prcomp(t(log2.cpm.norm), scale. = FALSE, retx = TRUE)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)
pca.res.df <- as_tibble(pca.res$x) %>%
  mutate(group = group, sampleLabels = sampleLabels)

pca.plot <- ggplot(pca.res.df, aes(x = PC1, y = PC2, color = group, label = sampleLabels)) +
  geom_point(size = 1) +
  stat_ellipse(aes(group = group), type = "norm") +
  # geom_label(size = 2.5, label.size = 0.15, alpha = 0.8) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA QIV1-HR NR both timepoints",
       subtitle = "allgenes, TMM normalized") +
  coord_fixed() +
  theme_bw()
ggsave('pca_HR_NR_bothtimepoints.png',plot = pca.plot, width = 15, height=10)

#===============================================================================================
#
QIV1_meta_filt = QIV1_meta_filt %>% dplyr::filter(Status %in% c("HR","NR"))

QIV1_all.num = QIV1_all.num[, (colnames(QIV1_all.num) %in% QIV1_meta_filt$SubjectIDNew)]

v1_cols <- grep("V1$", colnames(QIV1_all.num), value = TRUE)
v2_cols <- grep("V2$", colnames(QIV1_all.num), value = TRUE)


QIV1.num.V1 <- QIV1_all.num[, v1_cols]
QIV1.num.V2 <- QIV1_all.num[, v2_cols]
myDGEList <- DGEList(QIV1.num.V1)
myDGEList.norm <- calcNormFactors(myDGEList, method = "TMM")
log2.cpm.norm <- cpm(myDGEList.norm, log=TRUE)


library(tidyverse)


meta_long <- QIV1_meta_filt %>%
  dplyr::select(SubjectIDNew, SEX, Library_Batch, starts_with("Responder_")) %>%
  tidyr::pivot_longer(cols = starts_with("Responder_"),
                      names_to = "Strain",
                      values_to = "ResponderStatus") %>%
  mutate(Strain = str_remove(Strain, "^Responder_"))

samples_to_label <- c("QIV1KEM033V1", "QIV1KEM005V1", "QIV1KEM029V1","QIV1KEM084V1","QIV1KEM094V1")
samples_to_label <- c("QIV1KEM033V2", "QIV1KEM005V2", "QIV1KEM029V2","QIV1KEM084V2","QIV1KEM094V2")


sampleLabels <- colnames(log2.cpm.top50)
pca.res <- prcomp(t(log2.cpm.top50), scale. = FALSE, retx = TRUE)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)

pca.res.df <- as_tibble(pca.res$x) %>%
  mutate(SubjectIDNew = sampleLabels)


pca_long <- pca.res.df %>%
  left_join(meta_long, by = "SubjectIDNew") %>%
  mutate(ResponderStatus = factor(ResponderStatus,
                                  levels = c("Responder", "NonResponder")))
library(ggrepel)

pca.plot <- ggplot(pca_long, aes(x = PC1, y = PC2, color = ResponderStatus)) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = ResponderStatus), type = "norm") +
  geom_text_repel(
    data = pca_long %>% dplyr::filter(SubjectIDNew %in% samples_to_label),
    aes(label = SubjectIDNew),
    size = 3,
    max.overlaps = Inf
  ) +
  facet_wrap(~ Strain) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA: Responder Status per Strain V1",
       subtitle = "Top 50 genes(std.dev)log2 CPM, TMM normalized") +
  coord_fixed() +
  theme_bw()

ggsave("top50genesPCA_by_strain_responder_statusV1.png", plot = pca.plot, width = 16, height = 10)


#standard deviation for each gene
gene_sd = apply(log2.cpm.norm, 1, sd)
top50_genes = names(sort(gene_sd, decreasing = TRUE))[1:50]
top50 = sort(gene_sd, decreasing = TRUE)

#expr matrix of the top50
log2.cpm.top50 = log2.cpm.norm[top50_genes, ]


##################################################################
##################################################################
DGEList <- DGEList(QIV1_V2.num)
DGEList.norm <- calcNormFactors(DGEList, method = "TMM")
log2.cpm.norm <- cpm(DGEList.norm, log=TRUE)

meta_long <- QIV1_meta_filt %>%
  dplyr::select(SubjectIDNew, SEX, Library_Batch, starts_with("Responder_")) %>%
  tidyr::pivot_longer(cols = starts_with("Responder_"),
                      names_to = "Strain",
                      values_to = "ResponderStatus") %>%
  mutate(Strain = str_remove(Strain, "^Responder_"))

sampleLabels <- colnames(log2.cpm.norm)
pca.res <- prcomp(t(log2.cpm.norm), scale. = FALSE, retx = TRUE)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)

pca.res.df <- as_tibble(pca.res$x) %>%
  mutate(SubjectIDNew = sampleLabels)


pca_long <- pca.res.df %>%
  left_join(meta_long, by = "SubjectIDNew") %>%
  mutate(ResponderStatus = factor(ResponderStatus,
                                  levels = c("Responder", "NonResponder")),
         SEX = factor(SEX, levels = c("M","F")),
         Library_Batch = factor(Library_Batch))
library(ggrepel)

pca.plot = ggplot(pca_long, aes(x = PC1, y = PC2, color = ResponderStatus)) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = ResponderStatus), type = "norm") +
  facet_wrap(~ Strain) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA: V2 Colored by ResponderStatus",
       subtitle = "All genes log2 CPM, TMM normalized") +
  coord_fixed() +
  theme_bw()

ggsave("genesPCA_V2_ResponderStatus.png", plot = pca.plot, width = 16, height = 10)

pca_df = pca_long %>%
  mutate(IsOutlier = ifelse(abs(PC1) > 200 | abs(PC2) > 100, TRUE, FALSE))

#sample QIV1KEM043V1 also has low gene counts