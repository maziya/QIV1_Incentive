pca.res <- prcomp(QIV1_responder_BL, scale. = FALSE, retx = TRUE)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)

pca.res.df <- as_tibble(pca.res$x) %>%
  mutate(SubjectID = rownames(QIV1_responder_BL))



library(ggrepel)

pca.plot = ggplot(pca.res.df, aes(x = PC1, y = PC2)) +
  geom_point(size = 2) + geom_text_repel(aes(label = SubjectID), size = 3)+
  # stat_ellipse(aes(group = Library_Batch), type = "norm") +
  # facet_wrap(~ Strain) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA: baseline response") +
  coord_fixed() +
  theme_bw()

ggsave("PCA_baselinetitres.png", plot = pca.plot, width = 16, height = 10)


set.seed(123)
km <- kmeans(QIV1_responder_BL, centers = 4)  

km$cluster   # cluster assignment for each subject
pca.res <- prcomp(QIV1_responder_BL, scale. = TRUE)
pc_df <- as.data.frame(pca.res$x) %>%
  mutate(SubjectID = rownames(QIV1_responder_BL),
         Cluster = factor(km$cluster))

ggplot(pc_df, aes(x = PC1, y = PC2, color = Cluster, label = SubjectID)) +
  geom_point(size=3) +
  geom_text(vjust=-1) +
  theme_bw() +
  labs(title = "PCA with k-means clustering")
ggsave("PCA_baselinetitres.png", plot = pca.plot, width = 16, height = 10)


#========================================================================
#using only binary 1 and 0 for strain response and do a PCA
#=====================================================================
QIV1_res_pcadf = QIV1_responder[,c(1,22:25)]
QIV1_res_pcadf = QIV1_res_pcadf %>% mutate(across(starts_with("RES4.0."), as.integer))

mat <- QIV1_res_pcadf %>%
  column_to_rownames("SubjectID") %>%
  as.matrix()


pca.res <- prcomp(mat, scale. = TRUE)
pc.per <- round(100 * (pca.res$sdev^2 / sum(pca.res$sdev^2)), 1)
pca.df <- data.frame(pca.res$x,
                     SubjectID = rownames(mat)) %>%
  mutate(ResponsePattern = apply(mat, 1, paste, collapse = ""))

pca.plot = ggplot(pca.df, aes(x = PC1, y = PC2, color = ResponsePattern)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA on binary responses") +
  coord_fixed() +
  theme_bw()
ggsave("PCA_binaryresponse.png", plot = pca.plot, width = 16, height = 10)


#========use Foldchange or baseline to do PCA and color based on response pattern ========

mat2 = QIV1_responder_FC %>% as.matrix()
pca.res <- prcomp(mat2, scale. = TRUE)
pc.per <- round(100 * (pca.res$sdev^2 / sum(pca.res$sdev^2)), 1)
pca.df <- data.frame(pca.res$x,
                     SubjectID = rownames(mat)) %>%
mutate(ResponsePattern = apply(mat, 1, paste, collapse = ""))
pca.plot = ggplot(pca.df, aes(x = PC1, y = PC2, color = ResponsePattern)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA using foldchange") +
  coord_fixed() +
  theme_bw()
ggsave("PCA_foldchange_responsepattern.png", plot = pca.plot, width = 16, height = 10)
