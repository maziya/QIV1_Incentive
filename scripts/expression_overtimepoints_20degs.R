top20_genes <- all_genes.day0HRLR %>%
  arrange(P.Value) %>%
  slice_head(n = 20) %>%
  pull(target_id)


QIV1_all.num = QIV1_all.num[, (colnames(QIV1_all.num) %in% QIV1_meta_filt$SubjectIDNew)]

top20.QIV1_all.num = QIV1_all.num[top20_genes,]
top20.QIV1_all.num = top20.QIV1_all.num[,!colnames(top20.QIV1_all.num) %in% c('QIV1KEM001V1','QIV1KEM001V2','QIV1KEM063V2','QIV1KEM063V1')]

dge <- DGEList(top20.QIV1_all.num)
dge <- calcNormFactors(dge, method = "TMM")
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)


library(tidyverse)
norm_counts.df = as.data.frame(norm_counts)
counts_long <- norm_counts.df %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    SubjectID = str_sub(sample, 1, 10),
    Visit = str_sub(sample, 11, 13)
  )

plot_df <- counts_long %>%
  left_join(QIV1_meta_filt, by = c("SubjectID", "Visit"))

plot2 = ggplot(plot_df, aes(x = Visit, y = expression, group = SubjectID, color = Status)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Top20(pvalue<0.05) Day0 HRvsLR gene expression over Timepoints",
       x = "Visit", y = "Expression")
ggsave('QIV1_20genes_expressionday0_visits&responderstatus.png',plot = plot2,width = 15, height=15)
