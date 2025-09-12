QIV1_V1 <- QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V1")]
QIV1_V2 <-QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V2")]

colnames(QIV1_V1) <- gsub("V1$", "", colnames(QIV1_V1))
colnames(QIV1_V2) <- gsub("V2$", "", colnames(QIV1_V2))

adjMFC = data.frame(SubjectID=rownames(TR),adjMFC=TR$VT_decor)

adjMFC = as.data.frame(adjMFC)
adjMFC = as.data.frame(column_to_rownames(adjMFC, var="SubjectID"))

adjMFC = t(adjMFC)
common_cols = intersect(colnames(QIV1_V2), colnames(adjMFC))
adjMFC = as.matrix(adjMFC[, common_cols])
rownames(adjMFC) = "adjMFC"

QIV1_V2.adjMFC = rbind(QIV1_V2, adjMFC = adjMFC)

gene_rows = 1:(nrow(QIV1_V2.adjMFC) - 1)
adjMFC.row = nrow(QIV1_V2.adjMFC)

cor_vec = numeric(length(gene_rows))
pval_vec = numeric(length(gene_rows))
  
for (i in seq_along(gene_rows)) {
    gene_row = gene_rows[i]
    test = cor.test(as.numeric(QIV1_V2.adjMFC[gene_row, ]),as.numeric(QIV1_V2.adjMFC[adjMFC.row, ]),method ="spearman", exact = FALSE)
    cor_vec[i] = test$estimate
    pval_vec[i] = test$p.value
  }

padj_vec = p.adjust(pval_vec, method = "BH")
cor_results= data.frame(
    gene = rownames(QIV1_V2.adjMFC)[gene_rows],
    spearman_rho = cor_vec,
    p.value = pval_vec,
    adj.p.value = padj_vec)


#extract those genes which have a pval < 0.05

p_threshold = 0.05
df_sig = cor_results %>% dplyr::filter(cor_results$p.value < 0.05)

expr_long <- QIV1_V2[df_sig$gene, ] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "SubjectID", values_to = "expression") %>%
  left_join(
    data.frame(SubjectID = colnames(QIV1_V2), adjMFC = as.numeric(adjMFC)),
    by = "SubjectID"
  )
colnames(expr_long)[1] = "target_id"
expr_long = expr_long %>%
  left_join(protein_coding, by = "target_id")

#scatter plot adjMFC and the gene expression 
plot = ggplot(expr_long, aes(x = adjMFC, y = expression)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~HGNC_symbol, scales = "free_y") +
  xlab("adjMFC") +
  ylab("Gene expression") +
  ggtitle("Gene expression vs adjMFC (significant genes, p < 0.05)") +
  theme_bw() +
  theme(strip.text = element_text(size = 4),
    axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("adjMFC_correlationtogenes.png", plot = plot, width = 20, height = 15)
write.csv(df_sig, "genescorrelated_withadjMFC.csv", quote = FALSE, row.names = FALSE)
