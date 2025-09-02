
#extract HAI titre fold change for each QIV1 sample against each strain 
QIV1_responder_FC = QIV1_responder[,c(1,14,15,16,17)]
QIV1_responder_FC = column_to_rownames(QIV1_responder_FC,var = "SubjectID")
QIV1_responder_FC = t(QIV1_responder_FC)
QIV1_responder_FC = QIV1_responder_FC[,!colnames(QIV1_responder_FC) %in% c('QIV1KEM033','QIV1KEM019')]


#concatenate with the gene counts for all samples
#QIV1_V2 after removing the rowsums = 0
expr_FC = rbind(QIV1_V2.num,QIV1_responder_FC)

#compute spearman correlation coeff for each strain and each gene pair
gene_rows = 1:(nrow(expr_FC) - 4)
strain_rows = (nrow(expr_FC) - 3):nrow(expr_FC)
strain_names = rownames(expr_FC)[strain_rows]

cor_results =list()

for (strain_row in strain_rows) {
  strain_name = rownames(expr_FC)[strain_row]
  
  cor_vec = numeric(length(gene_rows))
  pval_vec = numeric(length(gene_rows))
  
  for (i in seq_along(gene_rows)) {
    gene_row = gene_rows[i]
    test = cor.test(expr_FC[gene_row, ],expr_FC[strain_row, ],method ="spearman",exact = FALSE)
    cor_vec[i] = test$estimate
    pval_vec[i] = test$p.value
  }
  padj_vec = p.adjust(pval_vec, method = "BH")
  cor_results[[strain_name]] = data.frame(
  gene = rownames(expr_FC)[gene_rows],
  spearman_rho = cor_vec,
  p.value = pval_vec,
  adj.p.value = padj_vec)
}

#extract those genes which have a pval < 0.05
sig_cor_results = list()
p_threshold = 0.05

for (strain in names(cor_results)) {
  df = cor_results[[strain]]
  df_sig = df[df$p.value < p_threshold,]  
  sig_cor_results[[strain]] = df_sig
  }

all_sig_genes = unique(unlist(lapply(sig_cor_results, function(x) x$gene)))

#rho_mat is a df with all genes that are sig based on pval 
# and their correlation coeff with each strain
rho_mat= matrix(NA, nrow = length(all_sig_genes), ncol = length(sig_cor_results))
rownames(rho_mat) =all_sig_genes
colnames(rho_mat) =names(sig_cor_results)

for (strain in names(sig_cor_results)) {
  df = sig_cor_results[[strain]]
  genes_present = df$gene[df$gene %in% rownames(rho_mat)]
  matching_rho = df$spearman_rho[match(genes_present, df$gene)]
  rho_mat[genes_present, strain] = matching_rho
}


library(pheatmap)
png("correlation_FC_geneexpr_heatmap.png", width = 1000, height = 800)
p= pheatmap(rho_mat,na_col = "grey",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         color = colorRampPalette(c("brown", "white", "purple"))(100),
         main = "Spearman Rho values for strain FC & gene expression in V2")
dev.off()

#some genes are correlated with more than 2 strains
genes_corrwithFC_sharedamongstrains = rowSums(!is.na(rho_mat)) >= 2
rho_mat_filtered = rho_mat[genes_corrwithFC_sharedamongstrains, ]

genes_shared_names = Gn$gene_name[which(Gn$target_id %in% genes_shared)]
