#V2 gene expression correlation with adjFC values for each strain in each 
#sample
#extract V2 gene exp matrix

QIV1_V1 = QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V1")]
QIV1_V2 = QIV1_adjusted[, endsWith(colnames(QIV1_adjusted), "V2")]

colnames(QIV1_V1) = gsub("V1$", "", colnames(QIV1_V1))
colnames(QIV1_V2) = gsub("V2$", "", colnames(QIV1_V2))
row_sums<- rowSums(QIV1_V2)
QIV1_V2<- QIV1_V2[which(row_sums != 0),]

#use the adjFC computed without max of all strains, instead for each strain
#an adjFC

adjFC = data.frame(adjFC=TR_decor)
adjFC = t(adjFC)
common_cols = intersect(colnames(QIV1_V2), colnames(adjFC))
adjFC = as.matrix(adjFC[, common_cols])

#bind the adjFC rows to the V2 matrix
QIV1_V2.adjFC = rbind(QIV1_V2, adjFC = adjFC)

gene_rows = 1:(nrow(QIV1_V2.adjFC) - 4)
strain_rows = (nrow(QIV1_V2.adjFC) - 3):nrow(QIV1_V2.adjFC)
strain_names = rownames(QIV1_V2.adjFC)[strain_rows]

#compute correlation between the adjFC and gene exp values
cor_results =list()

for (strain_row in strain_rows) {
  strain_name = rownames(QIV1_V2.adjFC)[strain_row]
  
  cor_vec = numeric(length(gene_rows))
  pval_vec = numeric(length(gene_rows))
  
  for (i in seq_along(gene_rows)) {
    gene_row = gene_rows[i]
    test = cor.test(QIV1_V2.adjFC[gene_row, ],QIV1_V2.adjFC[strain_row, ],method ="spearman",exact = FALSE)
    cor_vec[i] = test$estimate
    pval_vec[i] = test$p.value
  }
  padj_vec = p.adjust(pval_vec, method = "BH")
  cor_results[[strain_name]] = data.frame(
    gene = rownames(QIV1_V2.adjFC)[gene_rows],
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

#rho_mat is a df with all genes that are significant based on pval 
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
#only positively correlated genes 
rho_mat_positive <- rho_mat[apply(rho_mat[,1:4], 1, function(x) all(x[!is.na(x)] > 0)), ]

library(pheatmap)
png("correlation_adjFC_geneexpr_heatmap.png", width = 1000, height = 800)
p= pheatmap(rho_mat,na_col = "grey",
            cluster_rows = FALSE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            color = colorRampPalette(c("brown", "white", "purple"))(100),
            main = "Spearman Rho values for strain adjFC & gene expression in V2")
dev.off()

#some genes are correlated with more than 2 strains
genes_corrwithFC_sharedamongstrains = rowSums(!is.na(rho_mat)) >= 2
rho_mat_filtered = rho_mat[genes_corrwithFC_sharedamongstrains, ]
shared_ids = intersect(protein_coding$target_id, rownames(rho_mat_filtered))
genes_shared_names = protein_coding$HGNC_symbol[match(shared_ids, protein_coding$target_id)]
write.csv(genes_shared_names, "genes_corr_adjFC_across2ormore_Strains.csv", quote = FALSE, row.names = FALSE)

#get the HGNC_symbols
shared_ids = intersect(protein_coding$target_id, cor_results[["adjFC.Washington"]]$gene)
cor_results[["adjFC.Washington"]]$HGNC_symbol = protein_coding$HGNC_symbol[match(shared_ids, protein_coding$target_id)]

write.csv(cor_results[["adjFC.Washington"]],"corradjFC_genes_Washington.csv", quote = FALSE, row.names = FALSE)
