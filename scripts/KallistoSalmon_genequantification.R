library(rhdf5)
library(tximport)

#from the gencode v46 basic annotation file extracted transcript and gene information
#ENST and ENSG ids and gene symbols
Tx_gene = read.table('/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/transcript_gene_map.txt')
colnames(Tx_gene)[1]="target_id"
colnames(Tx_gene)[2]="gene_id"
colnames(Tx_gene)[3]="gene_name"

Tx = Tx_gene[,c(1,2)]
#list of sample IDs to create full paths
sample = read.table('/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/directories.txt')
colnames(sample)[1]="sampleIDs"

#to quantify genes using kallisto outputs
path <- file.path(sample$sampleIDs, "abundance.h5")
all(file.exists(path)) #should be TRUE

Txi_gene_kallisto = tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = FALSE)
QIVI_kallisto_counts = Txi_gene_kallisto$counts
colnames(QIVI_kallisto_counts) = sample$sampleIDs

write.csv(QIVI_kallisto_counts,'QIVI_kallisto_counts.csv', quote = FALSE)

#to quantify genes using salmon outputs
path <- file.path(sample$sampleIDs, "quant.sf")
all(file.exists(path)) #should be TRUE

Txi_gene_salmon = tximport(path, 
                    type = "salmon", 
                    tx2gene = Tx, 
                    txOut = FALSE, 
                    countsFromAbundance = "lengthScaledTPM",
                    ignoreTxVersion = FALSE)

QIVI_salmon_counts = Txi_gene_salmon$counts
colnames(QIVI_salmon_counts) = sample$sampleIDs

write.csv(QIVI_salmon_counts,'QIVI_salmon_counts.csv', quote = FALSE)


###comparing Bioaster salmon count table and rerun count table with salmon v1.10.3

QIV1_bio_salmon = read.delim("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/BA088-INCENTIVE_eukaryota_count_table.txt",
                       header = TRUE,check.names = FALSE, row.names = 1)
QIV1_salmon = read.delim("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/QIVI_salmon_counts.csv",
                             header = TRUE, sep=",", check.names = FALSE, row.names = 1)
QIV1_kallisto = read.delim("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/QIVI_kallisto_counts.csv",
                           header = TRUE, sep=",", check.names = FALSE, row.names = 1)


common_genes = intersect(rownames(QIV1_salmon), rownames(QIV1_bio_salmon))
QIV1_salmon_common = QIV1_salmon[common_genes, ]
QIV1_bio_salmon_common = QIV1_bio_salmon[common_genes, ]

plot_df = data.frame(
  df1_count = as.vector(as.matrix(QIV1_salmon_common)),
  df2_count = as.vector(as.matrix(QIV1_bio_salmon_common))
)

#pseudocount 1 added for genes with zero counts and log transform counts
plot_df$df1_adj = log(plot_df$df1_count +1)
plot_df$df2_adj = log(plot_df$df2_count +1)
cor_val =cor(plot_df$df1_adj, plot_df$df2_adj, method = "pearson")

plot_df$A = (plot_df$df1_adj + plot_df$df2_adj) / 2
plot_df$M = plot_df$df1_adj - plot_df$df2_adj

library(ggplot2)
p1 = ggplot(plot_df, aes(x = df1_adj, y = df2_adj)) +
  geom_hex(bins = 80) +     
  scale_fill_viridis_c(trans = "log10")+
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  xlab("Gene Counts from Salmon v.1.10.3") +
  ylab("Gene Counts from Kallisto v0.50.1 run") +
  ggtitle(paste0("Comparison of gene counts\nPearson r = ", round(cor_val, 3)))

ggsave("scatter_counts_salmon_kallisto_comparison.png", plot = p1, width = 8, height = 8, dpi = 300)

p1 = ggplot(plot_df, aes(x = A, y = M)) +
  geom_point(alpha = 0.3, size = 1.2, color = "darkred") + 
  xlab("Average log2 expression (Salmon & BioAster_Salmon)") +
  ylab("log2(Salmon) - log2(BioAster_Salmon)") +
  ggtitle("MA Plot: Salmon vs BioAster_Salmon") +
  theme_bw()

ggsave("MAplot_salmon_BioAster_Salmon.png", plot = p1, width = 8, height = 6, dpi = 300)
