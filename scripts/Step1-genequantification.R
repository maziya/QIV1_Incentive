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

write.csv(QIVI_kallisto_counts,'QIV1_kallisto_counts.csv', quote = FALSE)

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

write.csv(QIVI_salmon_counts,'QIV1_salmon_counts.csv', quote = FALSE)
