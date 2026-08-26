library(readxl)
library(dplyr)
library(tibble)
library(stringr)

# =============================
# METADATA PREPROCESSING
# =============================

# Define samples to remove for various QC reasons
exclude_samples <- c(
  'QIV1KEM085',                                 # No data available for V2
  'QIV1KEM033',                                 # Non-responder to all strains
  'QIV1KEM037',                                 # No serology data
  'QIV1KEM019', 'QIV1KEM055', 'QIV1KEM064',     # No consent
  'QIV1KEM077', 
  'QIV1KEM001', 'QIV1KEM079', 'QIV1KEM063',     # Outliers / Low read counts
  'QIV1KEM030', 'QIV1KEM040',
  'QIV1KEM023', 'QIV1KEM014', 'QIV1KEM043'      # High percentage of rRNA
)


QIV1_meta_filt <- read_excel("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data_24082026/QIV1_India_MetaData_Transcriptomics.xlsx", sheet = 1) %>%
  dplyr::select(SubjectID = 2, 3, 6, 7, 30) %>%
  mutate(SubjectIDNew = paste0(SubjectID, Visit)) %>%
  filter(!SubjectID %in% exclude_samples)

# covariates
QIV1_meta_filt_cond <- QIV1_meta_filt %>%
  mutate(
    Visit = factor(Visit, levels = c("V1", "V2")),
    AGE = as.numeric(AGE),
    SEX = factor(SEX),
    Library_Batch = factor(Library_Batch)
  )

# =============================
# COUNT MATRIX PREPROCESSING
# =============================

# Load Kallisto counts
QIV1count <- read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data_24082026/QIV1_kallisto_counts.csv", 
                      row.names = 1, check.names = FALSE)
protein_coding = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/protein_coding_ensemble_hgnclist.csv")


# Clean Ensembl IDs, merge with protein_coding, and filter duplicates
QIV1_all_mtx_filtered <- QIV1count %>%
  rownames_to_column(var = "target_id") %>%
  mutate(target_id = str_remove(target_id, "\\..*")) %>%          # Remove Ensemble version details
  inner_join(protein_coding, by = 'target_id') %>%                
  mutate(rowsums = rowSums(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>% 
  group_by(HGNC_symbol) %>%
  slice_max(order_by = rowsums, n = 1, with_ties = FALSE) %>%      # Keep only the row with the max counts per gene
  ungroup() %>%
  column_to_rownames(var = 'target_id') %>%
  dplyr::select(where(is.numeric), -rowsums)                              # Keep only numeric count columns

# Shift Visit labels (V2 -> V1, V3 -> V2)
colnames(QIV1_all_mtx_filtered) <- colnames(QIV1_all_mtx_filtered) %>%
  str_replace("V2$", "V1") %>%
  str_replace("V3$", "V2")


# Convert to matrix and subset by metadata subject IDs
QIV1_all.num <- as.matrix(QIV1_all_mtx_filtered)

# Ensure count matrix columns exactly match the metadata rows
matched_cols <- intersect(colnames(QIV1_all.num), QIV1_meta_filt_cond$SubjectIDNew)
QIV1_all.num <- QIV1_all.num[, matched_cols, drop = FALSE]

# Remove rows with a sum of zero across all samples
QIV1_all.num <- QIV1_all.num[rowSums(QIV1_all.num) > 0, ]

QIV1_all.num = QIV1_all.num[,QIV1_meta_filt_cond$SubjectIDNew] 

# Verify alignment (optional but recommended safety check)
if(!identical(colnames(QIV1_all.num), QIV1_meta_filt_cond$SubjectIDNew)) {
  warning("Column names in count matrix do not match metadata rows perfectly!")
}

# Split the final matrix into V1 and V2
QIV1_V1.num <- QIV1_all.num[, endsWith(colnames(QIV1_all.num), "V1")]
QIV1_V2.num <- QIV1_all.num[, endsWith(colnames(QIV1_all.num), "V2")]

# save matrices
write.csv(QIV1_all.num, "QIV1_V1&V2_countmatrix.csv", quote = FALSE)
write.csv(QIV1_V1.num, "QIV1_V1_countmatrix.csv", quote = FALSE)
write.csv(QIV1_V2.num, "QIV1_V2_countmatrix.csv", quote = FALSE)
write.csv(QIV1_meta_filt_cond, "QIV1_metadata_processed.csv", quote = FALSE, row.names = FALSE)
