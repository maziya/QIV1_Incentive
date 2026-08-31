library(mixOmics)
library(tidyverse)

prepare_diablo_data <- function(
    ssgsea_file, 
    meta_data_path, 
    timepoint, 
    olink_data_path, 
    proteomics_data_path, 
    olink_mapping_path, 
    proteomics_mapping_path
) {
  
  # ==========================================
  # 1. Read and Process Transcriptomics (BTM)
  # ==========================================
  ssgsea_scores <- read.csv(ssgsea_file, row.names = 1)
  olink_data       <- read.csv(olink_data_path)
  proteomics_data  <- read.csv(proteomics_data_path)
  olink_mapping    <- read.csv(olink_mapping_path)
  proteomics_mapping <- read.csv(proteomics_mapping_path)
  meta_data = read.csv(meta_data_path)
  meta_data <- meta_data %>% 
    filter(Visit == timepoint)
  
  # Variance filtering (Keep top 80%)
  btm_variance <- apply(ssgsea_scores, 1, var)
  variance_threshold <- quantile(btm_variance, 0.2)
  ssgsea_scores_fil <- ssgsea_scores[btm_variance > variance_threshold, ]
  
  # Reshape and filter for valid metadata samples
  BTM_wide <- as.data.frame(ssgsea_scores_fil) %>%
    rownames_to_column(var = "feature") %>%
    pivot_longer(cols = -feature, names_to = "sample", values_to = "value") %>%
    mutate(view = "BTM") %>%
    filter(sample %in% meta_data$SubjectIDNew) %>%
    dplyr::select(-view) %>% 
    pivot_wider(names_from = feature, values_from = value) %>%
    column_to_rownames(var = "sample") %>%
    as.matrix()
  
  # ==========================================
  # 2. Initial Filtering of Proteomics and Olink
  # ==========================================
  valid_samples <- rownames(BTM_wide)
  
  # Filter using the dynamic timepoint string (e.g., "V1" or "V2")
  Olink_filt <- olink_data %>% dplyr::select(Protein, ends_with(timepoint))
  Olink_filt <- Olink_filt[, c("Protein", intersect(colnames(Olink_filt), valid_samples))]
  
  proteomics_filt <- proteomics_data %>% dplyr::select(Protein, ends_with(timepoint))
  proteomics_filt <- proteomics_filt[, c("Protein", intersect(colnames(proteomics_filt), valid_samples))]
  
  # ==========================================
  # 3.Map Olink IDs & Transpose
  # ==========================================
  idx <- match(Olink_filt$Protein, olink_mapping$uniprotid)
  Olink_filt <- t(column_to_rownames(Olink_filt, var = 'Protein'))
  
  olink_names <- olink_mapping$HGNC_symbol[idx]
  olink_names[is.na(olink_names)] <- colnames(Olink_filt)[is.na(olink_names)]
  colnames(Olink_filt) <- make.unique(olink_names)
  
  # ==========================================
  # 4.  Map Proteomics IDs & Transpose
  # ==========================================
  proteomics_filt <- t(column_to_rownames(proteomics_filt, var = 'Protein'))
  
  # Drop proteins not found in mapping dict
  keep <- colnames(proteomics_filt) %in% proteomics_mapping$uniprotid
  proteomics_filt <- proteomics_filt[, keep]
  
  idx2 <- match(colnames(proteomics_filt), proteomics_mapping$uniprotid)
  prot_names <- proteomics_mapping$HGNC_symbol[idx2]
  prot_names[is.na(prot_names)] <- colnames(proteomics_filt)[is.na(prot_names)]
  colnames(proteomics_filt) <- make.unique(prot_names)
  
  # ==========================================
  # 5. Final Alignment and Output
  # ==========================================
  proteomics_filt <- proteomics_filt[rownames(BTM_wide), ]
  Olink_filt <- Olink_filt[rownames(BTM_wide), ]
  
  X_list <- list(
    Transcriptomics = BTM_wide, 
    Proteomics = proteomics_filt, 
    Olink = Olink_filt
  )

  return(X_list)
}

QIV1_meta_filt_cond_V1 <- read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/QIV1_metadata_processed.csv") %>%
  filter(Visit == "V1")

set.seed(42)
X <- prepare_diablo_data(
  ssgsea_file = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/ssgsea_salmon/ssgsea_scores_V1.csv",
  meta_data_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/QIV1_metadata_processed.csv",
  timepoint = "V1",
  olink_data_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/olink.csv",
  proteomics_data_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/proteomics_log.csv",
  olink_mapping_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/olink_mapping.csv",
  proteomics_mapping_path = "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/untargeteduniprot_ENSG_HGNC.csv"
)

library(mixOmics)
library(BiocParallel)

run_strain_diablo <- function(X_input, meta_data, strain_name, timepoint) {
  
  cat("\n=======================================================\n")
  cat("Running DIABLO for:", strain_name, "| Timepoint:", timepoint, "\n")
  cat("=======================================================\n")
  
  # 1. Setup Y and Subset for HR/LR
  Y <- factor(meta_data[[strain_name]], levels = c("LR", "MR", "HR"))
  keep <- Y %in% c("HR", "LR")
  
  X_bin <- lapply(X_input, function(m) m[keep, ])
  Y_bin <- factor(as.character(Y[keep]), levels = c("LR", "HR")) 
  
  X_bin <- lapply(X_bin, function(block) {
    variances <- apply(block, 2, var, na.rm = TRUE)
    block[, variances > 1e-5, drop = FALSE]
  })
  # 2. Design Matrix & Tuning Grid
  design <- matrix(0.1, ncol = length(X_bin), nrow = length(X_bin), 
                   dimnames = list(names(X_bin), names(X_bin)))
  diag(design) <- 0
  
  test.keepX <- list(
    Transcriptomics = c(5, 10, 15, 20, 25, 30),
    Proteomics      = c(5, 10, 15, 20, 25, 30),
    Olink           = c(5, 10, 15, 20, 25, 30)
  )
  
  # 3. Tuning
  tune.diablo <- tune.block.splsda(
    X          = X_bin,
    Y          = Y_bin,
    ncomp      = 3,
    design     = design,
    scale      = TRUE,
    validation = 'loo',        
    test.keepX = test.keepX,
    measure    = "BER",
    near.zero.var = TRUE,
    BPPARAM    = MulticoreParam(workers = parallel::detectCores() - 2)
  )
  
  # Extract optimal parameters for the first 2 components
  ncomp_opt <- 2
  keepX_opt <- lapply(tune.diablo$choice.keepX, function(x) x[1:ncomp_opt])
  
  # 4. Train Final Model
  diablo.model <- block.splsda(
    X      = X_bin,
    Y      = Y_bin,
    ncomp  = ncomp_opt,
    keepX  = keepX_opt,
    design = design,
    scale  = TRUE,
    near.zero.var = TRUE
  )
  
  # 5. Generate and Save Plots Dynamically
  prefix <- paste0("diablo_", timepoint, "_", strain_name)
  color_map <- c("LR" = "steelblue", "HR" = "tomato")
  
  png(paste0(prefix, "_plotIndiv.png"), width = 10, height = 8, units = "in", res = 300)
  plotIndiv(diablo.model, ind.names = FALSE, legend = TRUE, col.per.group = color_map)
  dev.off()
  
  png(paste0(prefix, "_plotLoadings_comp1.png"), width = 15, height = 8, units = "in", res = 300)
  plotLoadings(diablo.model, comp = 1, contrib = 'max')
  dev.off()
  
  png(paste0(prefix, "_circosPlot.png"), width = 10, height = 10, units = "in", res = 300)
  circosPlot(diablo.model, cutoff = 0.7)
  dev.off()
  
  png(paste0(prefix, "_cimDiablo.png"), width = 26, height = 20, units = "in", res = 300)
  cimDiablo(diablo.model, comp = 1, margins = c(20, 35), col.names = TRUE, color.Y = color_map)
  dev.off()
  
  # 6. Performance Evaluation
  perf.diablo <- perf(
    diablo.model, 
    validation = 'Mfold',
    folds = 4,
    nrepeat = 10, 
    dist = 'centroids.dist'
  ) 

  return(list(
    model = diablo.model,
    tune = tune.diablo,
    perf = perf.diablo,
    error_rates = perf.diablo$MajorityVote.error.rate
  ))
}

strains <- c("HongKong", "Victoria", "Phuket", "Washington")

results_V1 <- lapply(strains, function(strain) {
  run_strain_diablo(
    X_input = X, 
    meta_data = QIV1_meta_filt_cond_V1,
    strain_name = strain,
    timepoint = "V1"
  )
})
names(results_V1) <- strains

# Access specific model results later (e.g., to check HongKong error rates)
print(results_V1$HongKong$error_rates)