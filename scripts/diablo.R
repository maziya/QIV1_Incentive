library(mixOmics)

btm_fil = column_to_rownames(btm_fil, var = "BTM")
btm_fil = t(btm_fil)
proteomics_V1  = column_to_rownames(proteomics_V1, var = "Protein")
proteomics_V1  =t(proteomics_V1)
proteomics_V1 = proteomics_V1[rownames(proteomics_V1) %in% rownames(btm_fil),]

Olink_V1 = column_to_rownames(Olink_V1, var = "Protein")
Olink_V1 = t(Olink_V1)
Olink_V1 = Olink_V1[rownames(Olink_V1) %in% rownames(btm_fil),]

idx = match(colnames(Olink_V1), olink_ENSG_HGNC$uniprotid)
colnames(Olink_V1) = olink_ENSG_HGNC$HGNC_symbol[idx]

keep <- colnames(proteomics_V1) %in% untargeteduniprot_ENSG_HGNC$uniprotid
proteomics_V1 <- proteomics_V1[, keep]


idx2 = match(colnames(proteomics_V1), untargeteduniprot_ENSG_HGNC$uniprotid)
colnames(proteomics_V1) = untargeteduniprot_ENSG_HGNC$HGNC_symbol[idx2]


X <- list(Transcriptomics = btm_fil, Proteomics = proteomics_V1, Olink = Olink_V1)
Y <- factor(QIV1_meta_filt_cond_V1$Responder_HongKong, levels = c("LR", "MR", "HR"))

design <- matrix(0.1, ncol = 3, nrow = 3, dimnames = list(names(X),
                                                          names(X)))
diag(design) <- 0

test.keepX <- list(
  Transcriptomics = c(5, 10, 15, 20, 25, 30),
  Proteomics      = c(5, 10, 15, 20, 25, 30),
  Olink           = c(5, 10, 15, 20, 25, 30)
)
keep <- Y %in% c("HR", "LR")
X_bin <- lapply(X, function(m) m[keep, ])
Y_bin <- droplevels(Y[keep])  # 22 LR, 12 HR — manageable binary problem
table(Y_bin)



tune.diablo <- tune.block.splsda(
  X          = X_bin,
  Y          = Y_bin,
  ncomp      = 3,
  design     = design,
  validation = 'loo',       
  test.keepX = test.keepX,
  measure    = "BER",
  BPPARAM    = MulticoreParam(workers = parallel::detectCores() - 2))


tune.diablo$choice.ncomp
tune.diablo$choice.keepX
tune.diablo$error.rate

dev.off()
plot(tune.diablo)

ncomp_opt <- 2

keepX_opt <- list(
  Transcriptomics = c(30, 5),   # comp1 optimal, comp2 from choice.keepX
  Proteomics      = c(10, 5),
  Olink           = c(5, 20)
)

diablo.model <- block.splsda(
  X      = X_bin,
  Y      = Y_bin,
  ncomp  = 2,
  keepX  = keepX_opt,
  design = design
)


png("diablo_plotIndiv.png", width = 10, height = 8, units = "in", res = 300)
plotIndiv(diablo.model, ind.names = FALSE, legend = TRUE)
dev.off()


png("diablo_plotLoadings_comp1.png", width = 15, height = 8, units = "in", res = 300)
plotLoadings(diablo.model, comp = 1, contrib = 'max')
dev.off()


png("diablo_circosPlot.png", width = 10, height = 10, units = "in", res = 300)
circosPlot(diablo.model, cutoff = 0.7)
dev.off()
