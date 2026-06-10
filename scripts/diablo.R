library(mixOmics)

btm_fil = column_to_rownames(btm_fil, var = "BTM")
btm_fil = t(btm_fil)
proteomics_V1  = column_to_rownames(proteomics_V1, var = "Protein")
proteomics_V1  =t(proteomics_V1)
proteomics_V1 = proteomics_V1[rownames(proteomics_V1) %in% rownames(btm_fil),]

Olink_V1 = column_to_rownames(Olink_V1, var = "Protein")
Olink_V1 = t(Olink_V1)
Olink_V1 = Olink_V1[rownames(Olink_V1) %in% rownames(btm_fil),]


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

tune.diablo <- tune.block.splsda(X = X, Y = Y, ncomp = 5, design = design,
                                 validation = 'Mfold', folds = 5, nrepeat = 50,test.keepX = test.keepX,
                                 measure    = "BER")

diablo.model <- block.splsda(X = X, Y = Y, ncomp = 2, keepX =
                               list(Transcriptomics = 10, Proteomics = 10, Olink = 10), design =
                               design)

plotIndiv(diablo.model, ind.names = FALSE, legend = TRUE)
plotLoadings(diablo.model, comp = 1, contrib = 'max')
circosPlot(diablo.model, cutoff = 0.7)