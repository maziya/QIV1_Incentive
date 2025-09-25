QIV1_meta_filt$Responder_HongKong = "LR"
QIV1_meta_filt$Responder_HongKong[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$HongKong == "HR"]] = "HR"
QIV1_meta_filt$Responder_HongKong[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$HongKong == "MR"]] = "MR"
QIV1_meta_filt$Responder_HongKong = factor(QIV1_meta_filt$Responder_HongKong, levels = c("LR", "HR", "MR"))

QIV1_meta_filt$Responder_Victoria = "LR"
QIV1_meta_filt$Responder_Victoria[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Victoria == "HR"]] = "HR"
QIV1_meta_filt$Responder_Victoria[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Victoria == "MR"]] = "MR"
QIV1_meta_filt$Responder_Victoria = factor(QIV1_meta_filt$Responder_Victoria, levels = c("LR", "HR", "MR"))

QIV1_meta_filt$Responder_Phuket = "LR"
QIV1_meta_filt$Responder_Phuket[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Phuket == "HR"]] = "HR"
QIV1_meta_filt$Responder_Phuket[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Phuket == "MR"]] = "MR"
QIV1_meta_filt$Responder_Phuket = factor(QIV1_meta_filt$Responder_Phuket, levels = c("LR", "HR", "MR"))

QIV1_meta_filt$Responder_Washington = "LR"
QIV1_meta_filt$Responder_Washington[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Washington == "HR"]] = "HR"
QIV1_meta_filt$Responder_Washington[QIV1_meta_filt$SubjectID %in% ResponderGroups$SubjectID[ResponderGroups$Washington == "MR"]] = "MR"
QIV1_meta_filt$Responder_Washington = factor(QIV1_meta_filt$Responder_Washington, levels = c("LR", "HR", "MR"))


strains <- c("HongKong", "Victoria", "Phuket", "Washington")
results_list <- list()

for (strain in strains) {
  response_col <- paste0("Responder_", strain)
  interaction_col <- paste0("Status_Visit_", strain)
  
  QIV1_meta_filt[[response_col]] <- factor(QIV1_meta_filt[[response_col]], levels = c("LR", "HR", "MR"))
  QIV1_meta_filt$Visit <- factor(QIV1_meta_filt$Visit, levels = c("V1", "V2"))
  QIV1_meta_filt[[interaction_col]] <- interaction(QIV1_meta_filt[[response_col]], QIV1_meta_filt$Visit, drop = TRUE)
  
  #=====================================
  #Responder vs NonResponder both visits
  #=====================================
  
  design1 <- model.matrix(~0 + QIV1_meta_filt[[response_col]] + SEX + AGE + Visit, data = QIV1_meta_filt)
  colnames(design1)[1:3] <- c("LR", "HR", "MR")
  
  contrast1 <- makeContrasts(HRvsLR = HR-LR, MRvsLR = MR-LR, HRvsMR = HR-MR, levels = design1)
  
  DGE <- DGEList(QIV1_adjusted)
  DGE <- calcNormFactors(DGE, method = "TMM")
  print(identical(QIV1_meta_filt$SubjectIDNew,colnames(QIV1_adjusted)))
  v1 <- voom(DGE, design1, normalize.method = "quantile", plot = FALSE)
  corfit1 <- duplicateCorrelation(v1, design1, block = QIV1_meta_filt$SubjectID)
  v1 <- voom(DGE, design1, normalize.method = "quantile", block = QIV1_meta_filt$SubjectID, correlation = corfit1$consensus.correlation)
  
  fit1 <- lmFit(v1, design1, block = QIV1_meta_filt$SubjectID, correlation = corfit1$consensus.correlation)
  fit1 <- contrasts.fit(fit1, contrast1)
  fit1 <- eBayes(fit1, trend = TRUE)
  
  results_list[[paste0(strain, "_HRvsLR")]] <- topTable(fit1, coef = 1, number = Inf)
  results_list[[paste0(strain, "_MRvsLR")]] <- topTable(fit1, coef = 2, number = Inf)
  results_list[[paste0(strain, "_HRvsMR")]] <- topTable(fit1, coef = 3, number = Inf)
  
  
  #===========================================
  #V2 vs V1 in Res & NonRes and Res separately
  #============================================
  design2 <- model.matrix(~0 + QIV1_meta_filt[[interaction_col]] + SEX + AGE, data = QIV1_meta_filt)
  colnames(design2)[1:6] = c("LR.V1","HR.V1", "MR.V1", "LR.V2","HR.V2","MR.V2")
  
  contrast2 <- makeContrasts(
    V2vsV1_HR = HR.V2 - HR.V1,
    V2vsV1_LR = LR.V2 - LR.V1,
    V2vsV1_MR = MR.V2 - MR.V1,
    HRV2vsLRV2 = HR.V2 - LR.V2,
    HRV2vsMRV2 = HR.V2 - MR.V2,
    MRV2vsLRV2 = MR.V2 - LR.V2,
    levels = design2)
  
  v2 <- voom(DGE, design2, normalize.method = "quantile", plot = FALSE)
  corfit2 <- duplicateCorrelation(v2, design2, block = QIV1_meta_filt$SubjectID)
  v2 <- voom(DGE, design2, normalize.method = "quantile", block = QIV1_meta_filt$SubjectID, correlation = corfit2$consensus.correlation)
  
  fit2 <- lmFit(v2, design2, block = QIV1_meta_filt$SubjectID, correlation = corfit2$consensus.correlation)
  fit2 <- contrasts.fit(fit2, contrast2)
  fit2 <- eBayes(fit2, trend = TRUE)
  
  results_list[[paste0(strain, "_V2vsV1_HR")]] <- topTable(fit2, coef = 1, number = Inf)
  results_list[[paste0(strain, "_V2vsV1_LR")]] <- topTable(fit2, coef = 2, number = Inf)
  results_list[[paste0(strain, "_V2vsV1_MR")]] <- topTable(fit2, coef = 3, number = Inf)
  results_list[[paste0(strain, "_HRV2vsLRV2")]] <- topTable(fit2, coef = 4, number = Inf)
  results_list[[paste0(strain, "_HRV2vsMRV2")]] <- topTable(fit2, coef = 5, number = Inf)
  results_list[[paste0(strain, "_MRV2vsLRV2")]] <- topTable(fit2, coef = 6, number = Inf)
}
