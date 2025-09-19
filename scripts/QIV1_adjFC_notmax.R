#===============================================================================
# Modified from John Tsang et al. 2014
# Compute responder groups per viral strain (not max across strains)
#===============================================================================

#Standardization with median/MAD
standardize_titre <- function(x) {
  (x - median(x)) / mad(x, constant = 1)  
}
compute_VT <- function(df) {
  df_std = as.data.frame(lapply(df, standardize_titre))
  rownames(df_std) = rownames(df)
  mat = as.matrix(df_std)
  return(mat)  
}

#Inverse normal transform 
InvNorm <- function(x) {
  if (is.data.frame(x) || is.matrix(x)) {
    res = apply(x, 2, function(col) {
      n = sum(!is.na(col))
      ranks = rank(col, ties.method = "average", na.last = "keep")
      qnorm((ranks - 0.5) / n)
    })
    res = as.data.frame(res)
    rownames(res) <- rownames(x)
    return(res)
  }
}
#======================================================
#Standardize baseline and foldchange titres per strain
#======================================================
QIV1_responder = read.table("/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv",sep = ",", header = TRUE)

#Baseline titres(day0)
QIV1_responder_BL = QIV1_responder[, c(1:5)]
QIV1_responder_BL = column_to_rownames(QIV1_responder_BL, var = "SubjectID")

#Fold changes(day28/day0)
QIV1_responder_FC = QIV1_responder[, c(1, 14:17)]
QIV1_responder_FC = column_to_rownames(QIV1_responder_FC, var = "SubjectID")

TV = compute_VT(QIV1_responder_BL)  
TR = compute_VT(QIV1_responder_FC)   

#===============================================
#Apply inverse normal transform strain-wise
#===============================================
TR_INT = InvNorm(TR)
TR_INT = as.matrix(TR_INT)
#TR is now TR_INT
#===============================================================================
#Decorrelate each strain separately Bin subjects by baseline titre per strain
#then standardize foldchange within each bin
#===============================================================================
colnames(TR_INT) = c("HongKong","Victoria", "Phuket","Washington")
colnames(TV) = c("HongKong","Victoria", "Phuket","Washington")
 
TR_decor = TR_INT
for (strain in colnames(TR_INT)) {
  #compute quantile breaks
  qs <- quantile(TV[, strain], probs = seq(0, 1, 0.25), na.rm = TRUE)
  bins <- cut(TV[, strain], breaks = qs, include.lowest = TRUE, labels = FALSE)
  
  #initialize column in TR_decor
  TR_decor[, strain] <- NA
  
  for (b in unique(bins)) {
    idx <- which(bins == b)
    # assign directly to the data frame using row indices
    TR_decor[idx, strain] <- standardize_titre(TR_INT[idx, strain])
  }
}


#==================================
#Assign responder groups per strain
#===================================
ResponderGroups = as.data.frame(matrix(nrow = nrow(TR_decor), ncol = ncol(TR_decor)))
rownames(ResponderGroups) <- rownames(TR_decor)
colnames(ResponderGroups) <- colnames(TR_decor)

for (strain in colnames(TR_decor)) {
  #compute 20th and 80th percentiles for the strain
  q20 = quantile(TR_decor[, strain], probs = 0.2, na.rm = TRUE)
  q80 = quantile(TR_decor[, strain], probs = 0.8, na.rm = TRUE)
  
  #assign responder group per strain
  ResponderGroups[, strain] = cut(
    TR_decor[, strain],
    breaks = c(-Inf, q20, q80, Inf),
    labels = c("LR", "MR", "HR"),
    include.lowest = TRUE
  )
}
