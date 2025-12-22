#based on TGSig paper for creating a baseline signature
#Use anova to find genes that vary between subjects but do not have within subject
#variation at day0 and day3

meta = QIV1_meta_filt_cond[,1:6]

#normalize and log transform the batch corrected gene expr matrix
row_sums<- rowSums(QIV1_adjusted)
QIV1_adjusted<- QIV1_adjusted[which(row_sums != 0),]

DGEList = DGEList(QIV1_adjusted)
DGEList.norm = calcNormFactors(DGEList, method = "TMM")
cpm.norm = edgeR::cpm(DGEList.norm, log = TRUE)
expr = cpm.norm

ngenes = nrow(expr)

meta$SubjectIDNew <- factor(meta$SubjectIDNew)
meta$SubjectID <- factor(meta$SubjectID)
ss  = matrix(NA, ngenes, 2,
              dimnames = list(rownames(expr), c("ISV", "WSV")))
pv  = rep(NA_real_, ngenes)
names(pv) <- rownames(expr)

for (i in seq_len(ngenes)) {
  tmp = meta
  tmp$value = as.numeric(expr[i, as.character(tmp$SubjectIDNew)])

  #ANOVA
  fit = aov(value ~ SubjectID, data = tmp)
  tab = summary(fit)[[1]]
  #extract sum of squares
  ss_subject = tab["SubjectID", "Sum Sq"]
  ss_total = sum(tab[, "Sum Sq"], na.rm = TRUE)
  ss_resid = ss_total - ss_subject
  
  ss[i, "ISV"] = ss_subject
  ss[i, "WSV"] = ss_resid
  pv[i] <- tab["SubjectID", "Pr(>F)"]
}

#multiple testing correction
pv_bh = p.adjust(pv, method = "BH")

#normalize to fractions of total variance
ssn = ss / rowSums(ss, na.rm = TRUE)

#select those genes with high inter subject variation
#0.75 value taken from TGSig paper
high_tsm_isv = as.data.frame(ssn) %>% filter(ISV >= 0.75)