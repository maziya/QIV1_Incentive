#===============================================================================
#from John Tsang et.al 2014 paper
#to make individual viral titres and their post‐vaccination responses comparable 
#so that the maximum is meaningful, we first standardized day 0 titre for indivi-
#dual i as TVi(0) max{TViHongKong, TViVictoria...} 
#===============================================================================
standardize_titre <- function(x) {
  (x - median(x)) / mad(x, constant = 1)  
}
compute_VT <- function(df) {
  #use median to standardize titre data
  df_std = as.data.frame(lapply(df, standardize_titre))
  rownames(df_std) = rownames(df)
  mat = as.matrix(df_std)
  
  #max value for each subject
  VT_max = apply(mat, 1, max)
  
  #strains with max value
  VT_top_strain = apply(mat, 1, function(r) {
    cols = colnames(mat)[r == max(r)]
    paste(cols, collapse = "; ")
  })
  #save 
  data.frame(
    VT_max = VT_max,
    VT_top_strain = VT_top_strain,
    row.names = rownames(df_std)
  )
}

#for each sample day 0 baseline values against 4 strains in QIV1_responder_BL
#for each sample fold change day28/day0 values against 4 strains in QIVI_responder_F
QIV1_responder = read.table("/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv", sep=',',
                            header = TRUE)
QIV1_responder_BL = QIV1_responder[,c(1:5)]
QIV1_responder_BL = column_to_rownames(QIV1_responder_BL,var="SubjectID")

#TR is computed using the fold change FCi(v)-median(FCv)/mad(FCv)
QIV1_responder_FC = QIV1_responder[,c(1,14:17)]
QIV1_responder_FC = column_to_rownames(QIV1_responder_FC,var="SubjectID")
TV = compute_VT(QIV1_responder_BL)  
TR = compute_VT(QIV1_responder_FC)  


#inverse normal transformation to TRi values 
INT <- function(x) {
  n = sum(!is.na(x))
  ranks = rank(x, ties.method = "average", na.last = "keep")
  qnorm((ranks - 0.5) / n)
}

TR_INT <- data.frame(
  VT_INT = INT(TR$VT_max),
  row.names = rownames(TR)
)

plotdf = cbind(TV,TR_INT)

#we plotted TVi(0) and TRi values against each other across subjects, and as expected, we
#saw a strong inverse correlation between them

INTplot = ggplot(plotdf,aes(x= VT_max, y=VT_INT))+geom_point(size = 2) +
  xlab("TV response at day0 standardized") +
  ylab("TR fold change standardized")+theme_bw()+theme(axis.title = element_text(size = 16))+
  theme(axis.text = element_text(size = 14))

ggsave("QIV1_day0_std_FC_Invesenormaltransformed.png", plot = INTplot, width = 16, height = 10)

#plot histogram for VT_max values
png("QIV1_hist_VTmax.png", width = 1600, height = 1000, res = 150)
hist(plotdf$VT_max, main = "Histogram of VT_max", xlab = "VT_max", col = "steelblue", border = "white")
dev.off()


#we binned the subjects based on their TVi(0) values and let TRj be the
#vector of TR values from subjects in bin j. For each individual, we next computed the decorrelated
#response TR_decor = TRi-median(TRj)/mad(TRj), where subject i belongs to bin j

TV$bin = cut(plotdf$VT_max,breaks = quantile(plotdf$VT_max, probs = seq(0, 1, 0.25), na.rm = TRUE),
             include.lowest = TRUE, labels = FALSE)

TR$bin = TV$bin[match(rownames(TR), rownames(TV))]
TR$VT_decor = NA
for (i in unique(TR$bin)) {
  idx = which(TR$bin == i)
  TR_binned = TR$VT_max[idx]    
  TR$VT_decor[idx] = standardize_titre(TR_binned)
}

#calculate percentile
q20 = quantile(TR$VT_decor, probs = 0.2, names = FALSE)
q80 = quantile(TR$VT_decor, probs = 0.8, names = FALSE)

#categorize into low mid and high responders
TR$responder_group = cut(
  TR$VT_decor,
  breaks = c(-Inf, q20, q80, Inf),
  labels = c("LR", "MR", "HR"),
  include.lowest = TRUE
)


#validation of decorrelation 

#continuous
cor_test1 = cor.test(TV$VT_max, TR$VT_decor, method = "spearman")

#Discretized considering low=1, mid=2, high=3
group_num = as.numeric(TR$responder_group)
cor_test2 = cor.test(TV$VT_max, group_num, method = "spearman")

virus_cols = colnames(QIV1_responder_BL)
pvals = sapply(virus_cols, function(v) {
  cor.test(QIV1_responder_BL[[v]], TR$VT_decor, method = "spearman")$p.value
})

kruskal_test = kruskal.test(TR$VT_decor ~ TR$VT_top_strain)

responder_groups = data.frame(
  SubjectID = rownames(TR),
  ResponderGroup = TR$responder_group)
