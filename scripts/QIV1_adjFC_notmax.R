#===============================================================================
#Modified from John Tsang et al. 2014
#Compute responder groups per viral strain (not max across strains)
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
#De-correlate each strain separately bin subjects by baseline titre per strain
#then standardize foldchange within each bin
#===============================================================================
colnames(TR_INT) = c("HongKong","Victoria","Phuket","Washington")
colnames(TV) = c("HongKong","Victoria","Phuket","Washington")
 
TR_decor = TR_INT
for (strain in colnames(TR_INT)) {
  
  #compute quantile breaks
  qs = quantile(TV[, strain], probs = seq(0, 1, 0.25), na.rm = TRUE)
  bins = cut(TV[, strain], breaks = qs, include.lowest = TRUE, labels = FALSE)
  
  #initialize column in TR_decor
  TR_decor[, strain] <- NA
  
  for (b in unique(bins)) {
    idx = which(bins == b)
    TR_decor[idx, strain] = standardize_titre(TR_INT[idx, strain])
  }
}

#====================================================
#Assign responder groups per strain based on quantile
#====================================================
ResponderGroups = as.data.frame(matrix(nrow = nrow(TR_decor), ncol = ncol(TR_decor)))
rownames(ResponderGroups) <- rownames(TR_decor)
colnames(ResponderGroups) <- colnames(TR_decor)

for (strain in colnames(TR_decor)) {
  #compute 20th and 80th percentiles for the strain
  q20 = quantile(TR_decor[, strain], probs = 0.2, na.rm = TRUE)
  q80 = quantile(TR_decor[, strain], probs = 0.8, na.rm = TRUE)
  # q50 = quantile(TR_decor[, strain], probs = 0.5, na.rm = TRUE)
  
  #assign responder group per strain
  ResponderGroups[, strain] = cut(
    TR_decor[, strain],
    breaks = c(-Inf, q20, q80, Inf),
    labels = c("LR", "MR", "HR"),
    # breaks = c(-Inf, q50, Inf),
    # labels = c("LR", "HR"),
    include.lowest = TRUE
  )
}


### Plot the adjFC vs standardized baseline to check if decorrelation works
###

TV_df =as.data.frame(TV) %>% 
  mutate(SubjectID = rownames(TV))
TR_df = as.data.frame(TR_decor) %>% 
  mutate(SubjectID = rownames(TR_decor))

TV_long = TV_df %>%
  pivot_longer(cols = -SubjectID,
    names_to = "Strain",
    values_to = "TV")

TR_long = TR_df %>%
  pivot_longer(cols = -SubjectID,
    names_to = "Strain",
    values_to = "TR_decor")

plot_df = left_join(TV_long, TR_long, by = c("SubjectID", "Strain"))

p = ggplot(plot_df, aes(x = TV, y = TR_decor)) +
  geom_jitter(size = 0.5, width = 0.04, height = 0.04, color = "steelblue") +
  geom_smooth(method = "lm", color = "black", fill = "grey50", alpha = 0.05) + 
  stat_cor(method = "spearman",label.x.npc = "left", label.y.npc = "top", size = 3)+
  facet_wrap(~Strain, scales = "free") + 
  labs(title = "Decorrelated HAI Response vs Baseline (Per Strain)",
    x = "Standardized Baseline titre",
    y = "adjFC") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(size = 10), 
        panel.background = element_rect(fill = "white"),  
        axis.text.x = element_text(hjust = 1, size = 10)) 

ggsave("Decorrelated_FoldChange_vs_Baseline.png", plot = p, width = 10, height = 6, dpi = 300)

# Histogram for the adjFC bins to check if top 20 and lower 20 quantiles are
# good to define HR and LR

qs_long = TR_long %>%
  group_by(Strain) %>%
  summarise(
    `20th Quantile` = quantile(TR_decor, 0.2, na.rm = TRUE),
    `80th Quantile` = quantile(TR_decor, 0.8, na.rm = TRUE),
    `50th Quantile` = quantile(TR_decor, 0.5, na.rm = TRUE),
  ) %>%
  pivot_longer(-Strain, names_to = "Quantile", values_to = "q_value")


p = ggplot(TR_long, aes(x = TR_decor)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(data = qs_long,
             aes(xintercept = q_value, color = Quantile, linetype = Quantile),
             size = 1) +
  facet_wrap(~Strain, scales = "fixed") +
  theme_classic(base_size = 12) +
  labs(x = "adjFC", y = "Count", color = "Quantiles", linetype = "Quantiles") +
  scale_color_manual(values = c("20th Quantile" = "blue",
                                "80th Quantile"  = "orange",
                                "50th Quantile"  = "purple")) +
  scale_linetype_manual(values = c("20th Quantile" = "dashed",
                                   "80th Quantile" = "dashed",
                                   "50th Quantile" = "dotted"))
ggsave("Histogram_adjFC_bins.png", p, width = 10, height = 8, dpi = 300)
