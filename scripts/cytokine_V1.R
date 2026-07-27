cyt = c('IL9','IL9R','IL18','IL17A','IL17B','IL17C','IL17D','IL17F','IL17RA',
        'IL17RB','IL17RB','IL17RC','IL17RD','IL17RE','IL17REL','IL1RA','CXCL1',
        'IFNB1', 'IL1A', 'IL1B')
cyt_orr = c('IL18','CXCL1','IFNB1','IL17A','IL9','IL1A')
ensg = protein_coding[protein_coding$HGNC_symbol %in% cyt_orr,]


cyt_orr = QIV1_V1.tmm[ensg$target_id,]
rownames(cyt_orr) = ensg$HGNC_symbol

cyt_orr = QIV1_V1.num[ensg$target_id,]

annotation_col1 <- HeatmapAnnotation(
  Responder_HK = annotation_col$Responder_HongKong,
  col = list(Responder_HK = c(HR = "#D55E00", LR = "#E69F00", MR = "#56B4E9"))
)

cyt_heat = Heatmap(
  cyt_orr,
  name = "cytokine expression ",               
  col = heatmap_colors,              
  top_annotation = annotation_col1,
  show_row_names =TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  heatmap_legend_param = list(title = "Gene expression (raw)",
                              title_gp = gpar(fontsize = 18),   
                              labels_gp = gpar(fontsize = 16))
)
pdf("cyt_V1_raw_HK.pdf", width = 14, height = 15)
draw(cyt_heat)

dev.off()

QIV1_responder = read.table("/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv", sep=',',)
QIV1_responder = QIV1_responder %>% filter(SubjectID %in% QIV1_meta_filt_cond$SubjectID)
QIV1_resp_FC = QIV1_responder[,c(1,14:17)]

library(dplyr)
QIV1_resp_FC$maxFC <- apply(QIV1_resp_FC[, 2:5], 1, max)

QIV1_resp_FC = t(QIV1_resp_FC)
colnames(QIV1_resp_FC)= QIV1_resp_FC[1,]
df = rbind(cyt_orr,QIV1_resp_FC)
df = df[-(7),]
mode(df) <- "numeric"
df = df[-(7:10),]
genes  <- c("IL1A", "CXCL1", "IL9", "IL17A", "IFNB1", "IL18")
titres <- c("FC.A.HongKong", "FC.A.Victoria", "FC.B.Phuket", "FC.B.Washington")

titres = c('maxFC')
library(Hmisc)
full_spear <- rcorr(t(df), type = "spearman")
r_mat_spear <- full_spear$r[genes, titres]
p_mat_spear <- full_spear$P[genes, titres]

p_adj <- matrix(p.adjust(as.vector(p_mat_spear), method = "BH"),
                nrow = nrow(p_mat_spear), dimnames = dimnames(p_mat_spear))


#within each responder group checking the corr of cytokines to maxFC values 
group= QIV1_meta_filt_cond$ResponderGroup
strata <- levels(factor(group))

results <- data.frame()

for (g in genes) {
  for (t in titres) {
    for (s in strata) {
      idx <- which(group == s)
      x <- as.numeric(df[g, idx])
      y <- as.numeric(df[t, idx])
      mask <- !is.na(x) & !is.na(y)
      
      if (sum(mask) >= 5) {
        ct <- cor.test(x[mask], y[mask], method = "spearman")
        results <- rbind(results, data.frame(
          gene = g, titre = t, stratum = s,
          rho = ct$estimate, p = ct$p.value, n = sum(mask)
        ))
      }
    }
  }
}

# BH FDR correction across all gene-titre-stratum tests
results$fdr <- p.adjust(results$p, method = "BH")

library(tidyr)

wide_rho <- results %>% select(gene, titre, stratum, rho) %>%
  pivot_wider(names_from = stratum, values_from = rho, names_prefix = "rho_")

wide_p <- results %>% select(gene, titre, stratum, fdr) %>%
  pivot_wider(names_from = stratum, values_from = fdr, names_prefix = "fdr_")

wide <- merge(wide_rho, wide_p, by = c("gene", "titre"))
