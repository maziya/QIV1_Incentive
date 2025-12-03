h = ggplot(HR_V2vsV1_DEG_results, aes(x = P.Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  theme_bw() +
  labs(x = "P-value", y = "Count", title = "HR V2 vs V1")

ggsave("pvalhist_HR.png",plot = h)


heat_mat = QIV1_bulkcountsCIBERSORT_estimates 
annotation_colors <- list(
  Visit = c(V1 = "orange", V2 = "purple"),
  Status = c(HR = "deeppink4", NR = "blue")
)

annotation_col = QIV1_meta_filt %>%
  dplyr::select(SubjectIDNew, Visit, Status) %>%
  column_to_rownames("SubjectIDNew")
annotation_col$Status <- factor(annotation_col$Status)
annotation_col$Visit <- factor(annotation_col$Visit)
png("heatmap_CIBERSORTproportions_HRNR.png",  width = 1500, height = 1000, res = 300)
propheatmap = pheatmap(t(heat_mat),annotation_col= annotation_col,
                       annotation_colors = annotation_colors,
              cluster_rows = TRUE,
              cluster_cols = TRUE,fontsize = 6,
              fontsize_row = 5,legend = TRUE,fontsize_legend = 5,
              show_colnames = FALSE)
dev.off()



QIV1_bulkcountsCIBERSORT_estimates = read_csv("QIV1_bulkcountsCIBERSORT_estimates.csv")
QIV1_bulkcountsCIBERSORT_estimates= QIV1_bulkcountsCIBERSORT_estimates[,c(1,2,5,6,7,12,14,23)]
colnames(QIV1_bulkcountsCIBERSORT_estimates)[1] = "SubjectID"
QIV1_bulkcountsCIBERSORT_estimates$Visit = gsub(".*(V[0-9]+)$", "\\1", QIV1_bulkcountsCIBERSORT_estimates$SubjectID)
QIV1_bulkcountsCIBERSORT_estimates$SubjectID <-
  sub("(V[0-9]+)$", "", QIV1_bulkcountsCIBERSORT_estimates$SubjectID)
cibersort_respon = QIV1_bulkcountsCIBERSORT_estimates %>%
  inner_join(respondergroups_adjFC, by = "SubjectID")

cibersort_long <- cibersort_respon %>%
  pivot_longer(
    cols = `B cells naive`:`Neutrophils`,
    names_to = "CellType",
    values_to = "Proportion"
  )

cibersort_plotdata <- cibersort_long %>%
  pivot_longer(
    cols = HongKong:Washington,
    names_to = "Strain",
    values_to = "Responder"
  ) %>%
  filter(Responder %in% c("HR", "LR"))

cibersort_plotdata <- cibersort_plotdata %>%
  mutate(Cell_Visit = paste(CellType, Visit, sep = "_"))

cibersort_plotdata <- cibersort_plotdata %>%
  mutate(ResVisit = paste(Responder, Visit, sep = "_"))

fill_colors <- c(
  "HR_V1" = "#6a0dad",  # dark purple
  "HR_V2" = "#b19cd9",  # light purple
  "LR_V1" = "#8b4513",  # dark brown
  "LR_V2" = "#cd853f"   # light brown
)

fill_colors <- c(
  "V1" = "#6a0dad",  # dark purple
  "V2" = "#8b4513"  # dark brown
)

library(ggpubr)
comp_list <- combn(unique(cibersort_plotdata$ResVisit), 2, simplify = FALSE)

plot = ggplot(cibersort_plotdata, aes(x = CellType, y = Proportion, fill = ResVisit)) +
  geom_boxplot(aes(group = interaction(CellType, ResVisit)),
               position = position_dodge(width = 0.8),
               width = 0.6,size=0.2,
               alpha = 0.8) +
  facet_wrap(~ Strain, scales = "free_y") +
  stat_compare_means(aes(group = ResVisit),       
    method = "wilcox.test",   
    label = "p.signif",
    comparisons = comp_list,
    hide.ns = TRUE) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type",
    y = "Cell Proportion",
    fill = "Responder_Visit",
    title = "CIBERSORT Cell Proportions by Visit and Responder Group") +
  scale_fill_manual(values = fill_colors)

ggsave(filename = "CIBERSORT_cell_proportionsplot1.png", plot = plot, width = 12, height = 6, dpi = 300)

comparisons = list(c("V1", "V2"))

plot = ggplot(cibersort_plotdata, aes(x = CellType, y = Proportion, fill = Visit)) +
  geom_boxplot(aes(group = interaction(CellType, Visit)),
               position = position_dodge(width = 0.8),
               width = 0.6,size=0.2,
               alpha = 0.8) +
  facet_wrap(~ Strain, scales = "free_y") +
  stat_compare_means(aes(group = Visit),       
                     method = "wilcox.test",  
                     comparisons = comparisons,
                     label = "p.signif",       
                     hide.ns = TRUE) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type",
       y = "Cell Proportion",
       fill = "Visit",
       title = "CIBERSORT Cell Proportions by Visit") +
  scale_fill_manual(values = fill_colors)

ggsave(filename = "CIBERSORT_cell_proportionsplot_visit.png", plot = plot, width = 12, height = 6, dpi = 300)

