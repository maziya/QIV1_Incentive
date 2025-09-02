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
