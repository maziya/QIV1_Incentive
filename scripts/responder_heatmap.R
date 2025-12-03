library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(pheatmap)


QIV1_responder = read.table("/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv", sep=',',
                            header = TRUE)
responder = QIV1_responder[,c(1,14,15,16,17,27)]
responder = tibble::column_to_rownames(responder,var = "SubjectID")
mat = as.matrix(responder[,1:4])
annotation_row = data.frame(TotResp4.0 = as.factor(responder$TotResp4.0))
rownames(annotation_row) = rownames(responder)
mat = log2(mat)
mat = t(scale(t(mat)))

ann_colors = list(
  TotResp4.0 = c(
    "0" = "white","1" = "yellow",
    "2" = "orange","3" = "pink","4" = "purple"))
png("heatmap_serologicalHI_foldchange_logscaled2.png",  width = 1500, height = 3000, res = 300)
sero_map = pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         fontsize_col = 10,
         fontsize_row = 8,
         cellwidth = 25,
         cellheight = 5,show_rownames = FALSE,
         color = colorRampPalette(c("lightblue", "darkblue"))(100),
         main = "HAI FC (log2 transformed and row scaled)")
dev.off()

#for histogram
matdf = log(mat)
matdf = as.data.frame(matdf)

df_long = matdf %>% 
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Foldchange")
png("serologicalHI_logfoldchange_histogram_2.png",  width = 1500, height = 1000, res = 300)
ggplot(df_long, aes(x = Foldchange)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 10) +
  facet_wrap(~ Sample, scales = "fixed")+ theme_minimal() + ggtitle("Log Fold change histogram by strain")+
  theme(
  plot.title = element_text(size = 8),     
  axis.title = element_text(size = 6),      
  axis.text  = element_text(size = 6),strip.text = element_text(size = 6))
dev.off()

#with density
ggplot(df_long, aes(x = Foldchange)) +
  geom_histogram(aes(y = after_stat(density)), fill = "steelblue", color = "white", bins = 10) +
  geom_density(color = "blue", linewidth = 1) +
  facet_wrap(~ Sample) +
  theme_minimal()+ggtitle("Log Fold change histogram by strain")+
  theme(
    plot.title = element_text(size = 8),     
    axis.title = element_text(size = 6),      
    axis.text  = element_text(size = 6),strip.text = element_text(size = 6))


#### inverse triangle relationship of baseline titre vs log fold change
df_long = QIV1_responder %>%
  select(SubjectID, 1:9) %>%                      
  column_to_rownames("SubjectID") %>%
  as.matrix() %>%
  log() %>%
  as.data.frame() %>%
  rownames_to_column("SubjectID") %>%
  pivot_longer(-SubjectID, names_to = "Condition", values_to = "Value") %>%
  separate(Condition, into = c("Timepoint", "Strain"), sep = "\\.", extra = "merge") %>%
  pivot_wider(names_from = Timepoint, values_from = Value) %>%
  mutate(Ratio = PostVac / BL)

p2 = ggplot(df_long, aes(x = BL, y = Ratio)) +
  geom_jitter(size = 0.5, width = 0.04, height = 0.04, color = "steelblue") +
  geom_smooth(method = "lm", color = "black", fill = "grey50", alpha = 0.05) +
  stat_cor(method = "spearman",label.x.npc = "left", label.y.npc = "top", size = 3)+
  facet_wrap(~Strain, scales = "free") +
  xlab("Log Baseline (BL)") +
  ylab("Log of PostVac / BL ratio") + ggtitle("HAI Assay Response to Each Viral Strain")+
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(size = 10), 
      panel.background = element_rect(fill = "white"),  
      axis.text.x = element_text(hjust = 1, size = 10)) 

ggsave("HAI_response_baseline_vs_logfoldchange1.png", plot = p2, width = 10, height = 6, dpi = 300)

#===================================================================
#boxplots (raincloud) for diabetic and sex against HAI response
#===================================================================

QIV1_meta_master = inner_join(QIV1_metadata_medicalcondition, QIV1_responder, by = "SubjectID")
QIV1_meta_master = QIV1_meta_master %>% filter(Visit == "V2")
QIV1_condition = QIV1_meta_master %>% select(SubjectID, SEX, AGE, condition, FC.A.HongKong, FC.A.Victoria, FC.B.Phuket, FC.B.Washington)
QIV1_condition <-QIV1_condition %>% distinct(SubjectID, .keep_all = TRUE)

QIV1_condition = QIV1_condition %>% 
  mutate(diabetic = case_when(
    condition == "Type 2 diabetes mellitus" ~ "Diabetic",
    TRUE ~ "Non_diabetic"
  ))
QIV1_condition = QIV1_condition %>% 
  mutate(hypertension = case_when(
    condition == "Hypertension" ~ "Hypertensive",
    TRUE ~ "Non_hypertensive"
  ))
QIV1_condition <-QIV1_condition %>% distinct(SubjectID, .keep_all = TRUE)
QIV1_condition = QIV1_condition %>% 
  rowwise() %>%
  mutate(maxFC = max(c_across(c(FC.A.HongKong,
                                FC.A.Victoria,
                                FC.B.Phuket,
                                FC.B.Washington)), na.rm = TRUE)) %>%
  ungroup()


QIV1_condition <- QIV1_condition %>%
  mutate(SEX= factor(SEX, levels = c("M", "F")))
comparisons <- list(levels(QIV1_condition$SEX))

library(ggdist)
library(ggpubr)

plot = ggplot(QIV1_condition, aes(x = SEX, y = maxFC, fill = SEX)) +
  stat_halfeye(
    adjust = 0.5,
    width = 0.6,
    justification = -0.2,
    .width = 0,
    point_colour = NA,
    alpha = 0.4
  ) +
  stat_dots(
    side = "left",
    adjust = 0.5,
    width = 0.6,
    justification = 1.2,
    alpha = 0.6,
    dotsize = 0.3
  ) +
  geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    alternative = "two.sided",
    label = "p.format",
    label.y = max(QIV1_condition$maxFC) * 1.1
  ) +
  labs(
    title = "HAI response", subtitle = "Wilcox test",
    x = "Gender",
    y = "Max HAI response (Fold-Change) across 4 strains"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))
