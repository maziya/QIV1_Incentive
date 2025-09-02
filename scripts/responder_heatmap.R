library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)


QIV1_responder = read.table("/home/maziya/INCENTIVE/Codes/2024-05-10_Responder-Non-Responder_v6.1/QIV1-Responder.csv", sep=',',
                            header = TRUE)
responder = QIV1_responder[,c(1,14,15,16,17,27)]
responder = tibble::column_to_rownames(responder,var = "SubjectID")
mat = as.matrix(responder[,1:4])
annotation_row = data.frame(TotResp4.0 = as.factor(responder$TotResp4.0))
rownames(annotation_row) = rownames(responder)

ann_colors = list(
  TotResp4.0 = c(
    "0" = "white","1" = "yellow",
    "2" = "orange","3" = "pink","4" = "purple"))
png("heatmap_serologicalHI_foldchange.png",  width = 1500, height = 3000, res = 300)
sero_map = pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         fontsize_col = 10,
         fontsize_row = 8,
         cellwidth = 25,
         cellheight = 5,show_rownames = FALSE,
         color = colorRampPalette(c("lightblue", "blue", "darkblue"))(100),
         main = "Serological HI fold change")
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

p = ggplot(df_long, aes(x = BL, y = Ratio)) +
  geom_jitter(size = 0.5, width = 0.04, height = 0.04, color = "steelblue") +
  geom_smooth(method = "lm", color = "black", fill = "grey50", alpha = 0.05) + 
  facet_wrap(~Strain, scales = "free") +
  xlab("Log Baseline (BL)") +
  ylab("Log of PostVac / BL ratio") + ggtitle("HAI Assay Response to Each Viral Strain")+
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(size = 10), 
      panel.background = element_rect(fill = "white"),  
      axis.text.x = element_text(hjust = 1, size = 10)) 

ggsave("HAI_response_baseline_vs_logfoldchange.png", plot = p, width = 10, height = 6, dpi = 300)
