
library(readxl)
QIV1_elisa = read_excel("/home/maziya/INCENTIVE/Codes/QIV-1_THSTI.xlsx",sheet = 1,skip = 1, col_names = TRUE)

colnames_1 = c("SubjectID","Age","Gender","QIVdoses","D0_Victoria",
             "D28_Victoria","Fold-increase_VC","Responder",
             "D0_HongKong","D28_HongKong","Fold-increase_HK","Responder",
             "D0_Washington","D28_Washington","Fold-increase_WH","Responder",
             "D0_Phuket","D28_Phuket","Fold-increase_PH","Responder")

colnames(QIV1_elisa) = colnames_1
QIV1_elisa = QIV1_elisa[,-c(2,3,4,8,12,16,20)]

QIV1_elisa = QIV1_elisa %>% filter(SubjectID %in% QIV1_meta_filt_cond_V1$SubjectID)


QIV1_MN = read_excel("/home/maziya/INCENTIVE/Codes/QIV-1_THSTI.xlsx",sheet = 3,skip = 1, col_names = TRUE)
colnames(QIV1_MN) = colnames_1
QIV1_MN = QIV1_MN[,-c(2,3,4,8,12,16,20)]
QIV1_MN = QIV1_MN %>% filter(SubjectID %in% QIV1_meta_filt_cond_V1$SubjectID) %>%
  filter(if_all(everything(), ~ !is.na(.)))


subject_elisa_MN = intersect(QIV1_MN$SubjectID, QIV1_elisa$SubjectID)

QIV1_elisa = QIV1_elisa %>% filter(SubjectID %in% QIV1_MN$SubjectID)

QIV1_meta_filt_cond_V1 = QIV1_meta_filt_cond %>% filter (Visit == "V1")
QIV1_meta_filt_cond_V1 = QIV1_meta_filt_cond_V1 %>% filter (SubjectID %in% QIV1_elisa$SubjectID)
QIV1_meta_filt_cond_V1 = QIV1_meta_filt_cond_V1[,c(1,3,4,7,8,9,10,16,17)]

QIV1_MN = QIV1_MN %>% 
  rowwise() %>% 
  mutate(MN_maxfold = max(`Fold-increase_HK`,`Fold-increase_VC`,`Fold-increase_WH`,`Fold-increase_PH`))%>%
  ungroup()

QIV1_elisa = QIV1_elisa %>% 
  rowwise() %>% 
  mutate(elisa_maxfold = max(`Fold-increase_HK`,`Fold-increase_VC`,`Fold-increase_WH`,`Fold-increase_PH`))%>%
  ungroup()

QIV1_MN = left_join(QIV1_MN,QIV1_meta_filt_cond_V1, by = "SubjectID")

QIV1_responder = QIV1_responder %>% filter(SubjectID %in% QIV1_elisa$SubjectID)
QIV1_responder  = QIV1_responder  %>% 
  rowwise() %>% 
  mutate(HAI_maxfold = max(FC.A.HongKong,FC.A.Victoria,FC.B.Phuket,FC.B.Washington))%>%
  ungroup()

QIV1_MN = left_join(QIV1_MN,QIV1_responder, by = "SubjectID")
QIV1_MN = left_join(QIV1_MN, QIV1_elisa, by = "SubjectID")
QIV1_MN_1 = QIV1_MN[,c(1,14,15,16:22,49,62)]



library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# ---- Pivot long ----
QIV1_MN_long <- QIV1_MN_1 %>%
  pivot_longer(cols      = c(MN_maxfold, HAI_maxfold, elisa_maxfold),
               names_to  = "Assay",
               values_to = "MaxFold") %>%
  mutate(Assay = recode(Assay,
                        "MN_maxfold"    = "MN",
                        "HAI_maxfold"   = "HAI",
                        "elisa_maxfold" = "ELISA"),
         ResponderGroup = factor(ResponderGroup, levels = c("LR", "MR", "HR")),
         Assay          = factor(Assay, levels = c("MN", "HAI", "ELISA")))

assay_colors <- c("MN" = "#378ADD", "HAI" = "#C77CFF", "ELISA" = "#F8766D")

# # ---- Plot ----
# p <- ggplot(QIV1_MN_long, aes(x = Assay, y = MaxFold, fill = Assay)) +
#   geom_violin(alpha = 0.5, trim = FALSE) +
#   geom_jitter(aes(color = Assay), width = 0.1, size = 1.5, alpha = 0.6) +
#   
#   # pairwise Wilcoxon — spread brackets with step.increase
#   stat_compare_means(comparisons  = list(c("MN","HAI"), c("HAI","ELISA"), c("MN","ELISA")),
#                      method       = "wilcox.test",
#                      label        = "p.signif",
#                      step.increase = 0.12) +      
#   
#   facet_wrap(~ ResponderGroup, scales = "free_y") +
#   scale_fill_manual(values  = assay_colors) +
#   scale_color_manual(values = assay_colors) +
#   labs(title    = "Comparison of maxfold across assays by Responder Group",
#        caption  = "Kruskal-Wallis p-values — LR: p=0.014 | MR: p=0.02 | HR: p=0.0022",
#        x        = "Assay",
#        y        = "Max Fold Increase") +
#   theme_bw(base_size = 13) +
#   theme(plot.title    = element_text(face = "bold", hjust = 0.5),
#         plot.caption  = element_text(hjust = 0.5, size = 10, color = "grey40"),
#         strip.text    = element_text(face = "bold", size = 13),
#         legend.position = "none")
# 
# ggsave("assay_comparison_by_respondergroup.png", plot = p,
#        width = 12, height = 8, dpi = 600, bg = "white")


plot2 = ggplot(QIV1_MN_long,
       aes(Assay, MaxFold,
           group = SubjectID,
           color = ResponderGroup)) +
  
  geom_line(alpha = 0.15) +
  geom_point(alpha = 0.4) +
  labs(title    = "62 Subjects data from all 3 assays",
       subtitle = "Aggregate ResponderGroup used here",
       y = "Max(Foldchange (D28/D0))")+
  stat_summary(aes(group = ResponderGroup),
               fun = median,
               geom = "line",
               linewidth = 1.5) +
  
  stat_summary(aes(group = ResponderGroup),
               fun = median,
               geom = "point",
               size = 3) +
  
  facet_wrap(~ResponderGroup) +
  theme_bw()
ggsave("subject_3assay_response.png", plot = plot2,
       width = 12, height = 8, dpi = 600, bg = "white")


plot3 <- ggplot(QIV1_MN_long,
                aes(Assay, MaxFold, group = SubjectID)) +
  
  # subject trajectories
  geom_line(alpha = 0.08,
            linewidth = 0.4,
            color = "grey40") +
  
  # individual points
  geom_point(
    alpha = 0.7,
    size = 1,
    position = position_jitter(width = 0.04)) +
  
  # overall median trend
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linewidth = 1.8,
               color = "#A020F0") +
  
  stat_summary(aes(group = 1),
               fun = median,
               geom = "point",
               size =3,
               color = "#A020F0") +
  
  scale_y_log10() +
  
  labs(
    title = "Max Fold Change Across Assays for 62 Subjects",
    subtitle = "Individual subject trajectories with overall median trend",
    x = NULL,
    y = "Max Fold Change (log10 scale)") +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(linewidth = 0.8))

ggsave("subject_3assay_response_v4.png",
       plot = plot3,
       width = 9,
       height = 6,
       dpi = 600,
       bg = "white")



hk <- QIV1_HAI_MN %>%
  select(
    SubjectID,
    BL.A.HongKong,
    PostVac.A.HongKong,
    D0_HongKong,
    D28_HongKong
  ) %>%
  rename(
    HAI_BL = BL.A.HongKong,
    HAI_PostVac = PostVac.A.HongKong,
    MN_D0 = D0_HongKong,
    MN_D28 = D28_HongKong
  ) %>%
  mutate(Strain = "HongKong")


vc <- QIV1_HAI_MN %>%
  select(
    SubjectID,
    BL.A.Victoria,
    PostVac.A.Victoria,
    D0_Victoria,
    D28_Victoria
  ) %>%
  rename(
    HAI_BL = BL.A.Victoria,
    HAI_PostVac = PostVac.A.Victoria,
    MN_D0 = D0_Victoria,
    MN_D28 = D28_Victoria
  ) %>%
  mutate(Strain = "Victoria")


ph <- QIV1_HAI_MN %>%
  select(
    SubjectID,
    BL.B.Phuket,
    PostVac.B.Phuket,
    D0_Phuket,
    D28_Phuket
  ) %>%
  rename(
    HAI_BL = BL.B.Phuket,
    HAI_PostVac = PostVac.B.Phuket,
    MN_D0 = D0_Phuket,
    MN_D28 = D28_Phuket
  ) %>%
  mutate(Strain = "Phuket")

wh <- QIV1_HAI_MN %>%
  select(
    SubjectID,
    BL.B.Washington,
    PostVac.B.Washington,
    D0_Washington,
    D28_Washington
  ) %>%
  rename(
    HAI_BL = BL.B.Washington,
    HAI_PostVac = PostVac.B.Washington,
    MN_D0 = D0_Washington,
    MN_D28 = D28_Washington
  ) %>%
  mutate(Strain = "Washington")


combined <- bind_rows(hk, vc, ph, wh)

plot_df <- combined %>%
  pivot_longer(
    cols = -c(SubjectID, Strain),
    names_to = c("Assay", "Timepoint"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("BL", "D0", "PostVac", "D28")
    )
  )
plot2 = ggplot(plot_df,
               aes(Timepoint,
                   Value,
                   group = interaction(SubjectID, Assay),
                   color = Assay)) +
  
  geom_line(alpha = 0.2) +
  geom_point(size = 1.5, alpha = 0.7) +
  
  facet_wrap(~Strain) +
  
  scale_y_log10() +
  
  theme_bw(base_size = 14) +
  
  labs(
    title = "MN and HAI responses across strains",
    x = NULL,
    y = "Response / Titer (log10 scale)"
  ) +
  
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
  )
ggsave("subject_HAI_MN_response.png",
       plot = plot2,
       width = 9,
       height = 6,
       dpi = 600,
       bg = "white")

plot_df <- plot_df %>%
  mutate(
    Stage = case_when(
      Timepoint %in% c("BL", "D0") ~ "Baseline",
      Timepoint %in% c("PostVac", "D28") ~ "PostVaccination"
    )
  )


plot3 = ggplot(plot_df,
               aes(x = Assay,
                   y = Value,
                   group = interaction(SubjectID, Stage),
                   color = Stage)) +
  
  geom_line(alpha = 0.25,
            linewidth = 0.5) +
  
  geom_point(size = 2,
             alpha = 0.8) +
  
  facet_wrap(~Strain) +
  
  scale_y_log10() +
  
  theme_bw(base_size = 14) +
  
  labs(
    title = "HAI and MN ",
    x = NULL,
    y = "Response / Titer (log10 scale)"
  ) +
  
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
  )
ggsave("subject_HAI_MN_response_v2.png",
       plot = plot3,
       width = 9,
       height = 6,
       dpi = 600,
       bg = "white")
