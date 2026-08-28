library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scales) 
library(tibble) 

# Load datasets
QIV1_olink       <- read_excel("/home/maziya/INCENTIVE/Proteomics/NPX_data_Indian_cohort.xlsx", sheet = 1)
olink_uniprotid  <- read_excel("/home/maziya/INCENTIVE/Proteomics/olinkID_UniprotID.xlsx", sheet = 1)
QIV1_untargeted  <- read_excel("/home/maziya/INCENTIVE/Proteomics/proteomics_count_table_Indian_cohort.xlsx", sheet = 1)

# Subjects to remove
remove_subjects <- c("QIV1KEM085", "QIV1KEM019", "QIV1KEM055", 
                     "QIV1KEM037", "QIV1KEM077", "QIV1KEM064")
shared_prots    <- c('P00749', 'P01137', 'P09603', 'P80511')

# ==============================================================================
# Clean and Format Olink Data
# ==============================================================================
QIV1_olink_clean <- QIV1_olink %>%
  rename(SubjectID = USUBJID, Visit = `ParticipantVisit/Visit`) %>%
  mutate(Visit = ifelse(Visit == "V4", "V2", Visit)) %>%
  filter(!(SubjectID %in% remove_subjects)) %>%
  select(SubjectID, Visit, SEX = `DataSets/DM/SEX`, AGE = `DataSets/DM/AGE`, 39:ncol(.))

olink_matrix <- QIV1_olink_clean %>%
  mutate(SampleID = paste0(SubjectID, Visit)) %>%
  select(-SubjectID, -Visit, -any_of(c("SEX", "AGE", "Status"))) %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  as.matrix()


olink_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein") %>%
  write.csv("olink_clean.csv", quote = FALSE, row.names = FALSE)

olink_long <- QIV1_olink_clean %>%
  pivot_longer(
    cols = -c(SubjectID, Visit, any_of(c("SEX", "AGE", "Status"))), 
    names_to = "Protein", 
    values_to = "value"
  ) %>%
  pivot_wider(names_from = Visit, values_from = value) %>%
  mutate(V1 = as.numeric(V1), V2 = as.numeric(V2)) %>%
  drop_na(V1, V2)

# ==============================================================================
# 3. Clean and Format Untargeted LC-MS Data
# ==============================================================================
# Clean column names (Extract first Uniprot ID before | or ;)
clean_cols <- sub("[|;].*", "", colnames(QIV1_untargeted))
colnames(QIV1_untargeted) <- clean_cols

QIV1_untargeted_clean <- QIV1_untargeted %>%
  select(PartnerName, 27:ncol(.)) %>%
  rename(SampleID = PartnerName) %>%
  mutate(SampleID = gsub("^PLAS1-", "", SampleID)) %>%
  mutate(
    SubjectID = sub("V[12]$", "", SampleID),
    Visit = str_extract(SampleID, "V[12]$")
  ) %>%
  filter(!(SubjectID %in% remove_subjects))

untarg_matrix <- QIV1_untargeted_clean %>%
  select(-SubjectID, -Visit) %>%
  mutate(across(-SampleID, as.numeric)) %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%
  as.matrix()

untarg_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein") %>%
  write.csv("proteomics_clean.csv", quote = FALSE, row.names = FALSE)

untargeted_matrix_log <- log2(untarg_matrix + 1)
write.csv(untargeted_matrix_log, "proteomics_log.csv", quote = FALSE)

untarg_long <- QIV1_untargeted_clean %>%
  select(-SampleID) %>% 
  mutate(across(-c(SubjectID, Visit), as.numeric)) %>% 
pivot_longer(
  cols = -c(SubjectID, Visit), 
  names_to = "Protein", 
  values_to = "value"
) %>%
  pivot_wider(names_from = Visit, values_from = value) %>%
  mutate(V1 = as.numeric(V1), V2 = as.numeric(V2)) %>%
  drop_na(V1, V2)


olink_final_matrix <- QIV1_olink_clean %>%
  mutate(SampleID = paste0(SubjectID, Visit)) %>% 
  select(-SubjectID, -Visit, -`DataSets/DM/SEX`, -`DataSets/DM/AGE`, -Status) %>%
  pivot_longer(cols = -SampleID, names_to = "Protein", values_to = "NPX") %>%
  pivot_wider(names_from = SampleID, values_from = NPX) %>%
  column_to_rownames(var = "Protein") %>%
  as.matrix()

untarg_final_matrix <- QIV1_untargeted_clean %>%
  select(-SubjectID, -Visit, -Status) %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>% 
  as.matrix()


# ==============================================================================
# 4. V1 vs V2 Visualizations
# ==============================================================================
plot_theme <- list(
  geom_point(size = 3, alpha = 0.8),
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40"),
  facet_wrap(~ Protein, scales = "free", ncol = 2),
  scale_color_manual(values = c("Responder" = "#A020F0", "NonResponder" = "#ff7f0e")),
  theme_bw(base_size = 14),
  theme(
    panel.grid = element_line(color = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
)

# Plot 1: Olink
p1 <- ggplot(olink_long, aes(x = V1, y = V2, color = Status)) +
  plot_theme +
  scale_x_continuous("NPX (V1) Before Vaccination", breaks = pretty_breaks(n = 5)) +
  scale_y_continuous("NPX (V2) After Vaccination", breaks = pretty_breaks(n = 5)) +
  labs(title = "NPX Expression (V1 vs V2)", color = "Status")
ggsave('Olink_V1_V2_sharedprotwith_untargeted.png', plot = p1, width = 10, height = 10)

# Plot 2: Untargeted LC-MS
p2 <- ggplot(untarg_long, aes(x = V1, y = V2, color = Status)) +
  plot_theme +
  scale_x_continuous("Expression (V1) Before Vaccination", breaks = pretty_breaks(n = 5)) +
  scale_y_continuous("Expression (V2) After Vaccination", breaks = pretty_breaks(n = 5)) +
  labs(title = "Expression (V1 vs V2) from Untargeted LC-MS", color = "Status")
ggsave('UntargetedProt_V1_V2_sharedprotwith_olink.png', plot = p2, width = 10, height = 10)

# ==============================================================================
# 5. Cross-Platform Correlation (Untargeted vs Olink)
# ==============================================================================
merged_df <- inner_join(
  untarg_long, olink_long,
  by = c("SubjectID", "Protein"),
  suffix = c("_untarg", "_olink")
) %>%
  mutate(
    log2_V1_untarg = log2(V1_untarg + 1),
    log2_V1_olink  = log2(V1_olink + 1),
    log2_V2_untarg = log2(V2_untarg + 1),
    log2_V2_olink  = log2(V2_olink + 1)
  )

# Plot 3: Cross-Platform V1
p3 <- ggplot(merged_df, aes(x = log2_V1_untarg, y = log2_V1_olink)) +
  geom_point(aes(color = Status_untarg), alpha = 0.7, size = 2) +
  facet_wrap(~ Protein, scales = "free") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman") +
  scale_color_manual(values = c("Responder" = "#A020F0", "NonResponder" = "#ff7f0e")) +
  labs(
    title = "Untargeted LC-MS vs Olink (V1 Baseline)",
    x = "Log2(V1 Untargeted + 1)",
    y = "Log2(V1 Olink + 1)",
    color = "Status"
  ) +
  theme_bw()
ggsave('UntargetedProt_V1_OlinkProt_V1_logscale.png', plot = p3, width = 10, height = 10)