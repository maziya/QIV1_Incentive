library(MOFA2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tibble)

QIV1_meta_filt_cond_V1 <- QIV1_meta_filt_cond %>% 
  filter(Visit == "V1")

# Variance filtering
btm_variance <- apply(ssgsea_scores_V1, 1, var)
variance_threshold <- quantile(btm_variance, 0.2)
ssgsea_scores_fil <- ssgsea_scores_V1[btm_variance > variance_threshold, ]

# Reshape and format
BTM_long <- as.data.frame(ssgsea_scores_fil) %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = -feature, names_to = "sample", values_to = "value") %>%
  mutate(view = "BTM") %>%
  filter(sample %in% QIV1_meta_filt_cond_V1$SubjectIDNew)


# Proteomics
proteomics_long <- proteomics %>%
  dplyr::select(1, ends_with("V1")) %>%
  rename(feature = 1) %>% 
  pivot_longer(cols = -feature, names_to = "sample", values_to = "value") %>%
  mutate(view = "proteomics") %>%
  filter(sample %in% BTM_long$sample)

# Olink

Olink_long <- Olink %>%
  rename(feature = Protein) %>%
  dplyr::select(feature, ends_with("V1")) %>%
  pivot_longer(cols = -feature, names_to = "sample", values_to = "value") %>%
  group_by(feature) %>%
  mutate(value = as.numeric(scale(value, center = TRUE, scale = TRUE))) %>%  # z-score per protein
  ungroup() %>%
  mutate(view = "olink") %>%
  filter(sample %in% BTM_long$sample)



# Metabolomics
metabolomics_long <- metabolomics %>%
  dplyr::select(1, ends_with("V1")) %>%
  rename(feature = 1) %>%
  pivot_longer(cols = -feature, names_to = "sample", values_to = "value") %>%
  mutate(view = "metabolomics") %>%
  filter(sample %in% BTM_long$sample)



# ==========================================
# 1. MOFA Input Data
# ==========================================
# Bind views and select specific columns in the order MOFA expects
mofa_data <- bind_rows(BTM_long, Olink_long, metabolomics_long) %>%
  select(sample, feature, view, value) 

# ==========================================
# 2. Create MOFA Object & Plot Overview
# ==========================================
MOFAobject <- create_mofa(mofa_data)
print(MOFAobject)

ov <- plot_data_overview(MOFAobject)
ggsave("dataoverview_mofa.png", plot = ov, width = 6, height = 4, dpi = 300)

views <- unique(mofa_data$view)

# Plot histograms in a grid
par(mfrow = c(2, 2))
for(v in views){
  hist(
    mofa_data$value[mofa_data$view == v],
    main = paste("Histogram for", v),
    xlab = "Value",
    col = "lightblue",
    border = "white"
  )
}
par(mfrow = c(1, 1)) 

# ==========================================
# 4. Configure & Prepare MOFA Model Options
# ==========================================
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE 

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

print(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$seed <- 42

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# ==========================================
# 5. Run MOFA Training
# ==========================================
outfile <- file.path(getwd(), "QIVI_V1_BTM3Olinkmets_MOFA_salmon.hdf5")
set.seed(42)
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

# ==========================================
# 6. Extract & Custom Plot Variance Explained
# ==========================================
var_exp <- get_variance_explained(MOFAobject.trained)

df_long <- as.data.frame(var_exp$r2_per_factor[[1]]) %>%
  rownames_to_column(var = "Factor") %>%
  pivot_longer(
    cols = -Factor,
    names_to = "View",
    values_to = "VarianceExplained"
  ) %>%
  mutate(
    View = case_when(
      View == "BTM" ~ "BloodTranscriptionalModule",
      TRUE ~ View
    )
  ) %>% 
  filter(Factor %in% c("Factor1", "Factor2", "Factor3"))


factor_varplot <- ggplot(df_long, aes(x = Factor, y = VarianceExplained, fill = View)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Percentage variance explained by each view",
    x = "Factor",
    y = "Variance Explained (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave("factor_varplot_oldbtmolink.png", plot = factor_varplot, width = 6, height = 6, dpi = 600)

# ==========================================
# 7. Standard MOFA Variance Plots
# ==========================================
# Variance explained by factor
v <- plot_variance_explained(MOFAobject.trained)
ggsave("totalvariancebyfactor_mofa.png", plot = v, width = 6, height = 5, dpi = 300)
print(MOFAobject.trained@cache$variance_explained$r2_per_factor)

# Total variance explained across all factors
var_plot <- plot_variance_explained(MOFAobject.trained, plot_total = TRUE)[[2]]
ggsave("totalvariance_mofa.png", plot = var_plot, width = 5, height = 5, dpi = 300)
print(MOFAobject.trained@cache$variance_explained$r2_total)

# ==========================================
# Add metadata and plot Factors
# ==========================================
sample_metadata <- data.frame(
  sample = QIV1_meta_filt_cond_V1$SubjectIDNew,
  responder = QIV1_meta_filt_cond_V1$HongKong,
  visit = "Day0",
  age = QIV1_meta_filt_cond_V1$AGE, 
  sex = QIV1_meta_filt_cond_V1$SEX
)

samples_metadata(MOFAobject.trained) <- sample_metadata
meta <- samples_metadata(MOFAobject.trained)

factors <- MOFA2::get_factors(MOFAobject.trained, factors = "all")[[1]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample")

df <- meta %>%
  left_join(factors, by = "sample")

df2_long <- df %>%
  pivot_longer(
    cols = c("Factor1"), 
    names_to = "Factor",
    values_to = "Value"
  ) %>%
  mutate(responder = factor(responder, levels = c("LR", "MR", "HR")))
comparisons <- list(c("LR", "MR"), c("MR", "HR"), c("LR", "HR"))


p <- ggplot(df2_long, aes(x = responder, y = Value, fill = responder)) +
  geom_boxplot(alpha = 0.3, notch = FALSE) +
  geom_jitter(width = 0.1) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    size = 4,
    step.increase = 0.12
  ) +
  # facet_wrap(~Factor, scales = "free_y") +
  theme_classic() +
  labs(title = "Day0",
    x = "Responder aggregate",
    y = "MOFA factor 1 value"
  )+
  theme_classic()

ggsave("Factor1_agg_oldbtm_v1.png", p, width = 6, height = 5, dpi = 300)


#check correlation of factors
#they should not be correlated that would mean it is a poor model fit

png("factor_correlation.png", width = 6, height = 6, units = "in", res = 300)
plot_factor_cor(MOFAobject.trained)
dev.off()

p_cov_cor <- correlate_factors_with_covariates(
  MOFAobject.trained, 
  covariates = c('age', 'sex'),
  plot = "log_pval"
)
ggsave("factors_corr_covariates.png", plot = p_cov_cor, width = 6, height = 5)


#The sign of the weights indicates the direction of the effect: a positive weights 
#indicates that the feature has higher levels in the samples with positive factor 
#values, and vice-versa.

# Top drivers for Factor 1 in the BTM view
q_btm <- plot_weights(
  MOFAobject.trained, 
  view = "BTM", 
  factor = 1,
  nfeatures = 10, 
  scale = TRUE, 
  abs = FALSE, 
  text_size = 2.5
)
ggsave("factor1_BTM_weights.png", plot = q_btm, width = 6, height = 5)



#which genes drive the factor based on weights and get an association of features
#and factors in a scatterplot
q_scatter <- plot_data_scatter(
  MOFAobject.trained,
  view = "olink",
  factor = 1,
  features = 12,
  color_by = "responder",
  add_lm = TRUE,
  dot_size = 1
)
ggsave("factor1_olink_weightsscatter.png", plot = q_scatter, width = 10, height = 10)



clusters <- cluster_samples(MOFAobject.trained, k = 2, factors = 1)
df2 <- as.data.frame(get_factors(MOFAobject.trained, factors = "all")[[1]]) %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_metadata, by = "sample") %>%
  mutate(Cluster = factor(clusters$cluster))

c1 <- ggplot(df2, aes(x = Factor1, y = Factor2, color = responder, shape = Cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 12) +
  labs(
    title = "MOFA Sample Clusters (k = 2)",
    x = "Factor 1",
    y = "Factor 2",
    shape = "Cluster",
    color = "Responder"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggsave("factor12_clustersamples_mofa.png", plot = c1, width = 8, height = 6)


set.seed(12345)
MOFAobject_umap <- run_umap(MOFAobject.trained)

p_umap <- plot_dimred(MOFAobject_umap, method = "UMAP", color_by = "Factor1", dot_size = 3)
ggsave("mofa_umap_factor1.png", plot = p_umap, width = 6, height = 5)

# Factor 1 vs 2 Scatter 
category.colors <- c("HR" = "darkturquoise", "LR" = "magenta3", "MR" = "#FFCC99")

p_factors <- plot_factors(MOFAobject.trained, factors = c(1,2), color_by = "responder", dot_size = 4) + 
  stat_ellipse(aes(color = color_by), geom = "polygon", alpha = 0.25) +
  scale_fill_manual(values = category.colors) +
  scale_color_manual(values = category.colors)
ggsave("factor1_2_scatter_ellipses.png", plot = p_factors, width = 7, height = 5)




weights <- MOFA2::get_weights(MOFAobject.trained, as.data.frame = TRUE)
weights_factor1 <- weights %>% filter(factor == "Factor1")

weights_other <- weights_factor1 %>% filter(view != "olink")
weights_olink <- weights_factor1 %>% 
  filter(view == "olink") %>%
  rename(uniprotid = feature) %>%
  left_join(olink_ENSG_HGNC, by = "uniprotid") %>%
  mutate(feature = coalesce(HGNC_symbol, uniprotid)) %>% 
  select(feature, factor, value, view)

weights_factor1_clean <- bind_rows(weights_other, weights_olink)


# Extract Top 10 Positive and Negative Drivers per view
top_features <- weights_factor1_clean %>%
  group_by(view) %>%
  arrange(desc(value)) %>%
  slice(c(1:10, (n() - 9):n())) %>% 
  ungroup()

q3 <- ggplot(top_features, aes(x = reorder(feature, value), y = value, fill = value > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~view, scales = "free_y") +
  scale_fill_manual(values = c("#6A3D9A", "#FF7F00"), guide = "none") +
  theme_classic() +
  labs(
    x = "",
    y = "MOFA weight",
    title = "Top drivers of Factor 1 by view Day0"
  ) +
  theme(
    text = element_text(family = "Arial"), 
    axis.text.y = element_text(size = 14),   
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(size = 14),
    strip.text = element_text(size = 16)
  )

ggsave("Factor1_top_drivers_by_view_oldbtmv1.png", plot = q3, width = 15, height = 10, dpi = 300)
