library(MOFA2)


transcriptomics_V2 <- transcriptomics[, endsWith(colnames(transcriptomics), "V2")]

Olink_V2 <- Olink %>%
            select(1, ends_with("V2"))
proteomics_V2 <- proteomics %>%
                select(1, ends_with("V2"))
metabolomics_V2 <- metabolomics %>%
  select(1, ends_with("V2"))



transcriptomics_V2 <- transcriptomics_V2 %>%
  tibble::rownames_to_column("Gene")
transcriptomics_long = transcriptomics_V2 %>%
  pivot_longer(cols = -Gene,
               names_to = "SubjectID",
               values_to = "Expression")

ssgsea_scores_fil = ssgsea_V2_BTMscores #taken from csv file of BTM scores V2
btm_variance <- apply(ssgsea_scores_fil, 1, var)
variance_threshold <- quantile(btm_variance, 0.4)
ssgsea_scores_fil <- ssgsea_scores_fil[btm_variance > variance_threshold, ]

btm_fil = as.data.frame(ssgsea_scores_fil)
btm_fil = rownames_to_column(btm_fil, var= "BTM")
BTM_long = btm_fil %>%
  pivot_longer(cols = -BTM,
               names_to = "sample",
               values_to = "value") %>%
  mutate(view = "BTM")

#WGCNA_V1 and WGCNA_V2 come from net_V1[["MEs"]] of wgcna runs
wgcna = t(WGCNA_V2)
wgcna = as.data.frame(wgcna)
wgcna_fil = rownames_to_column(wgcna, var= "Modules")
wgcna_long = wgcna_fil %>%
  pivot_longer(cols = -Modules,
               names_to = "sample",
               values_to = "value") %>%
  mutate(view = "WGCNA")

metabolomics = as.data.frame(metabolomics)
metabolomics <- metabolomics %>%
  tibble::rownames_to_column("Metabolite")

metabolomics_long = metabolomics_V2 %>%
  pivot_longer(cols = -Metabolite,
               names_to = "SubjectID",
               values_to = "Expression")

proteomics_long = proteomics_V2 %>%
  pivot_longer(cols = -Protein,
               names_to = "SubjectID",
               values_to = "Expression")

Olink_long = Olink_V2 %>%
  pivot_longer(cols = -Protein,
               names_to = "SubjectID",
               values_to = "Expression")

colnames(Olink_long)[1] = "Feature"
colnames(metabolomics_long)[1] = "Feature"
colnames(transcriptomics_long)[1] = "Feature"
colnames(proteomics_long)[1] = "Feature"


transcriptomics_long<- transcriptomics_long %>%
rename(sample = SubjectID,
value = Expression) %>%
mutate(view = "transcriptomics")

transcriptomics_long = transcriptomics_long %>% filter(sample %in% QIV1_meta_filt_cond_V2$SubjectIDNew)


proteomics_long<- proteomics_long %>%
  rename(sample = SubjectID,
  value = Expression) %>%
  mutate(view = "proteomics")

proteomics_long <- proteomics_long %>%
  filter(sample %in% BTM_long$sample)

Olink_long<- Olink_long %>%
  rename(sample = SubjectID,
         value = Expression)%>%
  mutate(view = "olink")

Olink_long <- Olink_long %>%
  filter(sample %in% BTM_long$sample)

metabolomics_long<- metabolomics_long %>%
  rename(sample = SubjectID,
         value = Expression) %>%
  mutate(view = "metabolomics") 

metabolomics_long <- metabolomics_long %>%
  filter(sample %in% BTM_long$sample)

colnames(BTM_long)[1] = "Feature"
colnames(wgcna_long)[1] = "Feature"
data = bind_rows(proteomics_long,Olink_long,metabolomics_long)
data = bind_rows(wgcna_long,Olink_long)

colnames(data)[1] ="feature"
data= data[, c("sample", "feature", "view","value")]

MOFAobject <- create_mofa(data)
print(MOFAobject)
ov = plot_data_overview(MOFAobject)
ggsave("dataoverview_mofa.png", plot = ov)

data_opts <- get_default_data_options(MOFAobject)

views <- unique(data$view)
par(mfrow=c(2,2))
for(v in views){
  hist(
    data$value[data$view == v],
    main = paste("Histogram for", v),
    xlab = "Value",
    col = "lightblue",
    border = "white"
  )
}

data_opts$scale_views = TRUE 

model_opts = get_default_model_options(MOFAobject)
model_opts$num_factors <- 7
print(model_opts)
train_opts = get_default_training_options(MOFAobject)
MOFAobject = prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


outfile = file.path(getwd(),"QIVI_V1_olinkuntargmetsMOFA.hdf5")
MOFAobject.trained = run_mofa(MOFAobject, outfile, use_basilisk=TRUE)


var_exp <- get_variance_explained(MOFAobject.trained)


df3 <- as.data.frame(var_exp$r2_per_factor[[1]])

df3$Factor <- rownames(df3)

df3_long <- df3 %>%
  pivot_longer(
    cols = -Factor,
    names_to = "View",
    values_to = "VarianceExplained"
  )
df3_long <- reshape2::melt(var_exp$r2_per_factor[[1]])

df3_long = df3_long %>%
  mutate(
    View = case_when(
      View == "BTM" ~ "BloodTranscriptionalModule",
      TRUE ~ View))
colnames(df3_long) <- c("Factor", "View", "VarianceExplained")

df3_long = df3_long %>% filter(Factor %in% c("Factor1", "Factor2","Factor3"))
factor_varplot = ggplot(df3_long,
                        aes(Factor,
                            VarianceExplained,
                            fill = View)) +
  
  geom_col(position = "dodge") +
  labs(title = "Percentage variance explained by each view")+
  coord_flip() +
  
  theme_classic()+
                theme(plot.title = element_text(size = 16),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10))

ggsave("factor_varplot.png", plot = factor_varplot, width = 6, height = 6, dpi = 600)

v = plot_variance_explained(MOFAobject.trained)
ggsave("totalvariancebyfactor_mofa.png", plot = v)
print(MOFAobject.trained@cache$variance_explained$r2_per_factor)

var_plot <- plot_variance_explained(MOFAobject.trained, plot_total = T)[[2]]
ggsave("totalvariance_mofa.png", plot = var_plot)

print(MOFAobject.trained@cache$variance_explained$r2_total)


#add metadata
sample_metadata = data.frame(sample=QIV1_meta_filt_cond1$SubjectIDNew,
                                           responder =QIV1_meta_filt_cond1$ResponderGroup,visit = QIV1_meta_filt_cond1$Visit,
                             age = QIV1_meta_filt_cond1$AGE, sex =QIV1_meta_filt_cond1$SEX)

sample_metadata = sample_metadata %>% mutate(visit = case_when(visit == "V2" ~ "Day3", visit == "V1" ~ "Day0" )) 
sample_metadata = sample_metadata %>% filter(visit == "Day0")


samples_metadata(MOFAobject.trained) = sample_metadata

# p1= plot_factor(MOFAobject.trained, 
#             factor = 3,
#             color_by = "responder",
#             dot_size = 3,
#             dodge = T,
#             legend = T,
#             add_violin = T,
#             violin_alpha = 0.25) 

meta <- samples_metadata(MOFAobject.trained)

df <- cbind(meta, factors)

comparisons <- list(
  c("LR","MR"),
  c("MR","HR"),
  c("LR","HR")
)
df_long <- pivot_longer(
  df,
  cols = c(Factor1, Factor2,Factor4),
  names_to = "Factor",
  values_to = "Value"
)
df2_long <- pivot_longer(
  df,
  cols = c(Factor2),
  names_to = "Factor",
  values_to = "Value"
)

df2_long$responder <- factor(df2_long$responder, levels = c("LR", "MR","HR"))

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
    x = "Responder Group",
    y = "MOFA factor 2 value"
  )+
  theme_classic()

ggsave("Factor2_aggregateres_v2.png", p, width = 6, height = 5, dpi = 300)


library(ggplot2)


s1 = plot_factor_cor(MOFAobject.trained)
s2 = plot_factors_vs_cov(MOFAobject.trained)

#The sign of the weights indicates the direction of the effect: a positive weights 
#indicates that the feature has higher levels in the samples with positive factor 
#values, and vice-versa.



q = plot_weights(MOFAobject.trained,view = "BTM", factors = 1, 
                 nfeatures = 10, scale = T, abs = F, text_size = 2.5)
print(q)
ggsave("factor1_BTM_weights.png",plot = q)

plot_top_weights(MOFAobject.trained,view = "BTM",
                 factor = 2,
                 nfeatures = 10,    
                 scale = T)

#which genes drive the factor based on weights and get an association of features
#and factors in a scatterplot
q1 = plot_data_scatter(MOFAobject.trained, 
                  view = "olink", 
                  factor = 2, 
                  features = 12,
                  color_by = "responder",
                  add_lm = T,
                  dot_size = 1)
ggsave("factor2_olink_weightsscatter.png",plot = q1, width = 10, height = 10)

p3 = plot_factor_cor(MOFAobject.trained) #check correlation of factors
#they should not be correlated that would mean it is a poor model fit

p4 = correlate_factors_with_covariates(MOFAobject.trained, covariates = c('AGE','SEX'),
                                  plot = "log_pval")
ggsave("factors_corr_covariates.png",plot = p4)

sample_group = samples_metadata(MOFAobject.trained)
rownames(sample_group) = sample_group[,1]
h1 = heatmap_plot = plot_data_heatmap(MOFAobject.trained, 
                                  view = "WGCNA",
                                  factor = 1, features = 20,
                                  cluster_rows = TRUE, cluster_cols = TRUE,
                                  annotation_samples = sample_group[,"responder",drop=F],
                                  show_rownames = TRUE, show_colnames = FALSE,
                                  scale = "row")
ggsave("factor1_WGCNA_responder_heatmap.png",plot = h1, width = 10, height = 10)

factors <- MOFA2::get_factors(MOFAobject.trained, factors = "all")[[1]]
clusters= cluster_samples(MOFAobject.trained,k=2,factors = 1)
samples_metadata$Cluster = factor(clusters$cluster)

df2 <- as.data.frame(factors)
df2$Cluster <- factor(clusters$cluster)
df2$visit = factor(sample_metadata$visit)
df2$responder = factor(sample_metadata$responder)

c1 = ggplot(df2, aes(x = Factor1, y = Factor2,
                    color = responder, shape = Cluster)) +
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
ggsave("factor12_clustersamples_mofa.png",plot = c1, width = 10, height = 10)


set.seed(12345)
MOFAobject_umap <- run_umap(MOFAobject.trained)
MOFAobject_tsne <- run_tsne(MOFAobject.trained)

p <- plot_dimred(
  MOFAobject_umap,
  method = "UMAP",
  color_by = "Factor1",
  dot_size = 3
)

ggsave("mofa_umap_factor1.png", p, width = 6, height = 5)
category.colors <- c(
  "HR" = "darkturquoise", 
  "LR" = "magenta3", "MR" = "#FFCC99")
plot_factors(MOFAobject.trained, 
             factors = c(1,2), 
             color_by = "responder", 
             dot_size = 4) + 
  scale_fill_manual(values=category.colors) +
  stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) +
  scale_color_manual(values=category.colors)


covariates = QIV1_meta_filt_cond_V2 %>% dplyr::select(SubjectIDNew,SEX,AGE)
covariates = column_to_rownames(covariates, var = "SubjectIDNew")
covariates <- covariates %>%
mutate(SEX = ifelse(SEX == "F", 1, 0))

plot_variance_explained_per_feature(MOFAobject.trained,view="WGCNA",split_by_factor = TRUE, features = 10, factors = 2)


library(dplyr)
library(ggplot2)
weights =MOFA2::get_weights(MOFAobject.trained, as.data.frame = TRUE)
weights_factor1 = weights %>% filter(factor == "Factor2")

weights_factor_olink = weights_factor1 %>% filter(view== "olink")
colnames(weights_factor_olink)[1] = "uniprotid"
weights_factor_olink = left_join(weights_factor_olink,olink_ENSG_HGNC, by = "uniprotid")
weights_factor_olink = weights_factor_olink[,-c(1,5)]
colnames(weights_factor_olink)[4] = "feature"
weights_factor_olink= weights_factor_olink[, c("feature", "factor", "value","view")]
weights_factor1 = weights_factor1[-c(45:135),]
weights_factor1 = rbind(weights_factor1, weights_factor_olink)



top_pos <- weights_factor1 %>%
  group_by(view) %>%
  arrange(desc(value), .by_group = TRUE) %>%
  slice(1:10)

top_neg <- weights_factor1 %>%
  group_by(view) %>%
  arrange(value, .by_group = TRUE) %>%
  slice(1:10)

top_features <- bind_rows(top_pos, top_neg)

q3 <- ggplot(top_features,
            aes(x = reorder(feature, value),
                y = value,
                fill = value > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~view, scales = "free_y") +
  scale_fill_manual(values = c("#6A3D9A","#FF7F00"), guide = "none") +
  theme_classic() +
  theme(text = element_text(family = "Arial"), axis.text.y = element_text(size = 16),   
        axis.title.y = element_text(size = 16),  
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 16)) +
  labs(
    x = "",
    y = "MOFA weight",
    title = "Top drivers of Factor 2 by view Day0"
  )

ggsave("Factor2_top_drivers_by_view_v1.png", plot = q3, width = 15, height = 10, dpi = 300)


ssGSEA_V1_BTMscores <- read_csv("~/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/baseline_predictors/ssGSEA_V1_BTMscores.csv")
colnames(ssGSEA_V1_BTMscores)[1] = "Feature"
ssGSEA_V1_BTMscores = pivot_longer(ssGSEA_V1_BTMscores, cols = -Feature,names_to = "sample", values_to = "value" )
colnames(ssGSEA_V1_BTMscores)[3] = "BTMscore"
ssGSEA_V1_BTMscores = ssGSEA_V1_BTMscores %>% filter(Feature %in% BTM_long_sel$Feature)
BTM_avg = ssGSEA_V1_BTMscores %>%
  group_by(Feature) %>%
  summarise(mean_BTMscore = mean(BTMscore))
colnames(BTM_avg)[1]= "feature"
BTM_avg  = left_join(BTM_avg,top_features)
ssGSEA_V1_BTMscores = ssGSEA_V1_BTMscores %>% filter(Feature %in% top_pos4$feature[1:3])

ssGSEA_V1_BTMscores = left_join(ssGSEA_V1_BTMscores,df2_long,by="sample")

scat = ggplot(ssGSEA_V1_BTMscores, aes(x = Value.x, y = BTMscore)) +
geom_point(size = 3) +
theme_classic()+ stat_cor(method = "spearman") +
labs(
x = "MOFA factor4 value",
y = "BTM scores"
)
ggsave("BTMvsfactor4_value_scatter.png", scat)


ssGSEA_V1_BTMscores = left_join(ssGSEA_V1_BTMscores, sample_metadata, by = "sample")
plot2 = ggplot(ssGSEA_V1_BTMscores, aes(x=responder, y= BTMscore, fill= responder)) + geom_boxplot(alpha = 0.3)+
  facet_wrap(~Feature)+
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    size = 4,
    step.increase = 0.12
  ) +theme_bw()+
  labs(title = "BTMscores of responder HK Day0", x= "responder", y="BTMscores")
ggsave("topBTM_vsresponder_HK.png", plot2, width = 10)
