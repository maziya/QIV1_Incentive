library(MOFA2)

transcriptomics <- transcriptomics %>%
  tibble::rownames_to_column("Gene")
transcriptomics_long = transcriptomics %>%
  pivot_longer(cols = -Gene,
               names_to = "SubjectID",
               values_to = "Expression")
metabolomics = as.data.frame(metabolomics)
metabolomics <- metabolomics %>%
  tibble::rownames_to_column("Metabolite")

metabolomics_long = metabolomics %>%
  pivot_longer(cols = -Metabolite,
               names_to = "SubjectID",
               values_to = "Expression")

colnames(Olink_long)[1] = "Feature"

transcriptomics_long<- transcriptomics_long %>%
rename(sample = SubjectID,
value = expression) %>%
mutate(view = "transcriptomics")

proteomics_long<- proteomics_long %>%
  rename(sample = SubjectID,
  value = Expression) %>%
  mutate(view = "proteomics")

proteomics_long <- proteomics_long %>%
  filter(sample %in% transcriptomics_long$sample)

Olink_long<- Olink_long %>%
  rename(sample = SubjectID,
         value = Expression) %>%
  mutate(view = "olink")

Olink_long <- Olink_long %>%
  filter(sample %in% transcriptomics_long$sample)

metabolomics_long<- metabolomics_long %>%
  rename(sample = SubjectID,
         value = Expression) %>%
  mutate(view = "metabolomics")

metabolomics_long <- metabolomics_long %>%
  filter(sample %in% transcriptomics_long$sample)

data = bind_rows(transcriptomics_long,proteomics_long,Olink_long,metabolomics_long)
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
print(model_opts)
train_opts = get_default_training_options(MOFAobject)
MOFAobject = prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


outfile = file.path(getwd(),"QIVIMOFA_model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk=TRUE)
v = plot_variance_explained(MOFAobject.trained)
ggsave("totalvariancebyfactor_mofa.png", plot = v)
print(MOFAobject.trained@cache$variance_explained$r2_per_factor)

var_plot <- plot_variance_explained(MOFAobject.trained, plot_total = T)[[2]]
ggsave("totalvariance_mofa.png", plot = var_plot)

print(MOFAobject.trained@cache$variance_explained$r2_total)


#add metadata
sample_metadata = data.frame(sample=QIV1_meta_filt$SubjectIDNew,
                                           responder =QIV1_meta_filt$Responder_HongKong,visit = QIV1_meta_filt$Visit)
samples_metadata(MOFAobject.trained) = sample_metadata

p1= plot_factor(MOFAobject.trained, 
            factor = 4:5,
            color_by = "responder",
            dot_size = 3,
            dodge = T,
            legend = T,
            add_violin = T,
            violin_alpha = 0.25)
library(ggplot2)
p1 = p1 + scale_color_manual(values = c("MR"="black", "LR"="blue", "HR"="purple"))+
  scale_fill_manual(values = c("MR"="black", "LR"="blue","HR"="purple"))
print(p)

ggsave("factor4&5_olink&transcriptomics_visit.png",plot = p)

#The sign of the weights indicates the direction of the effect: a positive weights 
#indicates that the feature has higher levels in the samples with positive factor 
#values, and vice-versa.

q = plot_weights(MOFAobject.trained,view = "olink", factors = 4, 
                 nfeatures = 10, scale = T, abs = F, text_size = 2.5)
print(q)
ggsave("factor4_olink_weights.png",plot = q)

plot_top_weights(MOFAobject.trained,view = "transcriptomics",
                 factor = 4,
                 nfeatures = 10,    
                 scale = T)

#which genes drive the factor based on weights and get an association of features
#and factors in a scatterplot
q1 = plot_data_scatter(MOFAobject.trained, 
                  view = "transcriptomics", 
                  factor = 1, 
                  features = 12,
                  color_by = "responder",
                  add_lm = T,
                  dot_size = 1)
ggsave("factor1_transcriptomics_responder_weightsscatter1.png",plot = q1, width = 10)

plot_factor_cor(MOFAobject.trained) #check correlation of factors
#they should not be correlated that would mean it is a poor model fit

ggsave("allfactors_Corr.png",plot = p3)

sample_group = samples_metadata(MOFAobject.trained)
rownames(sample_group) = sample_group[,1]
h1 = heatmap_plot = plot_data_heatmap(MOFAobject.trained, 
                                  view = "transcriptomics",
                                  factor = 1, features = 50,
                                  cluster_rows = TRUE, cluster_cols = TRUE,
                                  annotation_samples = sample_group[,"responder",drop=F],
                                  show_rownames = TRUE, show_colnames = FALSE,
                                  scale = "row")
ggsave("factor1_transcriptomics_responder_heatmap.png",plot = h1, width = 10, height = 10)
