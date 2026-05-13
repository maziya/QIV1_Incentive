# script to find proteins whose change in concentration 
# is related to HAI foldchange for each strain using a linear model

meta_V1 = QIV1_responder %>% filter(Visit == "V1")
meta_V2 = QIV1_responder %>% filter(Visit == "V2") 

meta_V1 = meta_V1[order(meta_V1$SubjectID), ]
meta_V2 = meta_V2[order(meta_V2$SubjectID), ]

prot_V1 = untargeted_matrix_untransformed[, meta_V1$IDVisit]
prot_V2 = untargeted_matrix_untransformed[, meta_V2$IDVisit]

prot_V1[] <- lapply(prot_V1, function(x) as.numeric(as.character(x)))
prot_V2[] <- lapply(prot_V2, function(x) as.numeric(as.character(x)))

prot_V1 <- as.matrix(prot_V1)
prot_V2 <- as.matrix(prot_V2)

colnames(prot_V1) <- gsub("V1$", "", colnames(prot_V1))
colnames(prot_V2) <- gsub("V2$", "", colnames(prot_V2))

delta_prot = prot_V2 - prot_V1
delta_prot = scale(delta_prot)

strains = c("FC.A.HongKong", "FC.A.Victoria", "FC.B.Phuket","FC.B.Washington")
all_results <- lapply(strains, function(strain) {
  
  lm_results <- lapply(1:nrow(delta_prot), function(i) {
    x <- delta_prot[i, ]
    df = data.frame(
      response = meta_V1[[strain]],
      delta = x,
      sex = meta_V1$Sex,
      batch = meta_V1$Batch,
      age = meta_V1$AGE
    )
    fit =lm(response ~ delta + sex + batch + age, data = df)
    coef <- summary(fit)$coefficients
    as.numeric(coef["delta", ])
  })
  
  lm_mat <- do.call(rbind, lm_results)
  lm_df <- as.data.frame(lm_mat)
  
  lm_df$protein <- rownames(delta_prot)
  lm_df$strain <- strain
  
  colnames(lm_df)[1:4] <- c("beta", "se", "t", "pval")
  lm_df$padj <- p.adjust(lm_df$pval, method = "BH")
  
  lm_df
})

final_df = do.call(rbind, all_results) 

library(ggplot2)
for (strain_name in unique(final_df$strain)) {
  
  lm_df <- final_df %>% filter(strain == strain_name)
  
  lm_df$log10p <- -log10(lm_df$pval)
  lm_df$significant <- lm_df$pval < 0.05
  
  top_hits <- lm_df %>% arrange(pval) %>% head(10)
  top_hits <- inner_join(top_hits, untargeteduniprot_ENSG_HGNC, by = "protein")
  
  plot <- ggplot(lm_df, aes(x = beta, y = log10p, color = significant)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c("grey", "blue")) +
    geom_text_repel(data = top_hits,
                    aes(label = HGNC_symbol),
                    size = 3) +
    theme_bw() +
    ggtitle(strain_name)+ xlab("effect size beta")
  
  ggsave(paste0("volcano_", strain_name, ".png"), plot = plot, height = 8, width = 8)
}

all_results_corr <- lapply(strains, function(strain) {
  
  cor_results <- lapply(1:nrow(delta_prot), function(i) {
    
    x = delta_prot[i, ]
    y = meta_V1[[strain]]
    
    test = cor.test(x, y, method = "spearman")  
    
    c(cor = test$estimate,
      pval = test$p.value)
  })
  
  cor_mat <- do.call(rbind, cor_results)
  cor_df <- as.data.frame(cor_mat)
  
  cor_df$protein <- rownames(delta_prot)
  cor_df$strain <- strain
  
  cor_df$padj <- p.adjust(cor_df$pval, method = "BH")
  
  cor_df
})

final_cor_df <- do.call(rbind, all_results_corr)

final_cor_df_fil = final_cor_df %>% filter(padj < 0.05)


##for metaboolomics 
valid_ids <- meta_V1$IDVisit[meta_V1$IDVisit %in% colnames(metaboliteconc_batchadjusted)]
valid_ids <- meta_V2$IDVisit[meta_V2$IDVisit %in% colnames(metaboliteconc_batchadjusted)]
metabolomics_V1 <- metaboliteconc_batchadjusted[, valid_ids]
metabolomics_V2 <- metaboliteconc_batchadjusted[, valid_ids]

meta_V1 = meta_V1 %>% filter(IDVisit %in% valid_ids)
meta_V2 = meta_V2 %>% filter(IDVisit %in% valid_ids)
  
metabolomics_V1[] <- lapply(metabolomics_V1, function(x) as.numeric(as.character(x)))
metabolomics_V2[] <- lapply(metabolomics_V2, function(x) as.numeric(as.character(x)))

metabolomics_V1 <- as.matrix(metabolomics_V1)
metabolomics_V2 <- as.matrix(metabolomics_V2)

colnames(metabolomics_V1) <- gsub("V1$", "", colnames(metabolomics_V1))
colnames(metabolomics_V2) <- gsub("V2$", "", colnames(metabolomics_V2))

delta_metabolomics = metabolomics_V2 - metabolomics_V1
delta_metabolomics = scale(delta_metabolomics)

strains = c("FC.A.HongKong", "FC.A.Victoria", "FC.B.Phuket","FC.B.Washington")
all_results <- lapply(strains, function(strain) {
  
  lm_results <- lapply(1:nrow(delta_metabolomics), function(i) {
  x <- delta_metabolomics[i, ]
  df = data.frame(
    response = meta_V1[[strain]],
    delta = x,
    sex = meta_V1$Sex,
    batch = meta_V1$Batch,
    age = meta_V1$AGE
  )
  fit =lm(response ~ delta + sex + batch + age, data = df)
  coef <- summary(fit)$coefficients
  as.numeric(coef["delta", ])
  })
  
  lm_mat <- do.call(rbind, lm_results)
  lm_df <- as.data.frame(lm_mat)
  
  lm_df$metabolite <- rownames(delta_metabolomics)
  lm_df$strain <- strain
  
  colnames(lm_df)[1:4] <- c("beta", "se", "t", "pval")
  lm_df$padj <- p.adjust(lm_df$pval, method = "BH")
  
  lm_df
})

final_df = do.call(rbind, all_results) 

library(ggplot2)
for (strain_name in unique(final_df$strain)) {
  
  lm_df <- final_df %>% filter(strain == strain_name)
  
  lm_df$log10p <- -log10(lm_df$pval)
  lm_df$significant <- lm_df$pval < 0.05
  
  top_hits <- lm_df %>% arrange(pval) %>% head(10)

  plot <- ggplot(lm_df, aes(x = beta, y = log10p, color = significant)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c("grey", "blue")) +
    geom_text_repel(data = top_hits,
                    aes(label = metabolite),
                    size = 3) +
    theme_bw() +
    ggtitle(strain_name)+ xlab("effect size beta")
  
  ggsave(paste0("volcano_metabolites_", strain_name, ".png"), plot = plot, height = 8, width = 8)
}



all_results_corr <- lapply(strains, function(strain) {
  
  cor_results <- lapply(1:nrow(delta_metabolomics), function(i) {
    
    x = delta_metabolomics[i, ]
    y = meta_V1[[strain]]
    
    test = cor.test(x, y, method = "spearman")  
    
    c(cor = test$estimate,
      pval = test$p.value)
  })
  
  cor_mat <- do.call(rbind, cor_results)
  cor_df <- as.data.frame(cor_mat)
  
  cor_df$metabolite <- rownames(delta_metabolomics)
  cor_df$strain <- strain

  cor_df$padj <- p.adjust(cor_df$pval, method = "BH")
  
  cor_df
})

final_cor_df <- do.call(rbind, all_results_corr)

final_cor_df_fil = final_cor_df %>% filter(pval < 0.05)


###########
# for each strain get the list of lm genes and beta
for (strain_name in strains) {
  
  df = prot_lm_FC %>% 
    filter(strain == strain_name)
  
  write.csv(df, paste0(strain_name, "_prot_lm_FC.csv"), row.names = FALSE, quote = FALSE)
}


metlist = c("Tyrosine","Asparagine","Lysine","Valine","Isoleucine","Leucine")
metlist_V1= metabolomics_V1[metlist,]
metlist_V2= metabolomics_V2[metlist,]

metlist_V1 = as.data.frame(metlist_V1)
metlist_V2 = as.data.frame(metlist_V2)

metlist_V1$Metabolite <- rownames(metlist_V1)
metlist_V2$Metabolite <- rownames(metlist_V2)

metlist_V1_long <- metlist_V1 %>%
  pivot_longer(
    cols = -Metabolite,
    names_to = "SubjectID",
    values_to = "Value"
  ) %>%
  mutate(Visit = "Day0")

metlist_V2_long <- metlist_V2%>%
  pivot_longer(
    cols = -Metabolite,
    names_to = "SubjectID",
    values_to = "Value"
  ) %>%
  mutate(Visit = "Day3")


plot_df <- bind_rows(metlist_V1_long, metlist_V2_long)
plot_df$Visit <- factor(plot_df$Visit,
                        levels = c("Day0", "Day3"))

p <- ggplot(plot_df,
            aes(x = Visit,
                y = Value)) +
  
  geom_line(aes(group = SubjectID),
            color = "grey60",
            alpha = 0.4) +
  
  geom_point(aes(color = Visit),
             size = 2,
             alpha = 0.8) +
  
  geom_boxplot(aes(fill = Visit),
               alpha = 0.4,
               width = 0.5,
               outlier.shape = NA) +
  facet_wrap(~Metabolite,
             scales = "free_y") +
  
  theme_classic(base_size = 14) +
  
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45,
                               hjust = 1)
  ) +
  
  labs(
    x = "",
    y = "Metabolite Concentrations"
  )

ggsave("mets_associated_response.png", plot = p, height= 10, width = 10, dpi = 600)
