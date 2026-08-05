#for every strain DEG list of HR vs LR at V1 load the DEG results csv file 
#and the baseline signatures gene list
HongKong_HRV1vsLRV1_DEG_results_allgenes = HongKong_HRV1vsLRV1_DEG_results %>%
  dplyr::select(logFC, HGNC_symbol)

# build_dfs <- function(allgenes_df, deg_df, signatures_df) {
#   sig_names <- colnames(signatures_df)
#   sig_list <- lapply(sig_names, function(sig) {
#     deg_df %>%
#       filter(HGNC_symbol %in% na.omit(signatures_df[[sig]])) %>%
#       dplyr::select(HGNC_symbol, logFC) %>%
#       rename(!!paste0("logFC_", sig) := logFC)
#   })
#   names(sig_list) <- sig_names
#   dfs <- c(list(allgenes = allgenes_df), sig_list)
#   return(dfs)
# }
build_dfs <- function(allgenes_df, deg_df, signatures_df) {
  sig_names <- colnames(signatures_df)
  
  sig_list <- lapply(sig_names, function(sig) {
    sig_genes <- na.omit(signatures_df[[sig]])
    
    # signature genes
    sig_df <- deg_df %>%
      filter(HGNC_symbol %in% sig_genes) %>%
      dplyr::select(HGNC_symbol, logFC) %>%
      rename(!!paste0("logFC_", sig) := logFC)
    
    # background = allgenes EXCLUDING this signature's genes
    bg_df <- allgenes_df %>%
      filter(!HGNC_symbol %in% sig_genes) %>%
      rename(!!paste0("logFC_bg_", sig) := logFC)
    
    list(sig = sig_df, bg = bg_df)
  })
  
  names(sig_list) <- sig_names
  
  sig_dfs <- purrr::map(sig_list, "sig")
  bg_dfs  <- purrr::map(sig_list, "bg") %>% setNames(paste0(sig_names, "_bg"))
  
  dfs <- c(sig_dfs, bg_dfs)
  return(dfs)
}


dfs = build_dfs(
  HongKong_HRV1vsLRV1_DEG_results_allgenes ,
  HongKong_HRV1vsLRV1_DEG_results,
  baseline_signatures_various
)
library(purrr)
library(tidyr)
library(rstatix)
library(ggpubr)

combined <- map2_dfr(dfs, names(dfs),
                     ~ mutate(.x, dataset = .y)) %>%
  pivot_longer(cols = starts_with("logFC"),
               names_to = "contrast",
               values_to = "logFC") %>%
  filter(is.finite(logFC))%>%
  mutate(abs_logFC = abs(logFC))



sig_names_only <- combined %>%
  filter(!grepl("_bg$", dataset)) %>%
  distinct(dataset) %>%
  pull(dataset)

wilcox_results <- sig_names_only %>%
  map_dfr(function(d) {
    bg_name <- paste0(d, "_bg")
    
    x <- combined %>% filter(dataset == bg_name) %>% pull(abs_logFC)
    y <- combined %>% filter(dataset == d) %>% pull(abs_logFC)
    
    if (length(x) < 2 || length(y) < 2) {
      return(tibble(group2 = d, estimate = NA, conf.low = NA, conf.high = NA,
                    p.value = NA, n1 = length(x), n2 = length(y)))
    }
    
    broom::tidy(wilcox.test(y, x, conf.int = TRUE, exact = FALSE)) %>%
      mutate(group2 = d, n1 = length(x), n2 = length(y))
  }) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  add_significance("p.adj")



#stepwise brackets for pvalue and * notation
pvals$y.position <- seq(
  max(combined$abs_logFC, na.rm = TRUE) * 1.02,
  by = 1,
  length.out = nrow(pvals)
)

plot <- ggplot(combined, aes(dataset, abs_logFC)) +
  geom_boxplot(outlier.shape = NA) +
  stat_pvalue_manual(
    pvals,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  ) +
  theme_pubr(base_size = 14) +
  labs(x = " Baseline Signature", y = "Gene abs logFC", title = "HRvsLR at Day0 logFC")
ggsave("HRvsLR_geneabslogFC_sig.png", plot = plot, height = 8, width = 15)


# compute p-values vs allgenes using signed logFC
pvals_1 <- combined %>%
  wilcox_test(logFC ~ dataset, ref.group = "allgenes") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() 
%>%
  filter(p.adj < 0.05)

#stepwise brackets for pvalue and * notation
start_y <- max(combined$logFC, na.rm = TRUE) * 1.05
pvals_1$y.position <- start_y + (seq_len(nrow(pvals_1)) - 1) * 0.3

plot1 <- ggplot(combined_sel1, aes(dataset, logFC)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_jitter(width = 0.2, alpha = 0.5) +
  stat_pvalue_manual(
    pvals_1,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  theme_pubr(base_size = 16) +
  labs(x = "Signature", y = "Gene logFC",title = "Washington HRvsLR at Day0 logFC")

ggsave("Washington_gene_logFC_sig2.png", plot = plot1, height = 8, width = 15)

