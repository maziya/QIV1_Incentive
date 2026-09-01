input_dir <- "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/strainwise_DEG_salmon/HRLR"    
output_dir <- "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/abslogFC_strainwise_signature" 
baseline_signatures_various = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/baseline_signatures_various.csv")

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

purrr::walk(csv_files, function(file_path) {
  
  base_name <- tools::file_path_sans_ext(basename(file_path))
  deg_df <- read_csv(file_path, show_col_types = FALSE)
  
  clean_deg <- deg_df %>% 
    filter(is.finite(logFC)) %>%
    mutate(abs_logFC = abs(logFC))
  
  # 3. Iterate over every signature directly
  wilcox_results <- colnames(baseline_signatures_various) %>%
    purrr::map_dfr(function(sig_name) {
      
      sig_genes <- na.omit(baseline_signatures_various[[sig_name]])
      
      test_dat <- clean_deg %>%
        mutate(group = if_else(HGNC_symbol %in% sig_genes, "Signature", "Background"))
      
      x <- test_dat %>% filter(group == "Background") %>% pull(abs_logFC)
      y <- test_dat %>% filter(group == "Signature") %>% pull(abs_logFC)
      
      n1 <- length(x)
      n2 <- length(y)
      
      if (n1 < 2 || n2 < 2) {
        return(tibble(signature = sig_name, statistic = NA, p = NA,
                      effsize = NA, magnitude = NA,
                      estimate = NA, conf.low = NA, conf.high = NA,
                      n_bg = n1, n_sig = n2))
      }
      
      test_res <- test_dat %>%
        rstatix::wilcox_test(abs_logFC ~ group, ref.group = "Background", exact = FALSE)
      
      eff_res <- test_dat %>%
        rstatix::wilcox_effsize(abs_logFC ~ group, ref.group = "Background") %>%
        dplyr::select(effsize, magnitude)
      
      ci_res <- wilcox.test(y, x, conf.int = TRUE, exact = FALSE)
      
      test_res %>%
        dplyr::select(statistic, p) %>% 
        bind_cols(eff_res) %>%
        mutate(
          signature = sig_name, 
          estimate = ci_res$estimate,
          conf.low = ci_res$conf.int[1],
          conf.high = ci_res$conf.int[2],
          n_bg = n1,
          n_sig = n2
        ) %>%
        relocate(signature)
    }) %>%
    # 4. Apply BH correction across all signatures for THIS specific CSV
    mutate(p.adj = p.adjust(p, method = "BH")) %>%
    add_significance("p.adj") %>%
    mutate(source_strain = base_name, .before = 1)
  
  output_file <- file.path(output_dir, paste0(base_name, "_wilcox_results.csv"))
  write_csv(wilcox_results, output_file)
  
})


#########################################################################################

input_dir <- "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/strainwise_DEG_salmon/HRLR"    
output_dir <- "/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/abslogFC_strainwise_alldegs_vs_1sig_atatime"
baseline_signatures_various = read.csv("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/results/results_Aug2026/baseline_signatures_various.csv")


csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# combine ALL CSV files into one master dataframe
# The set_names function uses the file name to create a 'source_strain' identifier column
combined_deg_data <- csv_files %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE) %>%
            filter(is.finite(logFC)) %>%
            mutate(abs_logFC = abs(logFC)), 
          .id = "source_strain")

# Iterate over every signature
master_results <- colnames(baseline_signatures_various) %>%
  purrr::map_dfr(function(sig_name) {
    
    # Extract genes for this specific signature
    sig_genes <- na.omit(baseline_signatures_various[[sig_name]])
    
    # Run the test for every strain for this ONE signature
    sig_results <- combined_deg_data %>%
      mutate(group = if_else(HGNC_symbol %in% sig_genes, "Signature", "Background")) %>%
      split(.$source_strain) %>%
      purrr::map_dfr(function(strain_dat) {
        
        x <- strain_dat %>% filter(group == "Background") %>% pull(abs_logFC)
        y <- strain_dat %>% filter(group == "Signature") %>% pull(abs_logFC)
        
        n1 <- length(x)
        n2 <- length(y)
        
        if (n1 < 2 || n2 < 2) {
          return(tibble(statistic = NA, p = NA, effsize = NA, magnitude = NA,
                        estimate = NA, conf.low = NA, conf.high = NA,
                        n_bg = n1, n_sig = n2))
        }
        
        test_res <- strain_dat %>%
          rstatix::wilcox_test(abs_logFC ~ group, ref.group = "Background", exact = FALSE)
        
        eff_res <- strain_dat %>%
          rstatix::wilcox_effsize(abs_logFC ~ group, ref.group = "Background") %>%
          dplyr::select(effsize, magnitude)
        
        ci_res <- wilcox.test(y, x, conf.int = TRUE, exact = FALSE)
        
        test_res %>%
          dplyr::select(statistic, p) %>%
          bind_cols(eff_res) %>%
          mutate(
            estimate = ci_res$estimate,
            conf.low = ci_res$conf.int[1],
            conf.high = ci_res$conf.int[2],
            n_bg = n1,
            n_sig = n2
          )
      }, .id = "source_strain")
    
    # 4. Apply BH correction across all strains for THIS ONE signature
    sig_results %>%
      mutate(p.adj = p.adjust(p, method = "BH")) %>%
      add_significance("p.adj") %>%
      mutate(signature = sig_name, .before = 1)
    
  })

output_file <- file.path(output_dir, "abslogFC_wilcoxresults_by_signature.csv")
write_csv(master_results, output_file)

