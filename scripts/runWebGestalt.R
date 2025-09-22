run_WebGestalt = function(csv_file_path,
                           output_csv_path,
                           gene_column,  #Column name containing gene symbols
                           rank_column, #Column name for ranking (p.value, logFC, etc.)
                           input_type = "ranked_list", 
                           enrichment_method = "GSEA",
                           fdr_cutoff = 0.05,
                           organism = "hsapiens",
                           project_dir,
                           project_name = "WebGestalt_Analysis") {
  
  
  input_df = read.csv(csv_file_path, stringsAsFactors = FALSE)
  ranked_df = data.frame(gene = input_df[[gene_column]],
      score = input_df[[rank_column]],
      stringsAsFactors = FALSE)
  interestGene = ranked_df
  
  #Run WebGestaltR
  webgestalt_result <- WebGestaltR::WebGestaltR(
    enrichMethod = enrichment_method,
    organism = organism,
    enrichDatabase = "pathway_Reactome",
    interestGene = interestGene,
    interestGeneType = "genesymbol",
    fdrThr = fdr_cutoff,
    topThr = 100,
    isOutput = TRUE,
    outputDirectory = project_dir,
    projectName = project_name
  )
  
  if (!is.null(webgestalt_result) && nrow(webgestalt_result) > 0) {
    final_results <- webgestalt_result %>%
      dplyr::filter(FDR < fdr_cutoff) %>%
      dplyr::arrange(FDR) %>%
      dplyr::mutate(
        pathway = geneSet,
        padj = FDR,
        leadingEdge = sapply(userId, function(x) ifelse(!is.na(x) & x != "", paste(strsplit(x, ";")[[1]], collapse = ";"), "")))
    
    write.csv(final_results, file = output_csv_path, row.names = FALSE)
    cat(paste("Found", nrow(final_results), "significant pathways\n"))
    return(final_results)
  } else {
    cat("No significant pathways found\n")
    return(data.frame())
  }
}

#folder containing csv files input and output
input_folder = "/home/maziya/INCENTIVE/Proteomics/Untargeted_DEP/Untargeted_strainwise_DEPs/ResponderV2vsV1/webgestalt_reactome"
output_folder = "/home/maziya/INCENTIVE/Proteomics/Untargeted_DEP/Untargeted_strainwise_DEPs/ResponderV2vsV1/webgestalt_reactome/results"

csv_files = list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
results_list = lapply(csv_files, function(file_path) {
  output_csv = file.path(output_folder, paste0(tools::file_path_sans_ext(basename(file_path)), "_WebGestalt.csv"))
  
  run_WebGestalt(
    csv_file_path = file_path,
    output_csv_path = output_csv,
    gene_column = "HGNC_symbol",   
    rank_column = "logFC",        
    enrichment_method = "GSEA",   
    fdr_cutoff = 0.05,
    organism = "hsapiens",
    project_dir = output_folder,      
    project_name = tools::file_path_sans_ext(basename(file_path))
  )
})
