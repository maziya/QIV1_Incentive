#map uniprot IDs to gene symbols 
library(biomaRt)
library(UniProt.ws)
library(readxl)


olink_ensg = read.table("/home/maziya/INCENTIVE/Proteomics/olink_idmapping_2025_07_24.tsv",sep='\t', header = TRUE)
colnames(olink_ensg)[1] = "uniprotid"
colnames(olink_ensg)[2] = "target_id"
olink_ensg$target_id = gsub("\\..*$","",olink_ensg$target_id)
protein_coding = read.csv("/home/maziya/INCENTIVE/Proteomics/protein_coding_ensemble_hgnclist.csv")

library(dplyr)
olink_ensg = olink_ensg %>% dplyr::inner_join(protein_coding, by = "target_id")
write.csv(olink_ensg,"olink_ENSG_HGNC.csv",row.names = FALSE,quote = FALSE)


#add HGNC gene symbols for the corresponding protein ids before doing GSEA
folderpath= "/home/maziya/INCENTIVE/Proteomics/OLINK/Strainwise_DEPs_OLINK_adjFC"
csv_files = list.files(folderpath, pattern = "\\.csv$", full.names = TRUE)
process_deg_file <- function(file_path, olink_ensg) {
  deg <- read.csv(file_path)
  colnames(deg)[7] <- "uniprotid"
  deg <- deg %>%
    dplyr::inner_join(olink_ensg, by = "uniprotid")
  write.csv(deg, file_path, quote = FALSE, row.names = FALSE)
  return(deg)
}

deg_list <- lapply(csv_files, process_deg_file, olink_ensg = olink_ENSG_HGNC)


#mapping untargeted uniprot ids to ensembl ids
untarg_ensg = read.table("/home/maziya/INCENTIVE/Proteomics/untargeted_idmapping_2025_08_22.tsv",sep='\t', header = TRUE)
colnames(untarg_ensg)[1] = "uniprotid"
colnames(untarg_ensg)[2] = "target_id"
untarg_ensg$target_id = gsub("\\..*$","",untarg_ensg$target_id)
untarg_ensg_uniq = untarg_ensg %>% dplyr::filter(!duplicated(uniprotid))
untarg_ensg_uniq = untarg_ensg_uniq %>% dplyr::filter(!duplicated(target_id))
protein_coding = read.csv("/home/maziya/INCENTIVE/Proteomics/protein_coding_ensemble_hgnclist.csv")

library(dplyr)
untarg_ensg_uniq = untarg_ensg_uniq %>% dplyr::inner_join(protein_coding, by = "target_id")
write.csv(untarg_ensg_uniq,"unatargeteduniprot_ENSG_HGNC.csv",row.names = FALSE,quote = FALSE)

#add HGNC gene symbols for the corresponding protein ids before doing GSEA
folderpath= "/home/maziya/INCENTIVE/Proteomics/Untargeted_DEP/Untargeted_strainwise_adjFClabel_DEPs"
csv_files = list.files(folderpath, pattern = "\\.csv$", full.names = TRUE)
process_deg_file <- function(file_path, untargeteduniprot_ENSG_HGNC) {
  deg <- read.csv(file_path)
  colnames(deg)[7] <- "uniprotid"
  deg$uniprotid = sub("\\|.*", "", deg$uniprotid)
  deg <- deg %>%
    dplyr::inner_join(untargeteduniprot_ENSG_HGNC, by = "uniprotid")
  write.csv(deg, file_path, quote = FALSE, row.names = FALSE)
  return(deg)
}

deg_list <- lapply(csv_files, process_deg_file, untargeteduniprot_ENSG_HGNC)

