library(readxl)
kemdata = read_excel("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/Post DBL_Incentive-QIV1-Updated_KEM reviewed_16 Sep 2024.xlsx", sheet = "MEDICAL HISTORY6")

kem_QIV1id = read_excel("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/Post DBL_Incentive-QIV1-Updated_KEM reviewed_16 Sep 2024.xlsx", sheet = "SUBJECT ENROLLMENT11") %>%
  filter(!is.na(`Enrollment Number`)) %>%     
  select(SUBJECTID, `Enrollment Number`)  

subjects_with_conditions <- kemdata %>%
  filter(!is.na(`MH Event Term as per MeDRA`)) %>%
  inner_join(kem_QIV1id, by = "SUBJECTID")%>%
  select(SUBJECTID, `Enrollment Number`, `MH Event Term as per MeDRA`) %>%
  rename(SubjectID = `Enrollment Number`)

kem_diabetes = subjects_with_conditions %>% filter(`MH Event Term as per MeDRA` == "Type 2 diabetes mellitus")
kem_hyper = subjects_with_conditions %>% filter(`MH Event Term as per MeDRA` == "Hypertension")

QIV1_meta_filt = left_join(QIV1_meta_filt,subjects_with_conditions, by = "SubjectID")
colnames(QIV1_meta_filt)[16] = "condition"


#extracting COVID vaccination details

kemdata = read_excel("/home/maziya/INCENTIVE/RNASeq/QIV1_DEG_Analysis/data/Post DBL_Incentive-QIV1-Updated_KEM reviewed_16 Sep 2024.xlsx", sheet = "From KEM HISTORY OF VACCINATION")

covid_vacc = kemdata %>%
  select(`Enrollment ID`,`Vaccine Name`,)%>% 
  rename(SubjectID =`Enrollment ID`, covid_vaccinated = `Vaccine Name`)
QIV1_meta_filt_cond = left_join(QIV1_meta_filt_cond,covid_vacc, by = "SubjectID") 

QIV1_meta_filt_cond = QIV1_meta_filt_cond %>% 
  mutate(
    covid_vaccinated = case_when(
    covid_vaccinated == "Not available" ~ "No_covid_vaccine",
    TRUE ~ covid_vaccinated))

