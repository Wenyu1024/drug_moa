# prepare first the whole dataset but with ids that can later be matched to ctd2.
# I will run first on gdsc independently. 
# Do I have to run all the analysis on GDSC? perhaps first only do the sensitivity fitting.
# and then if the reviewer asked, do the rest.

# to run the analysis I need to map cell id. basically map Sanger_Model_ID to depmap_iD
# using the sample_info_broad.csv file

setwd("~/cluster_scratch/gdsc")
gdsc1 <- readxl::read_excel("GDSC1_fitted_dose_response_25Feb20.xlsx")
gdsc2 <- readxl::read_excel("GDSC2_fitted_dose_response_25Feb20.xlsx")
sample_info <- read_csv("sample_info_broad.csv")


gdsc_sen <- bind_rows(gdsc1, gdsc2) %>% 
  select("CELL_LINE_NAME" , "SANGER_MODEL_ID","DRUG_ID" , "DRUG_NAME","PUTATIVE_TARGET","AUC","LN_IC50") %>% 
  group_by(SANGER_MODEL_ID , DRUG_ID) %>% 
  summarise(AUC= mean(AUC), LN_IC50= mean(LN_IC50)) %>% 
  ungroup() %>% 
  inner_join(
    y = sample_info %>% select(DepMap_ID, Sanger_Model_ID) %>% distinct() %>% drop_na(), 
    by= c("SANGER_MODEL_ID" = "Sanger_Model_ID")) %>%
  select(-SANGER_MODEL_ID) %>% 
  group_by(DRUG_ID) %>% 
  nest() %>% 
  ungroup()

gdsc_df <- 
  bind_rows(gdsc1, gdsc2) %>%
  select("DRUG_ID" , "DRUG_NAME","PUTATIVE_TARGET", "PATHWAY_NAME") %>%
  distinct() %>%
  separate_rows(PUTATIVE_TARGET) %>%
  group_by(DRUG_ID , DRUG_NAME, PATHWAY_NAME) %>%
  nest() %>%
  ungroup() %>%
  inner_join(gdsc_sen,by= "DRUG_ID") %>% 
  rename( target=data.x, sensitivity=data.y)

save(gdsc_df, file = "gdsc.RData")  
