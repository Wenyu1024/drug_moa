# prepare first the whole dataset but with ids that can later be matched to ctd2.
# I will run first on gdsc independently. 
# Do I have to run all the analysis on GDSC? perhaps first only do the sensitivity fitting.
# and then if the reviewer asked, do the rest.

# to run the analysis I need to map cell id. basically map Sanger_Model_ID to depmap_iD
# using the sample_info_broad.csv file

setwd("~/cluster_scratch/gdsc")
# gdsc1 <- readxl::read_excel("GDSC1_fitted_dose_response_25Feb20.xlsx")
gdsc2 <- readxl::read_excel("GDSC2_fitted_dose_response_25Feb20.xlsx")
sample_info <- read_csv("sample_info_broad.csv")
drug_info <- read_csv("Drug_list.csv")
id_mapping_gdsc <- read_tsv("~/cluster_scratch/gdsc/pubchem_cid_to_smiles.txt",
                            col_names = c("PubCHEM", "smiles"),col_types = "cc")

GDSC_binarytarget <- read_csv("~/cluster_scratch/prior/GDSC_binarytarget_ids.csv") %>% select_at(c(1,4,5))
colnames(GDSC_binarytarget) <- c("DRUG_ID", "gene_target", "DRUG_TYPE")

predictor_df <- read_csv("~/cluster_scratch/ces_21q1_io/ces1_21q1_imputed.csv")


#correct a target naming to make it easier for subsequent splitting
gdsc_sen <- gdsc2 %>% 
  select("CELL_LINE_NAME" , "SANGER_MODEL_ID","DRUG_ID" , "DRUG_NAME","AUC","LN_IC50") %>% 
  group_by(SANGER_MODEL_ID , DRUG_ID) %>% 
  summarise(AUC= mean(AUC), LN_IC50= mean(LN_IC50)) %>% 
  ungroup() %>% 
  inner_join(
    y = sample_info %>% select(DepMap_ID, Sanger_Model_ID) %>% distinct() %>% drop_na(), 
    by= c("SANGER_MODEL_ID" = "Sanger_Model_ID")) %>%
  select(-SANGER_MODEL_ID) %>% 
  # group_by(DRUG_ID) %>% 
  nest(sensitivity=-DRUG_ID)
  # ungroup()

gdsc_df <- 
  gdsc2 %>% 
  select("DRUG_ID" , "DRUG_NAME", "PATHWAY_NAME") %>%
  distinct() %>%
  left_join(GDSC_binarytarget , by= "DRUG_ID") %>% 
  separate_rows(gene_target) %>%
  # group_by(DRUG_ID , DRUG_NAME, PATHWAY_NAME, DRUG_TYPE) %>%
  nest(target= contains("gene_target")) %>%
  # ungroup() %>%
  inner_join(gdsc_sen,by= "DRUG_ID") 

# adding ids to gdsc drugs
# in case there are multiple PubCHEM ids, only the first one is kept
get_unique_id <- function(x, sep =","){
  tmp  <- x %>% str_split(pattern = sep,simplify = F) %>% unlist()
  y= unique(unlist(trimws(tmp)))
  # if (y == ""){y <-  NA}
  # if (length(y)>1){y <-  str_c(y,collapse = ",")}
  # print(y)
  y= y[1]
  return(y)
}


gdsc_df <- gdsc_df %>% 
  left_join(y = drug_info %>%  select(drug_id, PubCHEM) %>% distinct(),
            by = c("DRUG_ID"="drug_id")) %>% 
  left_join(id_mapping_gdsc)  %>% 
  mutate(n_forward_modelling= map_int(.x=sensitivity, 
                                      .f=function(x){sum(x %>% pull(DepMap_ID) %in% predictor_df$DepMap_ID)})) %>% 
  mutate(PubCHEM= map_chr(PubCHEM,get_unique_id) )


# Note that drugs with same structure 
# (same smiles may led to overestimate of chemical fingerprint based method)

# no need to adjust them mannually, just accept them for what it is
# gdsc_df$PubCHEM[gdsc_df$DRUG_NAME== "Buparlisib"] <- "16654980"
# gdsc_df$PubCHEM[gdsc_df$DRUG_NAME== "Ulixertinib"] <- "11719003"
# gdsc_df$PubCHEM[gdsc_df$DRUG_NAME== "Ulixertinib"] <- "11719003"


save(gdsc_df, file = "gdsc2.RData")  
