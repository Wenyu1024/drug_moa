library(tidyverse)
setwd("~/cluster_scratch/prism")
prism_drugs <- read_csv("screened_drugs.csv")
prism_sc_data <- read_csv("secondary-screen-dose-response-curve-parameters.csv")
predictor_df <- read_csv("~/cluster_scratch/ces_21q1_io/ces1_21q1_imputed.csv")

# note drug names should also be unique in prism yet there is two broadcpd id mapping to 
# drug "U-0126", since only after these will we have 1448 drugs as described in the data
# description
# prism_drugs %>% filter(name== "U-0126")
# these two drugs will be merged under the id K18787491, where K91701654 will be discarded
# in addition, there was some problems in the names and smiles 
# for BRD-K05674516, BRD-K35952844 and BRD-K97799481

prism_drugs <- prism_drugs %>% 
  separate(col = broad_id,into=c("BROAD_ID","id",NA),sep="-") %>% 
  unite("BROAD_ID",BROAD_ID:id, sep= "-",remove = T) %>% 
  distinct()

get_unique_smiles <- function(x){
  tmp  <- x %>% str_split(pattern = ",",simplify = F) %>% unlist()
  y= unique(trimws(tmp))
  # if (y == ""){y <-  NA}
  if (length(y)>1){y <-  str_c(y,collapse = ",")}
  # print(y)
  return(y)
}

# for (i in 1:20){get_unique_smiles(prism_drugs$smiles[i]) ;print(i)}

prism_drugs <-  prism_drugs %>% 
  mutate(smiles= map_chr(smiles,get_unique_smiles) )
# prism_drugs_smiles <- prism_drugs %>% 
#   select(smiles)
#     separate_rows(smiles, sep= ",") %>%
#     mutate_all(trimws) %>%
#     distinct()

prism_drugs$BROAD_ID[prism_drugs$BROAD_ID == "BRD-K91701654"] <- "BRD-K18787491"

prism_drugs$name[prism_drugs$BROAD_ID == "BRD-K05674516"] <- "sofosbuvir"
prism_drugs$moa[prism_drugs$BROAD_ID == "BRD-K05674516"] <- "RNA polymerase inhibitor"

prism_drugs$name[prism_drugs$BROAD_ID == "BRD-K35952844"] <- "gluceptate"
# prism_drugs$smiles[prism_drugs$BROAD_ID == "BRD-K35952844"] <- "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)C(O)=O"

prism_drugs$name[prism_drugs$BROAD_ID == "BRD-K97799481"] <- "theophylline"
# prism_drugs$smiles[prism_drugs$BROAD_ID == "BRD-K97799481"] <- "Cn1c2nc[nH]c2c(=O)n(C)c1=O"
prism_drugs$target[prism_drugs$BROAD_ID == "BRD-K97799481"] <- "ADORA1, ADORA2A, ADORA2B, ADORA3, HDAC2, PDE3A, PDE3B, PDE4A, PDE4B, PDE4C, PDE4D"
# If there are multiples smiles only the first one is kept
# prism_drugs_smiles <- prism_drugs %>% 
#   select(smiles) %>% 
#   separate_rows(smiles, sep= ",") %>% 
#   mutate_all(trimws) %>% 
#   distinct()


sensitivity <- prism_sc_data %>% 
  select(depmap_id,auc,ec50, broad_id) %>% 
  separate(col = broad_id,into=c("BROAD_ID","id",NA),sep="-") %>% 
  unite("BROAD_ID",BROAD_ID:id, sep= "-") 
sensitivity$BROAD_ID[sensitivity$BROAD_ID == "BRD-K91701654"] <- "BRD-K18787491"

sensitivity <- sensitivity %>% 
  group_by(BROAD_ID, depmap_id) %>% 
  summarize_all(median,na.rm=T)

# tmp <-prism_sc_data  %>% select(broad_id, name) %>%
#   separate(col = broad_id,into=c("BROAD_ID","id",NA),sep="-") %>%
#   unite("BROAD_ID",BROAD_ID:id, sep= "-") %>%
#   distinct()
# tmp$BROAD_ID[tmp$BROAD_ID == "K91701654"] <- "K18787491"
# tmp[duplicated(tmp$name),]

sc_druglist <- unique(sensitivity$BROAD_ID)
prism_drugs <- prism_drugs %>% 
  mutate(phase= case_when(
    BROAD_ID %in% sc_druglist ~ 2,
    TRUE ~ 1
  ))

# repurposinghub_drug <- read_csv("~/cluster_scratch/prism/drug_repurposing_hub_nobroadid.csv")
# repurposinghub_sample <- read_csv("~/cluster_scratch/prism/drug_repurposing_hub_sample.csv")
# repurposinghub_drug <- repurposinghub_drug %>% 
#   left_join(repurposinghub_sample) %>% distinct()
# rm(repurposinghub_sample)
# prism_drugs <- prism_drugs %>% 
#   left_join(repurposinghub_drug) 

# sensitivity <- prism_sc_data %>% 
#   select(name, broad_id) %>% 
#   distinct() %>% 
#   separate(col = broad_id,into=c("BROAD_ID","id",NA),sep="-") %>% 
#   unite(BROAD_ID, id, sep= "-") %>% 
#   distinct() 
# sensitivity[duplicated(sensitivity$name),]

length(unique(sensitivity$BROAD_ID)) #1448
length(unique(sensitivity$depmap_id)) #481

prism_data <- prism_drugs %>% 
  filter(BROAD_ID %in% sc_druglist) %>% 
  select(-screen_ids, -phase) %>%
  distinct() %>% 
  nest_join(y=sensitivity ,by= "BROAD_ID") %>% 
  mutate(n_forward_modelling= map_int(.x=sensitivity, 
                                      .f=function(x){sum(x %>% pull(depmap_id) %in% predictor_df$DepMap_ID)}))



save(prism_data, file = "prism_data.RData")
