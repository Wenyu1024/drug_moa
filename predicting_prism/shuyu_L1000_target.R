library(tidyverse)
repurposinghub_drug <- read_csv("~/cluster_scratch/prism/drug_repurposing_hub_nobroadid.csv")
repurposinghub_sample <- read_csv("~/cluster_scratch/prism/drug_repurposing_hub_sample.csv")

repurposinghub_drug1 <- repurposinghub_drug %>% 
  left_join(repurposinghub_sample) %>% distinct() %>% 
  separate(col = broad_id, 
           into=c("BRD","id", NA), 
           sep = "-", remove = T, extra= "drop" ) %>% 
  unite("broad_id", BRD:id, sep = "-") %>% 
  separate(col = deprecated_broad_id, 
           into=c("BRD","id", NA), 
           sep = "-", remove = T, extra= "drop" ) %>% 
  unite("deprecated_broad_id", BRD:id, sep = "-") %>% 
  select(-vendor_name) %>% 
  distinct()

write_csv(repurposinghub_drug1, "~/cluster_scratch/prism/drug_repurposing_hub_alldrugsanno.csv")
