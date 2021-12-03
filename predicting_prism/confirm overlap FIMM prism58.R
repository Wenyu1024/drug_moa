# check whether any non-oncology drugs are within the FIMM library.

FIMM_compounds_Annotation_ZIA <- read_csv("~/cluster_scratch/prior/FIMM_compounds_Annotation_ZIA.csv")
prism_nononcology_dfx <- read_csv("~/cluster_wrk/drug_moa/predicting_prism/prism_nononcology_df.csv")  
targetpanel_drh <- read_tsv("~/cluster_scratch/prior/repurposing_drugs_20200324.txt") %>% 
  select(-disease_area, -indication) %>% 
  drop_na()

PRISM_58inchikey <- read_csv("~/cluster_scratch/prism/drug_repurposing_hub_sample.csv") %>% 
  select(broad_id,pert_iname,smiles,InChIKey,pubchem_cid) %>%
  separate(col = broad_id,into=c("BROAD_ID","id",NA),sep="-") %>% 
  unite("BROAD_ID",BROAD_ID:id, sep= "-",remove = T) %>% 
  distinct() %>% 
  filter(BROAD_ID %in% prism_nononcology_df$drug) %>% 
  slice(-48) %>% 
  distinct()
tmp <- inner_join(PRISM_58inchikey, 
                  FIMM_compounds_Annotation_ZIA, 
                  by=c("InChIKey"= "standaryInchiKet")) %>% 
  inner_join(prism_nononcology_df, by= c("BROAD_ID" = "drug"))
write_csv(tmp, "~/cluster_wrk/drug_moa/predicting_prism/FIMM_prism58_overlap.csv")
  