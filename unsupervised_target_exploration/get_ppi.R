# there is no guarrente which network favors our method. 
# For PPI based network I should provide at least 2files: 
# the one similiar to the one used in the Sanger MSB paper, the other one based on a more stringent inclusion criteria
# In addition to that I can explore more.
# The output are still two column csv file.


# 1
library(tidyverse)
setwd("~/cluster_scratch/prior/PPI")

# download.file(url = "https://stringdb-static.org/download/protein.links.full.v1.5/9606.protein.links.full.v11.5.txt.gz",
#               destfile = "9606.protein.links.full.v11.5.txt.gz")
# 
# download.file(url= "https://stringdb-static.org/download/protein.physical.links.full.v11.5/9606.protein.physical.links.full.v11.5.txt.gz",
#               destfile= "9606.protein.physical.links.full.v11.5.txt.gz")
# 
# download.file(url= "https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz",
#               destfile= "9606.protein.info.v11.5.txt.gz")
# download.file(url= "https://ndownloader.figshare.com/files/18820808", "ppi.tar.gz")

PPI_proteininfo <- read_tsv("9606.protein.info.v11.5.txt.gz") %>% 
  rename(protein_id= `#string_protein_id`) %>% 
  rename(hgnc_symbol = preferred_name ) %>% 
  select(protein_id, hgnc_symbol) 
dup_protein <- PPI_proteininfo$protein_id[duplicated(PPI_proteininfo$protein_id)]
PPI_proteininfo <- PPI_proteininfo %>% filter(!protein_id %in% dup_protein) %>% drop_na()
  
  

# remove self loop, interaction <900, and duplicated edges.
# unconnected nodes are not included.
PPI_tibble <- read_delim("9606.protein.links.full.v11.5.txt.gz",delim = " ")
PPI_alllinks_simplified <- PPI_tibble %>% 
  select(protein1, protein2, combined_score) %>% 
  filter(combined_score >900) %>% 
  inner_join(PPI_proteininfo, by= c("protein1"= "protein_id")) %>%
  rename("gene1"= "hgnc_symbol") %>%
  inner_join(PPI_proteininfo, by= c("protein2"= "protein_id")) %>%
  rename("gene2"= "hgnc_symbol") %>% 
  filter(gene1!=gene2) %>%
  drop_na() %>%   
  rowwise() %>%
  mutate(grp = paste(sort(c(gene1, gene2)), collapse = "_")) %>%
  group_by(grp) %>%
  slice(1) %>%
  ungroup() %>%  # 115243 why the number is so much different
  select(gene1, gene2) 
  
write_csv(PPI_alllinks_simplified, "PPI_alllinks_simplified.csv")
rm(PPI_tibble)
############################################################################################
# newer version, only physical links
PPI_tibble <- read_delim("9606.protein.physical.links.full.v11.5.txt.gz",delim = " ")
PPI_physicallinks_simplified <- PPI_tibble %>% 
  select(protein1, protein2, combined_score) %>% 
  filter(combined_score >900) %>% 
  inner_join(PPI_proteininfo, by= c("protein1"= "protein_id")) %>%
  rename("gene1"= "hgnc_symbol") %>%
  inner_join(PPI_proteininfo, by= c("protein2"= "protein_id")) %>%
  rename("gene2"= "hgnc_symbol") %>% 
  filter(gene1!=gene2) %>% 
  drop_na() %>%   
  rowwise() %>%
  mutate(grp = paste(sort(c(gene1, gene2)), collapse = "_")) %>%
  group_by(grp) %>%
  slice(1) %>%
  ungroup() %>%  
  select(gene1, gene2) 

write_csv(PPI_physicallinks_simplified, "PPI_physicallinks_simplified.csv")

######################################################################################################
# replicate analysis of the Sanger paper
PPI_old <- read_delim("9606.protein.links.full.v10.5.txt.gz",delim = " ")
PPI_old_anno <- read_tsv("9606.protein.aliases.v10.5.txt.gz") %>% 
  arrange(string_protein_id) %>% 
  filter(str_detect(source, "BioMart_HUGO")  ) %>% 
  select(-source) %>% 
  rename(protein_id= string_protein_id) %>% 
  rename(hgnc_symbol= alias)
  
dup_protein <- PPI_old_anno$protein_id[duplicated(PPI_old_anno$protein_id)]
PPI_old_anno <- PPI_old_anno %>% filter(!protein_id %in% dup_protein) %>% drop_na()

PPI_old_simplified <- PPI_old %>% 
  select(protein1, protein2, combined_score) %>% 
  filter(combined_score >900) %>% 
  inner_join(PPI_old_anno, by= c("protein1"= "protein_id")) %>%
  rename("gene1"= "hgnc_symbol") %>%
  inner_join(PPI_old_anno, by= c("protein2"= "protein_id")) %>%
  rename("gene2"= "hgnc_symbol") %>% 
  select(gene1, gene2) %>% 
  filter(gene1!=gene2) %>%
  drop_na() %>%   
  rowwise() %>%
  mutate(grp = paste(sort(c(gene1, gene2)), collapse = "_")) %>% 
  group_by(grp) %>%
  slice(1) %>%
  ungroup() %>%  
  select(-grp) 

#205251 Now I got exactly the same number here
write_csv(PPI_old_simplified, "PPI_old_simplified.csv")

# protein_list <-  str_split(unique(c(PPI_tibble$protein1, PPI_tibble$protein2)),
#                            pattern = "\\.",simplify = T,n = 2)[,2]
# library(biomaRt)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=protein_list,mart= mart)
# 
# G_list <- G_list  %>%   mutate(protein_id= str_c("9606",ensembl_peptide_id,sep = "." ))
# write_csv(G_list, "PPI_ENSPID_to_hgnc_symbol.csv")
# G_list <- read_csv("PPI_ENSPID_to_hgnc_symbol.csv")
# PPI_tibble_enhanced <- PPI_tibble %>% 
#   inner_join(G_list, by= c("protein1"= "protein_id")) %>% 
#   rename("gene1"= "hgnc_symbol") %>% 
#   inner_join(G_list, by= c("protein2"= "protein_id")) %>% 
#   rename("gene2"= "hgnc_symbol") 
# 
# write_csv(PPI_tibble_enhanced, "PPI_enhanced.csv" )

#
# check how to simplify PPI with igraph
# https://stackoverflow.com/questions/40086031/remove-self-loops-and-vertex-with-no-edges
# and also use the score threshold to filter the list

# check get_edges function from https://github.com/EmanuelGoncalves/dtrace/blob/a81e379e6ca511cc013aa25b0c7309f8fd0a5f16/dtrace/DataImporter.py#L560 

# https://figshare.com/articles/dataset/Integration_of_pharmacological_and_CRISPR-Cas9_screens_informs_on_drug_mechanism_of_action_in_cancer_cells/10338413/1?file=18820808



# check igraph 
# library(igraph)
# tmp <- PPI_physicallinks_simplified %>% select_at(4:5) 
# tmp1 <- graph_from_data_frame(d = tmp)
# tmp1
# tmp2 <- simplify(tmp1)
# tmp2
# 
# tmp <- PPI_tibble %>% 
#   select(protein1, protein2, combined_score) %>% 
#   filter(combined_score >900) %>% 
#   inner_join(PPI_proteininfo, by= c("protein1"= "protein_id")) %>%
#   rename("gene1"= "hgnc_symbol") %>%
#   inner_join(PPI_proteininfo, by= c("protein2"= "protein_id")) %>%
#   rename("gene2"= "hgnc_symbol") %>% 
#   filter(gene1!=gene2) %>%
#   drop_na() %>%
#   select(4:5)
# 
# tmp1 <- graph_from_data_frame(d = tmp)
# tmp1
# is_simple(tmp1)
# tmp2 <- simplify(tmp1)
# tmp2
# is_simple(tmp2)
# 
# tmp3 <- PPI_alllinks_simplified <- PPI_tibble %>% 
#   select(protein1, protein2, combined_score) %>% 
#   filter(combined_score >900) %>% 
#   inner_join(PPI_proteininfo, by= c("protein1"= "protein_id")) %>%
#   rename("gene1"= "hgnc_symbol") %>%
#   inner_join(PPI_proteininfo, by= c("protein2"= "protein_id")) %>%
#   rename("gene2"= "hgnc_symbol") %>% 
#   filter(gene1!=gene2) %>%
#   drop_na() %>%   #230,524
#   rowwise() %>%
#   mutate(grp = paste(sort(c(gene1, gene2)), collapse = "_")) 
# 
# tmp3 <- tmp3 %>% mutate(duplicates = duplicated(grp))
# 
# sum(duplicated(tmp3$grp))