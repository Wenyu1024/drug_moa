library(tidyverse)
setwd("/scratch/project_2003466/depmap_21q1")
list.files()
cell_info <- read_csv("sampleinfo")
ceres<- read_csv("ceres")
cn <- read_csv("cn")
exp_seq <- read_csv("exp_seq")
mut <- read_csv("mut")
demeter2 <- read_csv("demeter2")

# array-based mRNA expression data is quite old. 
exp_array <- data.table::fread("/scratch/project_2003466/depmap_static/exp_array")


# setwd("/home/cloud-user/cluster_scratch/depmap_static//")
# ceres_sanger <- read_csv("ceres_sanger")
# demeter2 <- read_csv("demeter2")
# exp_array <- data.table::fread("exp_array") # read_table2 has its own problem of failling to identify empty colnames
#* "NCIH292_LUNG" at locations 964 and 992
#* 


#ceres:
# colnames(ceres_broad)[2:18120] <- str_split(colnames(ceres_broad)[2:18120],
#                                             n=2,pattern= " ",simplify = T)[,1]
# colnames(ceres_sanger) <-c("DepMap_ID", str_split(colnames(ceres_sanger)[2:17800],
#                                                   n=2,pattern= " ",simplify = T)[,1])
# 
# ceres_columns <- intersect(colnames(ceres_broad), colnames(ceres_sanger))
# sanger_unique <- setdiff(ceres_sanger$DepMap_ID,ceres_broad$DepMap_ID)

# use broad data if available
# ceres <- union(ceres_broad %>% select(one_of(ceres_columns)), 
#                ceres_sanger %>% 
#                  select(one_of(ceres_columns)) %>% 
#                  filter(DepMap_ID %in% sanger_unique)) 
# 
# rm(ceres_broad, ceres_sanger,ceres_columns, sanger_unique)
colnames(ceres)[2:18120] <- str_split(colnames(ceres)[2:18120],n=2,pattern= " ",simplify = T)[,1]


#demeter2
gene <- str_split(demeter2$X1, n=2,pattern= " ",simplify = T)[,1]
demeter2[,1] <- gene
tmp <- colnames(demeter2)
colnames(demeter2)[(tmp %in% "AZ521_STOMACH")] <- "AZ521_SMALL_INTESTINE"
colnames(demeter2)[(tmp %in% "GISTT1_GASTROINTESTINAL_TRACT")] <- "GISTT1_STOMACH"
colnames(demeter2)[(tmp %in% "SW527_BREAST")] <- "SW527_LARGE_INTESTINE"
demeter2 <- demeter2 %>% 
  rename(gene= X1) %>% 
  pivot_longer(cols = -gene,names_to= "CCLE_Name", values_to= "demeter2") %>% 
  inner_join(cell_info %>% select(CCLE_Name, DepMap_ID), by=c("CCLE_Name") ) %>% 
  select(-CCLE_Name) %>% 
  pivot_wider(names_from= "gene", values_from= "demeter2") 
rm(tmp,gene)

#mut
mut <- mut %>% 
  select(one_of( "Hugo_Symbol","DepMap_ID" )) %>% 
  rename(gene=Hugo_Symbol) %>% 
  group_by(gene, DepMap_ID) %>% 
  summarize(mut= n()) %>% 
  pivot_wider(names_from= "gene", values_from= "mut") %>% 
  replace(is.na(.), 0 )

#seq
#duplicate genes in seq "PINX1" "TBCE"
exp_seq_PINX1 <-  exp_seq %>% select(contains("PINX1 (54984)")) %>% 
  transmute(`PINX1 (54984)`= rowMeans(select(., contains("PINX1 (54984)")), na.rm = TRUE))

exp_seq_TBCE <-   exp_seq %>% select(contains("TBCE (6905)")) %>% 
  transmute(`TBCE (6905)`= rowMeans(select(., contains("TBCE (6905)")), na.rm = TRUE))

exp_seq <- exp_seq %>% 
  select(!contains("PINX1 (54984)")) %>% 
  select(!contains("TBCE (6905)")) %>% 
  bind_cols(exp_seq_PINX1) %>% 
  bind_cols(exp_seq_TBCE) 

colnames(exp_seq) <-c("DepMap_ID", str_split(colnames(exp_seq)[2:19181],
                                             n=2,pattern= " ",simplify = T)[,1])

rm(exp_seq_PINX1, exp_seq_TBCE)

#cn
colnames(cn) <-c("DepMap_ID", str_split(colnames(cn)[2:27563],
                                        n=2,pattern= " ",simplify = T)[,1])

#exp_array
colnames(exp_array)[993] <- "NCIH292_LUNG_1"
exp_array$NCIH292_LUNG <- (exp_array$NCIH292_LUNG +exp_array$NCIH292_LUNG_1)/2
exp_array <- exp_array %>% select(-NCIH292_LUNG_1 )
exp_array <- exp_array %>% select(-Name) %>% filter(Description != "")


exp_array1 <- bind_rows((exp_array %>% filter(Description != "TTL") ),
                        exp_array %>% filter(Description == "TTL") %>%
                          group_by(Description) %>%
                          summarise_all(mean, na.rm=T) %>%
                          ungroup())
colnames(exp_array1)[2:1037] <- str_split(colnames(exp_array1)[2:1037], pattern = "_",n = 2,simplify = T)[,1]
colnames(exp_array1)[which(colnames(exp_array1) %in% c("D341MED"))] <- "D341"
colnames(exp_array1)[which(colnames(exp_array1) %in% c("UMRC6"))] <- "UMRC6NEO"
colnames(exp_array1)[which(colnames(exp_array) %in% "TT_OESOPHAGUS") ] <- "TDOTT"

exp_array2 <- exp_array1 %>%
  rename(gene= Description) %>%
  pivot_longer(cols = -gene, names_to= "stripped_cell_line_name", values_to= "exp_array") %>%
  inner_join((cell_info %>% select(stripped_cell_line_name, DepMap_ID)),
             by= "stripped_cell_line_name" ) %>%
  select(-stripped_cell_line_name) %>%
  pivot_wider(names_from = gene, values_from = exp_array)
exp_array <- exp_array2
rm(exp_array1, exp_array2)


#generate ess_long as input to calculate ces
cell_ess <- intersect(ceres$DepMap_ID, unique(demeter2$DepMap_ID)) #474
cells1 <- intersect(cell_ess, exp_seq$DepMap_ID)
cells1 <- intersect(cells1, mut$DepMap_ID)
cells1 <- intersect(cells1, cn$DepMap_ID)
cells2 <- intersect(cells1, exp_array$DepMap_ID) # from 474 to 437
gene_ess <- setdiff(intersect(colnames(ceres), colnames(demeter2)), "DepMap_ID")

wide_to_long <- function(df, value, cell_list){
  df %>% 
    filter(DepMap_ID %in% cell_list) %>% 
    gather(key= "gene", value= !!value, 2:ncol(df)) %>% 
    filter(gene %in% gene_ess)
}

# save an RData for running the following process on server 
save.image("/scratch/project_2003466/ces_21q1_io/ces_prepare_21q1_intermediate.RData")

setwd("/scratch/project_2003466/ces_21q1_io/")
# load("ces_prepare_21q1_intermediate.RData")

df_long_original <- 
  wide_to_long(ceres, "ceres",cells2) %>% 
  inner_join(wide_to_long(demeter2, "demeter2",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_array, "exp_array",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_seq, "exp_seq",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(mut, "mut",cells2) ,by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(cn, "cn",cells2) ,by=c("DepMap_ID", "gene"))%>% 
  drop_na() %>% 
  mutate_all(.funs = as.vector)


df_long_scaled_noarray <- 
  wide_to_long(ceres, "ceres",cells1) %>% 
  inner_join(wide_to_long(demeter2, "demeter2",cells1),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_seq, "exp_seq",cells1),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(mut, "mut",cells1) ,by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(cn, "cn",cells1) ,by=c("DepMap_ID", "gene")) %>% 
  drop_na() %>% 
  group_by(DepMap_ID) %>% 
  mutate_at(.vars = 3:7, .funs = scale) %>% 
  ungroup() %>% 
  mutate_all(.funs = as.vector)


df_long_scaled <- 
  wide_to_long(ceres, "ceres",cells2) %>% 
  inner_join(wide_to_long(demeter2, "demeter2",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_array, "exp_array",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(exp_seq, "exp_seq",cells2),by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(mut, "mut",cells2) ,by=c("DepMap_ID", "gene")) %>% 
  left_join(wide_to_long(cn, "cn",cells2) ,by=c("DepMap_ID", "gene"))%>% 
  drop_na() %>% 
  group_by(DepMap_ID) %>% 
  mutate_at(.vars = 3:8,.funs = scale) %>% 
  ungroup() %>% 
  mutate_all(.funs = as.vector)


write_csv(df_long_original, "ces_input_21q1_original.csv")
write_csv(df_long_scaled, "ces_input_21q1.csv")
write_csv(df_long_scaled_noarray,"ces_input_noarray_21q1.csv")

