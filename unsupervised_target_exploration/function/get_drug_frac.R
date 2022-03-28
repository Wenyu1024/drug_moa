library(tidyverse)
get_drug_frac <- function(dataset, top_num, res_type = "target_summary", top_method="top_abs"){
  # The function calculates that for a given dataset, how many drug, the target of which is prioritized? 
  # top drug is selected by abosoulte values (top_abs) or both side(top_both)
  
  n_drug= dataset %>% filter(anno_type1 %in%  "Target") %>% select(drug) %>% distinct() %>% nrow()
  
  # target
  if(top_method == "top_both") {
    top_druggenepairs<-bind_rows(
      dataset %>%
        group_by(drug,imptype) %>% 
        top_n(wt=imp,  n= top_num) %>% 
        ungroup() %>% 
        mutate(impdef= "toppos gene"),
      dataset %>% 
        group_by(drug, imptype) %>% 
        top_n(wt=imp, n = -top_num) %>% 
        ungroup() %>% 
        mutate(impdef= "topneg gene"),
    )
  }
  
  
  if(top_method == "top_abs") {
    top_druggenepairs<-
      dataset %>%
      mutate(imp_abs= abs(imp)) %>% 
      group_by(drug,imptype) %>% 
      top_n(wt=imp_abs,  n= top_num) %>% 
      ungroup() %>%
      mutate(impdef= case_when(
        imp>0 ~"toppos gene",
        imp<0 ~"topneg gene"))
    
  }
  
  
  drug_shortest_path_tibble <- top_druggenepairs %>% 
    mutate(anno_type2 = case_when(
      anno_type!="Target" ~anno_type,
      TRUE ~ "0")) %>%    
    mutate(anno_type2= as.numeric(anno_type2) ) %>% 
    drop_na() %>% 
    group_by(drug, imptype) %>% 
    top_n(wt = anno_type2,n = -1) %>% 
    ungroup() %>% 
    select(imptype,drug, anno_type2,impdef) %>% 
    distinct(imptype,drug, anno_type2,.keep_all = TRUE)
  
  target_summary <- 
    top_druggenepairs %>% 
    filter(anno_type1 %in% "Target") %>% 
    select(drug, imptype, impdef) %>% 
    distinct() %>% 
    group_by(imptype,impdef) %>% 
    summarize(frac= round(n()*100/n_drug, 1),
              count= n()) %>% 
    mutate(top_n= top_num)
  
  
  if (res_type == "target_summary") {return(target_summary)} 
  if (res_type == "toppair") {return(top_druggenepairs)} 
  if (res_type == "druglevel") {return(drug_shortest_path_tibble)} 
  
}