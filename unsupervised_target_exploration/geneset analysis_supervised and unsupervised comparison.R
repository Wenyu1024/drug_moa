library(tidyverse)
library(furrr)
library(fgsea)

setwd("~/cluster_wrk/drug_moa")
#Step 1 generate gmt files based on the top genes from supervised
#and unsupervised analysis
# basically for supervised
pred_58_unsupervised <- read_csv("pred_58_unsupervised.csv")
pred_58_supervised <- read_csv("pred_58_supervised.csv")

# a function to generate a gmt object from a pred wide tibble
get_gmt <- function(pred_tibble, top_num){
  tmp <- pred_tibble%>% 
    pivot_longer(cols=-target, names_to = "drug", values_to= "pred") %>% 
    group_by(drug) %>% 
    mutate(gene_rank= rank(-pred,ties.method = "first")) %>% 
    top_n(n = top_num+1,wt = -gene_rank) %>% 
    arrange(gene_rank) %>% 
    ungroup() %>% 
    select(-pred) %>% 
    pivot_wider(names_from = gene_rank, values_from= target,names_prefix = "gene") 
  return(tmp)  
}

return_pathway_for_genefeaturematrix_par <- function(data,pathways){
  # res_pathway <- vector(length = nrow(data), mode = "list")
  res_pathway <- future_map(
    .options = furrr_options(seed = 1),
    .x = 1:(nrow(data)),
    .f = function(x){  
      rank <- (unlist(data[x,-1 ]))  # use the absolute value for the rank.
      fgseaRes <- fgsea(pathways, rank, minSize=15, 
                        maxSize=1000,scoreType="pos")
      return(fgseaRes)})
  return(res_pathway)
}

pathway_sup_30 <- get_gmt(pred_tibble = pred_58_supervised,top_num = 30)
gmt::r2gmt(pathway_sup_30, "pathway_sup_30.gmt")
pathway_sup_30  <- gmtPathways("pathway_sup_30.gmt")

pathway_sup_50 <- get_gmt(pred_tibble = pred_58_supervised,top_num = 50)
gmt::r2gmt(pathway_sup_50, "pathway_sup_50.gmt")
pathway_sup_50  <- gmtPathways("pathway_sup_50.gmt")

pathway_sup_100 <- get_gmt(pred_tibble = pred_58_supervised,top_num = 100)
gmt::r2gmt(pathway_sup_100, "pathway_sup_100.gmt")
pathway_sup_100  <- gmtPathways("pathway_sup_100.gmt")

pathway_sup_150 <- get_gmt(pred_tibble = pred_58_supervised,top_num = 150)
gmt::r2gmt(pathway_sup_150, "pathway_sup_150.gmt")
pathway_sup_150  <- gmtPathways("pathway_sup_150.gmt")

 # test association 
# fgseaRes <- fgsea(pathway_sup_30, rank, minSize=15, maxSize=500)

data_unsupervised <- pred_58_unsupervised %>% 
  pivot_longer(cols = -target,names_to="drug", values_to="value") %>% 
  pivot_wider(id_col=drug, names_from = target, values_from=value)

tmp <- 
  bind_cols(
    return_pathway_for_genefeaturematrix_par(
      data = data_unsupervised,
      pathways = pathway_sup_150) %>% 
      tibble()
    ,
    data_unsupervised %>% select(drug)
    ) %>% 
  unnest(cols = c(".")) #%>% 
  # select(pathway, drug, pval, padj) 

tmp$id <- tmp$pathway == tmp$drug
tmp1 <- tmp %>%   filter(id)

hist(tmp1$ES,main = "",xlab = "Enrichment score")

# get_ran_gmt <- function(gene_list= colnames(data_supervised)[-1], 
#                         pathway_num=10, 
#                         pathway_size=150){
#   set.seed(1)
#   pathway_tmp <- map(
#     .x = 1:pathway_num, 
#     .f = ~sample(x = gene_list,size = pathway_size))
#   # pathway_tmp <- matrix(unlist(pathway_tmp),nrow = pathway_num,ncol = pathway_size,byrow = T)
#   
#   names(pathway_tmp) <- str_c("pathway_random",1:pathway_num,sep = "")
#   return(pathway_tmp)
# }
# 
# pathway_ran_sup150_10000 <- get_ran_gmt(pathway_num = 1000)
# tmp2 <- 
#   bind_cols(
#     return_pathway_for_genefeaturematrix_par(
#       data = data_unsupervised,
#       pathways = pathway_ran_sup150_10000) %>% 
#       tibble()
#     ,
#     data_unsupervised %>% select(drug)
#   ) %>% 
#   unnest(cols = c(".")) %>% 
#   select(pathway, drug, pval, padj) 


# what if I randomly generate ranked list?
set.seed(0)
data_unsupervised_random <- map(
  .x = 1:1000, 
  .f = ~sample(x = 1:10623,size = 10623)
  )

data_unsupervised_random <- 
  matrix(unlist(data_unsupervised_random),
         nrow = 1000,
         ncol = 10623,
         byrow = T
        ) %>% 
  as_tibble()

data_unsupervised_random <- bind_cols(drug=str_c("drug_random",1:1000,sep = ""),data_unsupervised_random)
colnames(data_unsupervised_random) <- colnames(data_unsupervised)
tmp3 <- 
  bind_cols(
    return_pathway_for_genefeaturematrix_par(
      data = data_unsupervised_random,
      pathways = pathway_sup_50) %>% 
      tibble()
    ,
    data_unsupervised_random %>% select(drug)
  ) %>% 
  unnest(cols = c(".")) #%>% 
  # select(pathway, drug,pval, padj) 


# ks.test(tmp1$pval,tmp$pval)
# ks.test(tmp1$pval,tmp2$pval)
ks.test(tmp1$pval,tmp3$pval)
# ks.test(tmp2$pval,tmp3$pval)

# hist(tmp$pval)
# summary(tmp$pval)
par(mar=c(1,1,1,1))
tiff("~/cluster_wrk/drug_moa/unsupervised_target_exploration/fig_res/figs1.tif",
    width = 10,height = 10,res = 300,units = 'cm')
hist(tmp1$pval,xlab = "p-value",xlim = c(0,1),
     ylim= c(0,25),main = "",breaks = 20)
dev.off()
# fig1 <- tmp1 %>% ggplot() +
#   geom_histogram(aes(x= pval),bins = 10,binwidth = 0.1, )+
# 
#   theme_classic()+
#   xlab("P-value")
# fig1
summary(tmp1$pval)

# hist(tmp2$pval)
# summary(tmp2$pval)

par(mar=c(1,1,1,1))
tiff("~/cluster_wrk/drug_moa/unsupervised_target_exploration/fig_res/figs2.tif",
     width = 10,height = 10,res = 300,units = 'cm')
hist(tmp3$pval,xlab = "p-value",xlim = c(0,1),main = "",breaks = 20)
dev.off()
# fig2 <- tmp3 %>% ggplot(aes(x= pval)) +
#   geom_histogram(binwidth = 0.05,colour = "black", fill = "grey",stat="bin")+
#   scale_x_continuous(limits=c(0,1),breaks = seq(0, 1, by = 0.05))+
#   # theme_classic()+
#   xlab("P-value")
# fig2
summary(tmp3$pval)


# If I want to do form rank list using the supervised target prediction,
# I should also add gene which are not among the 300 candidate genes in the
# list. just assign them to the bottom (instead of random or middle)

# pathway_unsup_100 <- get_gmt(pred_tibble = pred_58_unsupervised,top_num = 100)
# gmt::r2gmt(pathway_unsup_100, "pathway_unsup_100.gmt")
# pathway_unsup_100  <- gmtPathways("pathway_unsup_100.gmt")
# 
# pathway_unsup_150 <- get_gmt(pred_tibble = pred_58_unsupervised,top_num = 150)
# gmt::r2gmt(pathway_unsup_150, "pathway_unsup_150.gmt")
# pathway_unsup_150  <- gmtPathways("pathway_unsup_150.gmt")
# 
# pathway_unsup_200 <- get_gmt(pred_tibble = pred_58_unsupervised,top_num = 200)
# gmt::r2gmt(pathway_unsup_200, "pathway_unsup_200.gmt")
# pathway_unsup_200  <- gmtPathways("pathway_unsup_200.gmt")
# 
# pathway_unsup_250 <- get_gmt(pred_tibble = pred_58_unsupervised,top_num = 250)
# gmt::r2gmt(pathway_unsup_250, "pathway_unsup_250.gmt")
# pathway_unsup_250  <- gmtPathways("pathway_unsup_250.gmt")
# 
# pathway_unsup_500 <- get_gmt(pred_tibble = pred_58_unsupervised,top_num = 500)
# gmt::r2gmt(pathway_unsup_500, "pathway_unsup_500.gmt")
# pathway_unsup_500  <- gmtPathways("pathway_unsup_500.gmt")
# 
# data_supervised <- pred_58_supervised %>% 
#   pivot_longer(cols = -target,names_to="drug", values_to="value") %>% 
#   pivot_wider(id_col=drug, names_from = target, values_from=value)
# 
# data_supervised %>% select(-1) %>% min() -0.93
# 
# data_supervisedcomplement <- data_unsupervised %>% 
#   select(!any_of(c(colnames(data_supervised)))) %>% 
#   mutate_all(.funs = function(x){(0/x)-1})
# 
# data_supervisedcomplemented <- bind_cols(data_supervised, data_supervisedcomplement)
# 
# tmp1 <-
#   # tmp %>%
#   bind_cols(
#     return_pathway_for_genefeaturematrix_par(
#       data = data_supervisedcomplemented,
#       pathways = pathway_unsup_250) %>%
#       tibble()
#     ,
#     data_supervisedcomplemented %>% select(drug)
#   ) %>%
#   unnest(cols = c(".")) %>%
#   select(pathway, drug, pval, padj)
# 
# tmp1$id <- tmp1$pathway == tmp1$drug
# tmp1 <- tmp1 %>%   filter(id)

# lets do a permutation!
# (randomly generate pathways (150 genes) from the 300 genes)



