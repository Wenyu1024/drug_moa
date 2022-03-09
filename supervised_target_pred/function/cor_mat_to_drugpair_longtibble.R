cor_mat_to_drugpair_longtibble <- function(mat, metric= "similiarity" ){
  # if (metric == "similiarity"){print("Similarity returned") } 
  if (metric == "rank_percentage"){print("Rank percentage returned")
    diag(mat) <- 0
    mat = apply(mat, 2, FUN = function(x){round((rank(x,ties.method = "max")*100)/length(x),digits = 3)})}
  mat[upper.tri(mat,diag = T)] <- -1000
  mat %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "drug1") %>% 
    pivot_longer(cols = -drug1, names_to="drug2", values_to= "value") %>% 
    filter(value != -1000)
}