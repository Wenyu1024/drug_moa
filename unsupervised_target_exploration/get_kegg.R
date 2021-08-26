# load("C:/Users/wenyu/Downloads/KG.Rdata")
kegg <- igraph::as_long_data_frame(KEGGgraph)
kegg <- kegg[,3:4]
colnames(kegg) <- c("gene1", "gene2")