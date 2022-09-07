########################################################################
# Title: functions.R
# Purpose: script to hold all important functions
# Author: Leah Keting
# Date last edited: 7 September 2022
########################################################################

library(tidyverse)
library(igraph)

projection_from_bipartite_df <- function(vertex_clique.df, nodes){
  bipartite.net <- graph_from_data_frame(vertex_clique.df, directed = FALSE)
  V(bipartite.net)$type <- V(bipartite.net)$name %in% nodes
  return(bipartite.projection(bipartite.net)$proj2)
}

net_gen_double_poisson <- function(n = 100, v = 4, u = 1){
  nodes <- str_c("v_",1:n, sep="")
  node_tri_deg <- rpois(n, lambda = v)
  node_link_deg <- rpois(n, lambda = u)
  tri_deg_sum <- sum(node_tri_deg)
  link_deg_sum <- sum(node_link_deg)
  i = 1
  while (tri_deg_sum %% 3 != 0 | link_deg_sum %% 2 != 0) {
    nodes <- c(nodes, str_c("v_", n+i))
    node_tri_deg <- c(node_tri_deg, rpois(1, v))
    node_link_deg <- c(node_link_deg, rpois(1, u))
    i = i + 1
    tri_deg_sum <- sum(node_tri_deg)
    link_deg_sum <- sum(node_link_deg)
  }
  three_cliques <- str_c("three_",1:(tri_deg_sum/3), sep = "") %>% rep(., each = 3) %>% sample()
  two_cliques <- str_c("two_",1:(link_deg_sum/2), sep = "") %>% rep(., each = 2) %>% sample()
  vertex_clique.df <- tibble(vertex = rep(nodes,node_tri_deg), clique = three_cliques) %>%
    add_row(tibble(vertex = rep(nodes,node_link_deg), clique = two_cliques))
  
  net <- projection_from_bipartite_df(vertex_clique.df, nodes)
}
