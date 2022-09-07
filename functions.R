########################################################################
# Title: functions.R
# Purpose: script to hold all important functions
# Author: Leah Keting
# Date last edited: 7 September 2022
########################################################################

library(tidyverse)
library(igraph)
library(doParallel)
library(tictoc)

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

# This function is used to simulate cascades

generate_cc_cascades <- function(follower.net = follower.net, follower.adj = follower.adj, p1 = 0.002, alpha = 0.002, total = 1000){#, seed_ = NULL
  all_cascades.df <- tibble(parent = character(), child = character(), generation = numeric(), ID = numeric(), exposures = numeric())
  for (j in 1:total) {
    active <- numeric()
    inactive <- numeric()
    removed <- numeric()
    
    vertex_names <- vertex_attr(follower.net)$name
    seed <- sample(vertex_names,1)
    active <- seed
    inactive <- vertex_names[! vertex_names %in% seed]
    
    exposures <- numeric(gorder(follower.net))
    names(exposures) <- vertex_names
    cascade.df <- tibble(parent = character(), child = character(), generation = numeric())
    generation <- 1
    while (length(active)>0) {
      new_active <- character()
      # shuffle active
      if(length(active)>1){
        active <- sample(active)
      }
      for (i in active) {
        followers <- vertex_names[which(follower.adj[,i]==1)]
        potential_adopters <- followers[followers %in% inactive]
        exposures[potential_adopters] <- exposures[potential_adopters] + 1
        if(length(potential_adopters)>0){
          # fix this, problem is with having n a vector
          adopters <- potential_adopters[runif(length(potential_adopters)) < pk(p = p1, alpha = rep(alpha, length(potential_adopters)), k = exposures[potential_adopters])]
          if(length(adopters)>0){
            new_active <- c(new_active, adopters)
            inactive <- inactive[! inactive %in% new_active]
            cascade.df <- cascade.df %>% add_row(parent = rep(i, length(adopters)), child = adopters, generation = rep(generation, length(adopters)))
          }
        }
      }
      generation <- generation + 1
      removed <- c(removed, active)
      active <- new_active
    }
    if(nrow(cascade.df)>0){
      all_cascades.df <- all_cascades.df %>% add_row(cascade.df %>% mutate(ID = rep(j, nrow(cascade.df)), exposures = exposures[cascade.df$child]))
    }
  }
  return(all_cascades.df)
}

# function for simulating the cascades in parallel
cascade_sim_par <- function(j, net, adj, p1, alpha, total = 1000){
  out.df <- generate_cc_cascades(follower.net = net, follower.adj = adj, p1=p1, alpha = alpha, total = total)
  out.df <- out.df %>% mutate(ID = str_c(j,ID,sep = "_"))
  return(out.df)
}
