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
library(cowplot)
library(latex2exp)
theme_set(theme_cowplot())

# want to calculate the cascade size distribution from a single clique
# pgf for each clique type
pk <- function(k, p = 0.05, alpha = 0.05) 1 - ( (1 - p) * (1-alpha)^(k-1))

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

# to find the EECC

EECC_generator <- function(net, max_cl_size = NaN){
  net.el <- net %>% as_edgelist()
  net.el <- tibble(vertex_1 = net.el[,1], vertex_2 = net.el[,2])
  
  C <- max_cliques(net)
  a <- length(C) # this is just to make sure all cliques have different IDs
  
  cliques <- tibble(vertex_1 = character(), vertex_2 = character(), cl_ID = NA)
  
  EECC <- tibble(vertex_1 = character(), vertex_2 = character(), cl_ID = NA)
  
  # make a data frame with the cliques as edgelists
  
  for (i in 1:length(C)) {
    vertices <- C[[i]]$name
    cliques <- cliques %>% add_row(net.el %>% filter(vertex_1 %in% vertices & vertex_2 %in% vertices) %>% unique() %>%
                                     mutate(cl_ID = i))
  }
  
  cliques <- cliques %>% group_by(vertex_1, vertex_2) %>% mutate(num_cls = n())
  
  if(is.na(max_cl_size) == FALSE){
    cliques <- cliques %>% group_by(cl_ID) %>% mutate(size = n_distinct(c(vertex_1, vertex_2)))
  }
  
  # concern: can the edge A - B appear as B - A, it seems that this doesn't happen
  
  overlap_IDs <- cliques[cliques$num_cls>1,]$cl_ID %>% unique()
  if(is.na(max_cl_size) == FALSE){
    higher_order_IDs <- cliques[cliques$size > max_cl_size,]$cl_ID %>% unique()
    overlap_IDs <- overlap_IDs[!overlap_IDs %in% higher_order_IDs]
  }else{
    higher_order_IDs <- numeric()
  }
  
  EECC <- EECC %>% add_row(cliques %>% filter((!cl_ID %in% overlap_IDs) & (!cl_ID %in% higher_order_IDs)) %>% select(c(vertex_1, vertex_2, cl_ID)))
  
  # cliques is C in the algorithm
  
  ordinary_cliques <- cliques %>% filter(cl_ID %in% overlap_IDs)
  if(is.na(max_cl_size) == FALSE){
    higher_order_cliques <- cliques %>% filter(cl_ID %in% higher_order_IDs)
  }
  
  j = 0
  while (nrow(ordinary_cliques)>0) {
    j <- j + 1
    c_net <- ordinary_cliques %>% graph_from_data_frame(directed = FALSE)
    # find the new cliques
    C <- max_cliques(c_net)
    
    ordinary_cliques <- tibble(vertex_1 = character(), vertex_2 = character(), cl_ID = NA)
    
    for (i in 1:length(C)) {
      vertices <- C[[i]]$name
      # the 100*j should be changed some way (if j > 100)
      ordinary_cliques <- ordinary_cliques %>% add_row(net.el %>% filter(vertex_1 %in% vertices & vertex_2 %in% vertices) %>% unique() %>%
                                                         mutate(cl_ID = (a*j + i)))
    }
    
    ordinary_cliques <- ordinary_cliques %>% group_by(vertex_1, vertex_2) %>% mutate(num_cls = n(), edge_id = cur_group_id()) %>% ungroup()
    
    clique_to_add_ordered <- ordinary_cliques %>% group_by(cl_ID) %>% mutate(num_edges = n(), shared_edge = num_cls>1, shared_edge_num = sum(shared_edge), rho = shared_edge_num/num_edges) %>%
      arrange(rho, desc(num_edges)) %>% select(vertex_1, vertex_2, cl_ID, edge_id)
    clique_to_add <- clique_to_add_ordered %>% filter(cl_ID == clique_to_add_ordered$cl_ID[1])
    EECC <- EECC %>% add_row(clique_to_add %>% select(c(vertex_1, vertex_2, cl_ID)))
    ordinary_cliques <- ordinary_cliques %>% filter(cl_ID != clique_to_add_ordered$cl_ID[1], !edge_id %in% clique_to_add$edge_id)
  }
  old_j <- j
  if(is.na(max_cl_size) == FALSE){
    # remove higher-order edges that are already in the EECC
    higher_order_cliques <- higher_order_cliques %>% anti_join(EECC, by = c("vertex_1", "vertex_2"))
    j = 0
    while (nrow(higher_order_cliques)>0) {
      j <- j + 1
      c_net <- higher_order_cliques %>% graph_from_data_frame(directed = FALSE)
      # find the new cliques
      C <- max_cliques(c_net)
      
      higher_order_cliques <- tibble(vertex_1 = character(), vertex_2 = character(), cl_ID = NA)
      
      for (i in 1:length(C)) {
        vertices <- C[[i]]$name
        # the 100*j should be changed some way (if j > 100)
        higher_order_cliques <- higher_order_cliques %>% add_row(net.el %>% filter(vertex_1 %in% vertices & vertex_2 %in% vertices) %>% unique() %>%
                                                                   mutate(cl_ID = (a*(j+old_j) + i)))
      }
      
      higher_order_cliques <- higher_order_cliques %>% group_by(vertex_1, vertex_2) %>% mutate(num_cls = n(), edge_id = cur_group_id()) %>% ungroup()
      
      clique_list <- higher_order_cliques %>% 
        group_by(cl_ID) %>% 
        mutate(num_edges = n(), size = n_distinct(c(vertex_1, vertex_2)), shared_edge = num_cls>1, shared_edge_num = sum(shared_edge), rho = shared_edge_num/num_edges) %>%
        arrange(rho, desc(num_edges))
      cliques_below_max <- clique_list %>% filter(size <= max_cl_size) %>%
        select(vertex_1, vertex_2, cl_ID, edge_id, size)
      if(nrow(cliques_below_max)>0){
        clique_to_add <- cliques_below_max %>% filter(cl_ID == cliques_below_max$cl_ID[1])
      }else{
        full_clique_to_add <- clique_list %>% filter(cl_ID == clique_list$cl_ID[1])
        clique_to_add_vertices <- c(full_clique_to_add$vertex_1, full_clique_to_add$vertex_2) %>% unique() %>% sample(size = max_cl_size)
        clique_to_add <- full_clique_to_add %>% filter(vertex_1 %in% clique_to_add_vertices & vertex_2 %in% clique_to_add_vertices)
      }
      higher_order_cliques <- higher_order_cliques %>% filter(!edge_id %in% clique_to_add$edge_id)
      EECC <- EECC %>% add_row(clique_to_add %>% select(c(vertex_1, vertex_2, cl_ID)))
    }
  }
  return(EECC)
}

# cascade size pgf

size_pgf <- function(x, pk, alpha, p, deg_pgf, dGdx, dGdy, tol = 10^(-4), verbose = FALSE, max_n = 1000){
#  if(n>0){
    p1 <- pk(k = 1, alpha = alpha, p = p)
    p2 <- pk(k = 2, alpha = alpha, p = p)
    # offspring dist from triangle
    g_r <- function(x,y) dGdy(x=x,y=y)/dGdy(x=1,y=1)
    # offspring dist from link
    g_q <- function(x,y) dGdx(x=x,y=y)/dGdx(x=1,y=1)
    # initial conditions
    G_prev <- c(1,1,1) # zero generations to grow
    G_new <- numeric(3)
    G_prev_old <- G_prev
    n <- 0
    while(Mod(sum(G_new - G_prev_old))>tol) {
      n <- n + 1
      #G_new <- numeric(3)
      G_new[1] <- (1-p1)^2 +2*p1*(1-p1)*x*G_prev[2]*g_r(G_prev[3],G_prev[1]) + (p1^2)*(x^2)*(g_r(G_prev[3],G_prev[1])^2)
      G_new[2] <- 1 - p2 + p2*x*g_r(G_prev[3],G_prev[1])
      G_new[3] <- 1-p1 + p1*x*g_q(G_prev[3],G_prev[1])
      G_prev_old <- G_prev
      G_prev <- G_new
      if (n > max_n) print("max iterations reached")
    }
    G_tilde <- x*deg_pgf(x = G_new[3],y = G_new[1])
#  }else if(n == 0) G_tilde <- x
  if(verbose) print(n)
  return(G_tilde)
}

# 1-D inversion

invert_pgf_via_ifft <- function(gen_fn, M = 10^5){
  x <- exp(2*pi*1i*(0:(M-1))/M)
  
  pdf <- Re(fft(gen_fn(x) %>% Conj(),inverse = TRUE))
  pdf[pdf < 0] <- 0
  pdf <- pdf/length(pdf)
  
  cas_dist <- tibble(cascade_size = 0:(length(pdf) - 1), prob = pdf) 
  
  return(cas_dist)
}

# generate a Newman-Miller network from any distribution

net_gen_joint_dist <- function(n = 100, joint_dist){
  nodes <- str_c("v_",1:n, sep="")
  node_deg_indices <- sample(1:nrow(joint_dist),size = n, replace = TRUE, prob = joint_dist$p)
  node_link_deg <- joint_dist$n_2[node_deg_indices]
  node_tri_deg <- joint_dist$n_3[node_deg_indices]
  tri_deg_sum <- sum(node_tri_deg)
  link_deg_sum <- sum(node_link_deg)
  i = 1
  while (tri_deg_sum %% 3 != 0 | link_deg_sum %% 2 != 0) {
    nodes <- c(nodes, str_c("v_", n+i))
    new_deg_index <- sample(1:nrow(joint_dist),size = 1, replace = TRUE, prob = joint_dist$p)
    node_tri_deg <- c(node_tri_deg, joint_dist$n_3[new_deg_index])
    node_link_deg <- c(node_link_deg, joint_dist$n_2[new_deg_index])
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
