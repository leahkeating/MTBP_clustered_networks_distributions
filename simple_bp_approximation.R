########################################################################
# Title: simple_bp_approximation.R
# Purpose: to approximate the network as tree-like and run the dynamics on it
# Author: Leah Keating
# Date last edited: 20 September 2022
########################################################################

rm(list = ls()) # tidy work space
gc()

source("functions.R")

# first, we need a network to extract the distribution from and
# run the dynamics on.

# we chose an Erdos-Renyi network because it should have very few
# triangles so the distributions should match well
net <- sample_gnp(n = 10000, p = 0.0005)
vertex_attr(net)$name <- str_c("v_",1:gorder(net))
adjacency <- as_adj(net)

# next, we get the degree distribution from the network

degree_dist <- net %>% degree() %>% tibble(k = .) %>%
  group_by(k) %>% summarise(n = n()) %>%
  ungroup() %>% mutate(p = n/sum(n))

# plot the degree distribution

degree_dist %>%
  ggplot(aes(x = k, y = p)) +
  geom_point()

# now we need the pgf for the degree distribution and the excess-degree distribution

# degree distribution
f_tilde <- function(x, degree_dist){
  return(sum(degree_dist$p*x^degree_dist$k))
}

# excess-degree distribution
f <- function(x, degree_dist){
  return(sum(degree_dist$k*degree_dist$p*(x^(degree_dist$k - 1))/sum(degree_dist$k*degree_dist$p)))
}

# here we have the iterative function for cascade size (as given in Appendix A)
cascade_size_simple_pgf <- function(z,p,N=10000, dist = degree_dist){
  R <- 1
  for (i in 1:(N-1)) {
    R_new <- 1 - p + p*z*f(R, degree_dist = dist)
    R <- R_new
  }
  R_tilde <- z*f_tilde(R, degree_dist = degree_dist)
  return(R_tilde)
}

# for sanity, check that the following returns a 1
cascade_size_simple_pgf(z = 1, p = 0.01)

# vectorise the pgf function
pgf_cascade_size_simple_vec <- function(x) sapply(x, cascade_size_simple_pgf, p = 0.05, dist = degree_dist)

# invert the pgf using 100 points, use more for more accuracy
theoretical_cascade_size_simple.df <- invert_pgf_via_ifft(pgf_cascade_size_simple_vec, M = 100)
theoretical_cascade_size_simple.df <- theoretical_cascade_size_simple.df %>% filter(cascade_size >0)

theoretical_cascade_size_simple.df %>%
  ggplot(aes(x = cascade_size, y = prob)) +
  geom_line() +
  scale_y_log10() +
  scale_x_log10()

# simulate cascades on this network to check that everything is working as it should.

cluster <- makeCluster(6)
registerDoParallel(cluster)

tic()
cascades.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(j=i, net = net, adj = adjacency, p1 = 0.05, alpha = 0.0, total = 1667)
toc()

stopImplicitCluster()

cascade_size_dist.df <- cascades.df %>% group_by(ID) %>%
  summarise(size = n_distinct(c(parent, child))) %>% group_by(size) %>%
  summarise(n = n()) %>% ungroup()
num_size_1 <- cascade_size_dist.df %>% summarise(1667*6 - sum(n)) %>% as.numeric()
cascade_size_dist.df <- cascade_size_dist.df %>% add_row(size = 1, n = num_size_1) %>% arrange(size)
cascade_size_dist.df <- cascade_size_dist.df %>% mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)

cascade_size_dist.df <- cascade_size_dist.df %>% rename(prob = p, cascade_size = size)

# plot theory with simulations

theoretical_cascade_size_simple.df %>% filter(cascade_size>0) %>%
  ggplot(., aes(x = cascade_size, y = prob)) +
  geom_line() +
  geom_point(data = cascade_size_dist.df) +
  #geom_line(data = theoretical_size_dist_simple.df, colour = "red") +
  scale_y_log10(limits = c(10^(-6),10^(0))) +
  scale_x_log10()
