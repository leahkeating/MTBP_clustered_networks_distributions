########################################################################
# Title: EECC_on_empirical_networks.R
# Purpose: to find the theoretical cascade-size distributions for real-
# world networks
# Author: Leah Keating
# Date last edited: 9 September 2022
########################################################################

rm(list = ls()) # tidy work space
gc()

# libraries, functions

source("functions.R")

# load data

net <- read_graph("power.gml", format = "gml") %>%
  set_vertex_attr("name", value = as.character(1:gorder(.)))

# this bit gets the lcc - only needed fo science co-authorship network
lcc_id <- which(components(net)$csize == max(components(net)$csize))
vertices_in_lcc <- V(net)[which(components(net)$membership == lcc_id)]
net <- induced_subgraph(net, vertices_in_lcc)
#----

V(net)$name <- as.character(1:gorder(net))
adjacency <- as_adj(net)

# or generate a synthetic network

net <- net_gen_double_poisson(n = 5000, v = 4, u = 1)
adjacency <- as_adj(net)

########

# get the edge-disjoint edge clique cover
tic()
EECC <- EECC_generator(net, max_cl_size = 3)
toc()
#EECC <- read_csv("synthetic_EECC.csv")
#EECC <- read_csv("net_science_EECC.csv")
#EECC <- read_csv("c_elegans_EECC.csv")
#EECC <- read_csv("power_grid_EECC.csv")

# check that we only have 2 and 3 cliques

EECC %>% group_by(cl_ID) %>%
  summarise(clique_size = n_distinct(c(vertex_1, vertex_2))) %>% arrange(desc(clique_size)) %>%
  group_by(clique_size) %>% summarise(n())

# here we are getting the clique distribution from the EECC, allowing us to find the 
# pgf for the Newman-Miller distribution

EECC <- EECC %>% group_by(cl_ID) %>% mutate(size = n_distinct(c(vertex_1, vertex_2))) %>%
  arrange(desc(size))

EECC_deg_dist <- EECC %>% ungroup() %>% select(-vertex_2) %>% rename(vertex = vertex_1) %>%
  add_row(EECC %>% ungroup() %>% select(-vertex_1) %>% rename(vertex = vertex_2)) %>% unique() %>%
  arrange(desc(size), cl_ID) %>%
  group_by(vertex) %>% summarise(n_3 = sum(size == 3), n_2 = sum(size == 2))

EECC_deg_prob <- EECC_deg_dist %>% #filter(n_3!=0 & n_2!=0) %>%
  group_by(n_3, n_2) %>% summarise(count = n()) %>%
  ungroup() %>% mutate(p = count/sum(count)) # can build the pgf from this

# Get a pgf from this

EECC_pgf <- function(x, y, deg_prob = EECC_deg_prob){
  return(sum((deg_prob$p)*(x^deg_prob$n_2)*(y^deg_prob$n_3)))
}
# need the partial derivatives for the cascade size pgf (and others)
EECC_dGdx <- function(x, y, deg_prob = EECC_deg_prob){
  return(sum((deg_prob$p)*(deg_prob$n_2)*(x^(deg_prob$n_2 - 1))*(y^deg_prob$n_3)))
}
EECC_dGdy <- function(x, y, deg_prob = EECC_deg_prob){
  return(sum((deg_prob$p)*(deg_prob$n_3)*(x^(deg_prob$n_2))*(y^(deg_prob$n_3 - 1))))
}

size_pgf(x = exp(0.325*1i), pk = pk, p = 0.05, alpha = 0.005, deg_pgf = EECC_pgf, dGdx = EECC_dGdx, dGdy = EECC_dGdy, tol = 10^(-5), verbose = TRUE)

pgf_size_vec <- function(x) sapply(x, size_pgf, tol = 10^(-5),pk = pk, p = 0.05, alpha = 0.0, deg_pgf = EECC_pgf, dGdx = EECC_dGdx, dGdy = EECC_dGdy)

theoretical_size_dist.df <- invert_pgf_via_ifft(pgf_size_vec, M = 100)
theoretical_size_dist.df <- theoretical_size_dist.df %>% filter(cascade_size>0)
theoretical_size_dist.df <- theoretical_size_dist.df %>% mutate(cdf = cumsum(prob), ccdf = 1-cdf)

#######################
# simulate cascades
#######################

cluster <- makeCluster(6)
registerDoParallel(cluster)

# simulation function to generate cascades
tic()
cascades.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(j = i, net = net, adj = adjacency, p = 0.05, alpha = 0.0, total = 166667)
toc()

stopImplicitCluster()

# find the distribution of cascade sizes from simulations
cascade_size_dist.df <- cascades.df %>% group_by(ID) %>%
  summarise(size = n_distinct(c(parent, child))) %>% group_by(size) %>%
  summarise(n = n()) %>% ungroup()
num_size_1 <- cascade_size_dist.df %>% summarise(166667*6 - sum(n)) %>% as.numeric()
cascade_size_dist.df <- cascade_size_dist.df %>% add_row(size = 1, n = num_size_1) %>% arrange(size)
cascade_size_dist.df <- cascade_size_dist.df %>% mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)

#cascade_size_dist.df <- read_csv("powergrid_cascade_size.csv")
#cascade_size_dist.df <- read_csv("science_coauthorship_cascade_size.csv")
#cascade_size_dist.df <- read_csv("cascade_size_cc.csv")
#cascade_size_dist.df <- read_csv("c_elegans_cascade_size.csv")
#cascade_size_dist.df <- read_csv("cascade_size_sc.csv")

#### try SED (for comparison)

# instead of the clique distribution we get the degree distribution
degree_dist <- net %>% degree() %>% tibble(k = .) %>%
  group_by(k) %>% summarise(n = n()) %>%
  ungroup() %>% mutate(p = n/sum(n))

# degree distribution
f_tilde <- function(x, degree_dist){
  return(sum(degree_dist$p*x^degree_dist$k))
}

# excess-degree distribution
f <- function(x, degree_dist){
  return(sum(degree_dist$k*degree_dist$p*(x^(degree_dist$k - 1))/sum(degree_dist$k*degree_dist$p)))
}

# here we have the iterative function for cascade size (as given in Appendix A)
cascade_size_simple_pgf <- function(z,p,tol = 10^(-5), dist = degree_dist, verbose = FALSE, max_n = 10000){
  R <- 1
  R_old <- R
  R_new <- 0
  N <- 0
  while(Mod((R_new - R_old))>tol) {
    N <- N+1
    R_new <- 1 - p + p*z*f(R, degree_dist = dist)
    R_old <- R
    R <- R_new
    if(N>max_n) {
      print("max iterations reached")
      break
      }
  }
  R_tilde <- z*f_tilde(R, degree_dist = degree_dist)
  if(verbose) print(N)
  return(R_tilde)
}

# just checking that it works
cascade_size_simple_pgf(z = 0.1, p = 0.005, dist = degree_dist, verbose = TRUE)

# create a function that takes in a vector of values and evaluates the pgf at those values
pgf_size_simple_vec <- function(x) sapply(x, cascade_size_simple_pgf, p = 0.05, dist = degree_dist)

# this function evaluates the pgf at M evenly spaced points around the unit circle in the complex plane,
# inverting it to return the distribution
theoretical_size_dist_simple.df <- invert_pgf_via_ifft(pgf_size_simple_vec, M = 100)

theoretical_size_dist_simple.df <- theoretical_size_dist_simple.df %>% filter(cascade_size>0)

theoretical_size_dist_simple.df <- theoretical_size_dist_simple.df %>%
  mutate(cdf = cumsum(prob), ccdf = 1-cdf)

theoretical_size_dist.df <- theoretical_size_dist.df %>% filter(ccdf>0 & prob>0)
theoretical_size_dist_simple.df <- theoretical_size_dist_simple.df %>% filter(ccdf>0 & prob>0)
require(scales)
cascade_size_dist.df %>% rename(cascade_size = size, prob = p) %>% filter(ccdf > 0) %>%
  ggplot(aes(x = cascade_size, y = prob)) +
  geom_line(data = theoretical_size_dist.df, colour = "black", size = 1) +
  geom_line(data = theoretical_size_dist_simple.df, colour = "red", size = 1, linetype = "dashed") +
  geom_point(colour = "black", fill = NA, size = 3, alpha = 1, 
             shape = 21) +
  scale_y_log10(limits = c(10^(-6),1), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1,10^1.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("probability") +
  xlab("cascade size")

######################################################################################
# Build a larger network from the distribution of the science co-authorship network
# (and others) as approximated by the EECC
# this gives us the plots in Figure C.1.
######################################################################################

EECC_deg_prob

# generate a network from a given Newman-Miller distribution
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

new_net_same_size <- net_gen_joint_dist(n = 4941, joint_dist = EECC_deg_prob)
same_size_adj <- as_adj(new_net_same_size)
# the one below is to make a larger network in the case of small networks
# where finite-size effects might be substantial
new_net_10000 <- net_gen_joint_dist(n = 1000, joint_dist = EECC_deg_prob)
adj_10000 <- as_adj(new_net_10000)
new_net %>% transitivity()

cluster <- makeCluster(6)
registerDoParallel(cluster)
# simulate cascades on these new networks
tic()
new_net_small_cascades.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(j = i, net = new_net_same_size, adj = same_size_adj, p = 0.04, alpha = 0.15, total = 166667)
toc()
tic()
new_net_large_cascades.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(j = i, net = new_net_10000, adj = adj_10000, p = 0.005, alpha = 0.005, total = 166667)
toc()

stopImplicitCluster()

# tidy this simulation information into their distributions
cascade_size_dist_small.df <- new_net_small_cascades.df %>% group_by(ID) %>%
  summarise(size = n_distinct(c(parent, child))) %>% group_by(size) %>%
  summarise(n = n()) %>% ungroup()
num_size_1 <- cascade_size_dist_small.df %>% summarise(166667*6 - sum(n)) %>% as.numeric()
cascade_size_dist_small.df <- cascade_size_dist_small.df %>% add_row(size = 1, n = num_size_1) %>% arrange(size)
cascade_size_dist_small.df <- cascade_size_dist_small.df %>% mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)
cascade_size_dist_small.df <- cascade_size_dist_small.df %>% rename(cascade_size = size, prob = p) %>% filter(ccdf > 0)

cascade_size_dist_large.df <- new_net_large_cascades.df %>% group_by(ID) %>%
  summarise(size = n_distinct(c(parent, child))) %>% group_by(size) %>%
  summarise(n = n()) %>% ungroup()
num_size_1 <- cascade_size_dist_large.df %>% summarise(166667*6 - sum(n)) %>% as.numeric()
cascade_size_dist_large.df <- cascade_size_dist_large.df %>% add_row(size = 1, n = num_size_1) %>% arrange(size)
cascade_size_dist_large.df <- cascade_size_dist_large.df %>% mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)
cascade_size_dist_large.df <- cascade_size_dist_large.df %>% rename(cascade_size = size, prob = p) %>% filter(ccdf > 0)

cascade_size_dist.df <- cascade_size_dist.df %>% rename(cascade_size = size, prob = p) %>% filter(ccdf > 0)

# for plotting label the networks that the simulations have come from
cascade_size_all_dists.df <- cascade_size_dist.df %>% mutate(network = "empirical") %>%
  add_row(cascade_size_dist_small.df%>% mutate(network = "4947 nodes")) #%>%
  add_row(cascade_size_dist_large.df%>% mutate(network = "1,000 nodes")) %>% arrange((network))

#cascade_size_all_dists.df <- read_csv("cascade_size_all_dists_c_elegans.csv")
#cascade_size_all_dists.df <- read_csv("cascade_size_all_dists_coauthorship.csv")

ggplot(data = NULL,aes(x = cascade_size, y = prob)) +
geom_line(data = theoretical_size_dist.df, colour = "black", size = 1) +
geom_line(data = theoretical_size_dist_simple.df, colour = "red", size = 1, linetype = "dashed") +
geom_point(data = cascade_size_all_dists.df, aes(colour = network), size = 2) +
scale_y_log10(limits = c(10^(-6),10^(0)), breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
scale_x_log10(limits = c(1,10^1.5), breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
ylab("probability") +
xlab("cascade size") +
theme(legend.text = element_text(size = 10), legend.position = c(0.5,0.8),
      legend.box.background = element_rect(fill='white'),
      legend.margin = margin(6, 6, 6, 6),
      legend.title = element_text(size = 11)) +
  scale_color_manual(values = c("#d95f02","#7570b3"))
  #scale_colour_brewer(palette = "Dark2")

########################################################
# Check that the clique cover is giving what we expect
########################################################

EECC_dist <- EECC_deg_prob %>% mutate(deg = 2*n_3 + n_2) %>%
  group_by(deg) %>% summarise(count = sum(count)) %>%
  ungroup() %>% mutate(prob = count/sum(count))

EECC_dist %>%
  ggplot(aes(x = deg, y = prob)) +
  geom_point()

