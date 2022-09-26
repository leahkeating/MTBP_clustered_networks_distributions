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

net <- read_graph("netscience.gml", format = "gml") %>%
  set_vertex_attr("name", value = as.character(1:gorder(.)))
V(net)$name <- as.character(1:gorder(net))
adjacency <- as_adj(net)

# or generate a synthetic network

net <- net_gen_double_poisson(n = 5000, v = 4, u = 1)
adjacency <- as_adj(net)

########

EECC <- EECC_generator(net, max_cl_size = 3)

# check that we only have 2 and 3 cliques

EECC %>% group_by(cl_ID) %>%
  summarise(clique_size = n_distinct(c(vertex_1, vertex_2))) %>% arrange(desc(clique_size)) %>%
  group_by(clique_size) %>% summarise(n())

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
EECC_dGdx <- function(x, y, deg_prob = EECC_deg_prob){
  return(sum((deg_prob$p)*(deg_prob$n_2)*(x^(deg_prob$n_2 - 1))*(y^deg_prob$n_3)))
}
EECC_dGdy <- function(x, y, deg_prob = EECC_deg_prob){
  return(sum((deg_prob$p)*(deg_prob$n_3)*(x^(deg_prob$n_2))*(y^(deg_prob$n_3 - 1))))
}

pgf_size_vec <- function(x) sapply(x, size_pgf, pk = pk, n = 1000, p = 0.02, alpha = 0.2, deg_pgf = EECC_pgf, dGdx = EECC_dGdx, dGdy = EECC_dGdy)

theoretical_size_dist.df <- invert_pgf_via_ifft(pgf_size_vec, M = 100)
theoretical_size_dist.df <- theoretical_size_dist.df %>% filter(cascade_size>0)
theoretical_size_dist.df <- theoretical_size_dist.df %>% mutate(cdf = cumsum(prob), ccdf = 1-cdf)

#######################
# simulate cascades
#######################

cluster <- makeCluster(6)
registerDoParallel(cluster)

tic()
cascades.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(j = i, net = net, adj = adjacency, p = 0.02, alpha = 0.2, total = 166667)
toc()

stopImplicitCluster()

cascade_size_dist.df <- cascades.df %>% group_by(ID) %>%
  summarise(size = n_distinct(c(parent, child))) %>% group_by(size) %>%
  summarise(n = n()) %>% ungroup()
num_size_1 <- cascade_size_dist.df %>% summarise(166667*6 - sum(n)) %>% as.numeric()
cascade_size_dist.df <- cascade_size_dist.df %>% add_row(size = 1, n = num_size_1) %>% arrange(size)
cascade_size_dist.df <- cascade_size_dist.df %>% mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)

#cascade_size_dist.df <- read_csv("powergrid_cascade_size.csv")
#cascade_size_dist.df %>% write_csv("science_coauthorship_cascade_size.csv")
#cascade_size_dist.df %>% write_csv("cascade_size_cc.csv")

#### try SED

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
cascade_size_simple_pgf <- function(z,p,N=1000, dist = degree_dist){
  R <- 1
  for (i in 1:(N-1)) {
    R_new <- 1 - p + p*z*f(R, degree_dist = dist)
    R <- R_new
  }
  R_tilde <- z*f_tilde(R, degree_dist = degree_dist)
  return(R_tilde)
}

pgf_size_simple_vec <- function(x) sapply(x, cascade_size_simple_pgf, N = 1000, p = 0.02, dist = degree_dist)

theoretical_size_dist_simple.df <- invert_pgf_via_ifft(pgf_size_simple_vec, M = 100)

theoretical_size_dist_simple.df <- theoretical_size_dist_simple.df %>% filter(cascade_size>0)

theoretical_size_dist_simple.df <- theoretical_size_dist_simple.df %>%
  mutate(cdf = cumsum(prob), ccdf = 1-cdf)

require(scales)
cascade_size_dist.df %>% rename(cascade_size = size, prob = p) %>% filter(ccdf > 0) %>%
  ggplot(aes(x = cascade_size, y = ccdf)) +
  geom_line(data = theoretical_size_dist.df, colour = "black", size = 1) +
  geom_line(data = theoretical_size_dist_simple.df, colour = "red", size = 1, linetype = "dashed") +
  geom_point(colour = "black", fill = "#fc03ec", size = 2, alpha = 1, 
             shape = 21) +
  scale_y_log10(limits = c(10^(-6),1), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1,10^1.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("ccdf") +
  xlab("cascade size")
