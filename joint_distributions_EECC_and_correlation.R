########################################################################
# Title: joint_distributions_EECC_and_correlation.R
# Purpose: finding the joint distribution for cascade size and cumulative
# depth, allowing us to calculate the EATD and Pearson's correlation.
# Author: Leah Keating
# Date last edited: 7 September 2022
########################################################################

source("functions.R")

# Newman uses a doubly poisson degree distribution in his paper
# v is the mean number of triangles per vertex (dummy variable y)
# u is the mean number of single links per vertex (dummy variable x)

double_poisson_pgf <- function(x, y, v = 4, u = 1) exp(u*(x-1) + v*(y - 1))
pois_dGds1 <- function(x, y, v = 4, u = 1) v*double_poisson_pgf(x, y, v, u) # partial derivative wrt s1
pois_dGds3 <- function(x, y, v = 4, u = 1) u*double_poisson_pgf(x, y, v, u) # partial derivative wrt s3

#########################################################################################
# For real-world networks we need to infer the Newman-Miller distribution using the EECC
#########################################################################################

net <- read_graph("netscience.gml", format = "gml") %>%
  set_vertex_attr("name", value = as.character(1:gorder(.)))
V(net)$name <- as.character(1:gorder(net))
adjacency <- as_adj(net)

EECC <- EECC_generator(net, max_cl_size = 3)

EECC %>% group_by(cl_ID) %>%
  summarise(clique_size = n_distinct(c(vertex_1, vertex_2))) %>% arrange(desc(clique_size)) %>%
  group_by(clique_size) %>% summarise(n())

EECC <- EECC %>% group_by(cl_ID) %>% mutate(size = n_distinct(c(vertex_1, vertex_2))) %>%
  arrange(desc(size))

EECC_deg_dist <- EECC %>% ungroup() %>% select(-vertex_2) %>% rename(vertex = vertex_1) %>%
  add_row(EECC %>% ungroup() %>% select(-vertex_1) %>% rename(vertex = vertex_2)) %>% unique() %>%
  arrange(desc(size), cl_ID) %>%
  group_by(vertex) %>% summarise(n_3 = sum(size == 3), n_2 = sum(size == 2))

EECC_deg_prob <- EECC_deg_dist %>%
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

# function to evaluate the pgf at specific x and y values

cd_size_pgf <- function(x, y, n = 100, p1 = 0.09, alpha = 0, deg_pgf, dGdx, dGdy){
  p2 <- pk(2, p = p1, alpha = alpha)
  # define functions
  
  # offspring dist from triangle
  g_r <- function(x,y) dGdy(x,y)/dGdy(1,1)
  # offspring dist from edge
  g_q <- function(x,y) dGdx(x,y)/dGdx(1,1)
  
  # seed pgf
  H1 <- function(x, y, h1, h2, h5) (1-p1)^2 + 2*p1*(1-p1)*x*y*g_r(h5, h1)*h2 + (p1^2)*((x*y)^2)*((g_r(h5, h1))^2)
  H2 <- function(x, y, h1, h5) 1-p2 + p2*x*y*g_r(h5, h1)
  H5 <- function(x, y, h1, h5) 1 - p1 + p1*x*y*g_q(h5, h1)
  
  # ICs
  H_IC <- function(x) c(1,1,1)
  H_old <- H_IC(x)
  H_new <- numeric(3)
  
  for (i in 1:(n-1)) { # the third entry is actually Hn,5(x)
    H_new[1] <- H1(x*(y^(n-i)), y, H_old[1], H_old[2], H_old[3])
    H_new[2] <- H2(x*(y^(n-i)), y, H_old[1], H_old[3])
    H_new[3] <- H5(x*(y^(n-i)), y, H_old[1], H_old[3])
    H_old <- H_new
  }
  
  H_new[1] <- H1(x, y, H_old[1], H_old[2], H_old[3])
  H_new[2] <- H2(x, y, H_old[1], H_old[3])
  H_new[3] <- H5(x, y, H_old[1], H_old[3])
  
  full_pgf <- x*deg_pgf(H_new[3], H_new[1])
  return(full_pgf)
}

# for sanity, we might want to check that the pgf returns a 1 when x=y=1
cd_size_pgf(1,1, deg_pgf = EECC_pgf, dGdx = EECC_dGdx, dGdy = EECC_dGdy)

########################################################################
# Evaluate the pgf at points in the complex plane and save to csv file
########################################################################

pgf_cd_size_mv <- function(x,y) mapply(cd_size_pgf, x, y, MoreArgs =  c(n = 1000, p1 = 0.09, alpha = 0.0, deg_pgf = EECC_pgf, dGdx = EECC_dGdx, dGdy = EECC_dGdy))

M <- 50
N <- 50

x_vals <- exp(-2*pi*1i*(0:(M-1))/M)
y_vals <- exp(-2*pi*1i*(0:(N-1))/N)

pgf_evaluated.df <- expand_grid(x_vals, y_vals)

pgf_evaluated <- pgf_cd_size_mv(pgf_evaluated.df$x_vals, pgf_evaluated.df$y_vals)

pgf_to_invert <- pgf_evaluated.df %>% mutate(P_ft = pgf_evaluated) %>%
  mutate(x_re = Re(x_vals), x_im = Im(x_vals), y_re = Re(y_vals), y_im = Im(y_vals),
         P_re = Re(P_ft), P_im = Im(P_ft)) %>%
  select(-c(x_vals, y_vals, P_ft))

pgf_to_invert %>% write_csv("joint_dist_fft.csv")

# next step is to run the "joint_pgf_inversion.py" script

joint_dist_size_cd <- read_delim("net_science_probabilities.csv", col_names = FALSE, delim= ' ') %>% as.matrix()

# moments of the distribution to find the correlation

# X - size, Y - cumulative depth

E_X <- sum(rowSums(joint_dist_size_cd)*(0:49))
E_Y <- sum(colSums(joint_dist_size_cd)*(0:49))
E_X_squared <- sum(rowSums(joint_dist_size_cd)*((0:49)^2))
E_Y_squared <- sum(colSums(joint_dist_size_cd)*((0:49)^2))
E_XY <- sum(t(joint_dist_size_cd*(0:49))*(0:49))

pearson_correlation <- (E_XY - (E_X*E_Y))/(sqrt(E_X_squared - (E_X^2))*sqrt(E_Y_squared - (E_Y^2)))

EATD <- sum(t(joint_dist_size_cd[2:50,]/(1:49))*(0:49))

################################
# compare to simulations
################################

pois_net <- net_gen_double_poisson(n = 10000, v = 4, u = 1)
pois_adj <- as_adj(pois_net)

# make a cluster with 6 cores - adapt this to your own computer
cluster <- makeCluster(6)
registerDoParallel(cluster)

# increase "total" for more simulations
tic()
pois_sims <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(i, net = net, adj = adjacency, p1 = 0.09, alpha = 0.0, total = 166667)
toc()

stopImplicitCluster()

size_cd_dist.df <- pois_sims %>%
  group_by(ID) %>%
  summarise(size = n_distinct(c(parent, child)), cd = sum(generation)) %>%
  group_by(size, cd) %>% summarise(N = n()) %>% ungroup()

n_zero <- 6*166667 - sum(size_cd_dist.df$N)

size_cd_dist.df <- size_cd_dist.df %>% add_row(tibble(size = 1, cd = 0, N = n_zero)) %>% arrange(size, cd) %>% mutate(p = N/sum(N))

EATD_sim <- size_cd_dist.df %>% mutate(average_depth = cd/size) %>% summarise(EATD = sum(p*average_depth))

pearson_correlation_sim <- size_cd_dist.df %>% uncount(N) %>% select(-p) %>% summarise(pearson_correlation = cor(size, cd))
