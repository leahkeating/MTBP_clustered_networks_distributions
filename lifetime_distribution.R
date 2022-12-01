########################################################################
# Title: lifetime_distribution.R
# Purpose: calculate the distribution of cascade lifetimes
# Author: Leah Keating
# Date last edited: 26 September 2022
########################################################################

rm(list = ls()) # tidy work space
gc()

source("functions.R")

double_poisson_pgf <- function(x, y, v = 2, u = 2) exp(u*(x-1) + v*(y - 1))
pois_dGds1 <- function(x, y, v = 2, u = 2) v*double_poisson_pgf(x, y, v, u) # partial derivative wrt s1
pois_dGds3 <- function(x, y, v = 2, u = 2) u*double_poisson_pgf(x, y, v, u) # partial derivative wrt s3

# the pgf for the number of particles in generation N:

particle_pgf <- function(z,p,alpha,N=1000, deg_pgf, dGdx, dGdy){
  if(N == 0){
    R_tilde <- z
  }else{
    g_r <- function(x,y) dGdy(x,y)/dGdy(1,1)
    # offspring dist from edge
    g_q <- function(x,y) dGdx(x,y)/dGdx(1,1)
    p1 <- pk(1,p,alpha)
    p2 <- pk(2,p,alpha)
    R <- c((1-p1)^2 + 2*p1*(1-p1)*z + (p1^2)*z^2, 1-p2+p2*z, 1-p1 + p1*z)
    R_new <- numeric(3)
    if(N>1){
      for (i in 2:N) {
        R_new[1] <- (1-p1)^2 + 2*p1*(1-p1)*R[2]*g_r(R[3],R[1]) + (p1^2)*(g_r(R[3],R[1])^2)
        R_new[2] <- 1 - p2 + p2*g_r(R[3],R[1])
        R_new[3] <- 1 - p1 + p1*g_q(R[3],R[1])
        R <- R_new
      }
    }
    R_tilde <- deg_pgf(R[3],R[1])
  }
  return(R_tilde)
}

particle_pgf(z = 0, p = 0.2, N = 7, alpha = 0.01, deg_pgf = double_poisson_pgf, dGdx = pois_dGds3, dGdy = pois_dGds1)

# the probability that there are no particles in generation in can be used to find the
# probability that the lifetime is N generations

prob_no_particles_in_gen_n <- function(N_max) sapply(N_max, particle_pgf, z = 0, p = 0.05, alpha = 0.2, deg_pgf = double_poisson_pgf, dGdx = pois_dGds3, dGdy = pois_dGds1)

lifetime_dist_theory.df <- prob_no_particles_in_gen_n(N_max = 0:20) %>%
  tibble(prob_no_particles = .) %>%
  mutate(n = 0:20) %>%
  # lag gives the "prob_no_particles" for the previous  generation
  mutate(omega_n = prob_no_particles - lag(prob_no_particles)) %>%
  filter(n>0) %>%
  rename(lifetime = n) %>%
  mutate(cdf = cumsum(omega_n), ccdf = 1-cdf)

lifetime_dist_theory.df %>%
  ggplot(aes(x = lifetime, y = omega_n)) +
  geom_line() +
  scale_y_log10()

#########################
# Compare to simulations
#########################

net <- net_gen_double_poisson(n = 10000, v = 2, u = 2)
adjacency <- as_adj(net)

# Simulate cascades

cluster <- makeCluster(6)
registerDoParallel(cluster)

tic()
cascades.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(j = i, net = net, adj = adjacency, p = 0.05, alpha = 0.2, total = 16667)
toc()

stopImplicitCluster()

lifetime_dist_sims.df <- cascades.df %>% group_by(ID) %>% summarise(lifetime = (max(generation)+1)) %>%
  group_by(lifetime) %>% summarise(n = n())
num_lifetime_1 <- lifetime_dist_sims.df %>% ungroup() %>% summarise(total = (16667*6-sum(n))) %>% as.numeric()
lifetime_dist_sims.df <- lifetime_dist_sims.df %>% add_row(lifetime = 1, n = num_lifetime_1) %>% arrange(lifetime)
lifetime_dist_sims.df <- lifetime_dist_sims.df %>% mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)

require(scales)
lifetime_dist_sims.df %>% rename(omega_n = p) %>%
  ggplot(aes(x = lifetime, y = omega_n)) +
  geom_line(data = lifetime_dist_theory.df, colour = "black", size = 0.75) +
  geom_point(colour = "black", fill = "#f56b02", size = 2, alpha = 1, 
             shape = 21) +
  scale_y_log10(limits = c(10^(-7),1), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("probability") +
  xlim(c(0,20)) +
  xlab("n") +
  ylab(TeX("$\\Omega_{n}$"))
