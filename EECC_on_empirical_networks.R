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
V(net)$name <- as.character(1:gorder(net))
adjacency <- as_adj(net)

EECC <- EECC_generator(net, max_cl_size = 3)

EECC %>% group_by(cl_ID) %>%
  summarise(clique_size = n_distinct(c(vertex_1, vertex_2))) %>% arrange(desc(clique_size)) %>%
  group_by(clique_size) %>% summarise(n())
