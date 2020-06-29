# Moisès Bernabeu
# Started Cocentaina February 4th, 2020
# Density plot

library(ggplot2)
library(ggtree)
library(phylosignal)
library(phytools)
library(phylotools)
library(ape)
library(tibble)
library(tidyr)
library(Biostrings)
library(treeio)
library(grid)
library(cowplot)
library(dplyr)
library(stringr)
require(reshape2)

setwd("../../real_seqs/selected_trees/model_trees/")

#----Density trees----
#trees import
trees <- read.newick("nt_sim_model_trees_sorted.nwk")

#rename
# taxa <- data.frame(old = read.csv("../../../preliminar/names/raw_seqs_tree_names.txt", sep = "\t")[,1],
#                       new = read.csv("../../../preliminar/names/raw_seqs_tree_names.txt", sep = "\t")[,3])
# for (i in 1:length(trees)) trees[[i]] <- rename_taxa(trees[[i]], taxa, key = old, value = new)

#group trees
trees_p <- list(trees[[1]] %>% fortify %>% mutate(tree="a"))
for (i in 2:length(trees)) trees_p[[i]] <- trees[[i]] %>% fortify %>% mutate(tree="b")
trees_p[[length(trees_p)+1]] <- trees_p[[1]]
trees_p[[1]] <- NULL

#plot density with selected model in black
den <- ggdensitree(trees_p, aes(colour = tree), align.tips = TRUE) +
  # geom_tiplab(colour='black', align = TRUE, size = 1.5) +
  # xlim(0,1) +
  scale_color_manual(values = alpha(c("black", "steelblue"), c(1, .06))) +
  theme(legend.position = "none")
den

