setwd("OneDrive - Universitat de Valencia/RECERCA/PROJECTS/TFG/rubevo/cured/convergence")

library(ape)

x <- read.tree("01_aa_nobp.nwk")
y <- read.tree("04_aat_sim_nobp.nwk")

comparePhylo(x, y, plot = T, force.rooted = T,
             use.edge.length = F)
