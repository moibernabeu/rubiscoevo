library(ape)
library(phytools)
library(phylotools)
library(phylosignal)
library(ade4)
library(adephylo)
library(phylobase)
library(treeio)
library(magrittr)
library(dplyr)
library(stringi)

treename <- "aat_sim"
phy <- read.newick("simulated/aat_sim_model_trees_sorted.nwk")
for (i in 1:length(phy)) {
  phy_null <- which(phy[[i]]$edge.length==0)
  phy[[1]]$edge.length[phy_null]<-0.00000001
}
for (i in 1:length(phy)) phy[[i]]$node.label <- NULL

physig <- data.frame("Cmean" = 1, "pval" = 1, "I" = 1,
                      "pval" = 1, "Lambda" = 1, "pval" = 1,
                      "K" = 1, "pval" = 1, "K.star" = 1, "pval" = 1)

for (i in 1:length(phy)) {
  t <- phylo4d(phy[[i]])
  tipLabels(t) <- stri_rand_strings(length(t@label), 3, pattern = "[A-Za-z]")
  tipData(t)$rand <- rnorm(length(t@label))
  tipData(t)$BM <- rTraitCont(as(t, "phylo"), model = "BM")
  a <- phyloSignal(t)
  physig <- rbind(physig, c(a$stat$Cmean[2], a$pvalue$Cmean[2], a$stat$I[2], a$pvalue$I[2],
                              a$stat$Lambda[2], a$pvalue$Lambda[2], a$stat$K[2], a$pvalue$K[2],
                              a$stat$K.star[2], a$pvalue$K.star[2]))
}; physig <- physig[-1,]; head(physig)

write.csv(physig, file = paste(treename, "_phylosig.csv", sep = ""))
