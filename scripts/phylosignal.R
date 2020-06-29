# Moisès Bernabeu
# Started Cocentaina April 4th, 2020
# Phylogenetic signal

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

# trees import
aa_trees <- read.newick("aa/no_trimmed/model_trees.nwk")
for (i in 1:length(aa_trees)) {
  aa_null_branches <- which(aa_trees[[i]]$edge.length==0)
  aa_trees[[1]]$edge.length[aa_null_branches]<-0.00000001
}
for (i in 1:123) aa_trees[[i]]$node.label <- NULL

aat_trees <- read.newick("aa/trimmed/model_trees.nwk")
for (i in 1:length(aat_trees)) {
  aat_null_branches <- which(aat_trees[[i]]$edge.length==0)
  aat_trees[[1]]$edge.length[aat_null_branches]<-0.00000001
}
for (i in 1:123) aat_trees[[i]]$node.label <- NULL

nt_trees <- read.newick("nt/no_trimmed/model_trees.nwk")
for (i in 1:length(nt_trees)) {
  nt_null_branches <- which(nt_trees[[i]]$edge.length==0)
  nt_trees[[1]]$edge.length[nt_null_branches]<-0.00000001
}
for (i in 1:88) nt_trees[[i]]$node.label <- NULL

ntt_trees <- read.newick("nt/trimmed/model_trees.nwk")
for (i in 1:length(ntt_trees)) {
  ntt_null_branches <- which(ntt_trees[[i]]$edge.length==0)
  ntt_trees[[1]]$edge.length[ntt_null_branches]<-0.00000001
}
for (i in 1:88) ntt_trees[[i]]$node.label <- NULL

# phylogenetic signal table
# amino acid
aa_sig <- data.frame("Cmean" = 1, "pval" = 1, "I" = 1,
                   "pval" = 1, "Lambda" = 1, "pval" = 1,
                   "K" = 1, "pval" = 1, "K.star" = 1, "pval" = 1)
for (i in 1:123) {
  t <- phylo4d(aa_trees[[i]]) 
  # random
  tipData(t)$rand <- rnorm(length(aa_trees))
  # signal
  tipData(t)$BM <- rTraitCont(as(t, "phylo"), model = "BM")
  a <- phyloSignal(t)
  aa_sig <- rbind(aa_sig, c(a$stat$Cmean[2], a$pvalue$Cmean[2], a$stat$I[2], a$pvalue$I[2],
                            a$stat$Lambda[2], a$pvalue$Lambda[2], a$stat$K[2], a$pvalue$K[2],
                            a$stat$K.star[2], a$pvalue$K.star[2]))
}
write.table(file = "aa/no_trimmed/aa_phylogenetic_signals.csv", x = aa_sig, sep = ";", dec = ".")

aat_sig <- data.frame("Cmean" = 1, "pval" = 1, "I" = 1,
                     "pval" = 1, "Lambda" = 1, "pval" = 1,
                     "K" = 1, "pval" = 1, "K.star" = 1, "pval" = 1)
for (i in 1:123) {
  t <- phylo4d(aat_trees[[i]]) 
  # random
  tipData(t)$rand <- rnorm(length(aat_trees))
  # signal
  tipData(t)$BM <- rTraitCont(as(t, "phylo"), model = "BM")
  a <- phyloSignal(t)
  aat_sig <- rbind(aat_sig, c(a$stat$Cmean[2], a$pvalue$Cmean[2], a$stat$I[2], a$pvalue$I[2],
                              a$stat$Lambda[2], a$pvalue$Lambda[2], a$stat$K[2], a$pvalue$K[2],
                              a$stat$K.star[2], a$pvalue$K.star[2]))
}
write.table(file = "aa/trimmed/aat_phylogenetic_signals.csv", x = aat_sig, sep = ";", dec = ".")

# nucleotides
nt_sig <- data.frame("Cmean" = 1, "pval" = 1, "I" = 1,
                     "pval" = 1, "Lambda" = 1, "pval" = 1,
                     "K" = 1, "pval" = 1, "K.star" = 1, "pval" = 1)
for (i in 1:123) {
  t <- phylo4d(nt_trees[[i]]) 
  # random
  tipData(t)$rand <- rnorm(length(nt_trees[[i]]))
  # signal
  tipData(t)$BM <- rTraitCont(as(t, "phylo"), model = "BM")
  a <- phyloSignal(t)
  nt_sig <- rbind(nt_sig, c(a$stat$Cmean[2], a$pvalue$Cmean[2], a$stat$I[2], a$pvalue$I[2],
                            a$stat$Lambda[2], a$pvalue$Lambda[2], a$stat$K[2], a$pvalue$K[2],
                            a$stat$K.star[2], a$pvalue$K.star[2]))
}
write.table(file = "nt/no_trimmed/nt_phylogenetic_signals.csv", x = nt_sig, sep = ";", dec = ".")

ntt_sig <- data.frame("Cmean" = 1, "pval" = 1, "I" = 1,
                     "pval" = 1, "Lambda" = 1, "pval" = 1,
                     "K" = 1, "pval" = 1, "K.star" = 1, "pval" = 1)
for (i in 1:123) {
  t <- phylo4d(ntt_trees[[i]]) 
  # random
  tipData(t)$rand <- rnorm(length(ntt_trees[[i]]))
  # signal
  tipData(t)$BM <- rTraitCont(as(t, "phylo"), model = "BM")
  a <- phyloSignal(t)
  ntt_sig <- rbind(ntt_sig, c(a$stat$Cmean[2], a$pvalue$Cmean[2], a$stat$I[2], a$pvalue$I[2],
                              a$stat$Lambda[2], a$pvalue$Lambda[2], a$stat$K[2], a$pvalue$K[2],
                              a$stat$K.star[2], a$pvalue$K.star[2]))
}
write.table(file = "nt/trimmed/ntt_phylogenetic_signals.csv", x = ntt_sig, sep = ";", dec = ".")
