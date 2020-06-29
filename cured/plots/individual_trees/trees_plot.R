# Moisès Bernabeu
# Started February 4th, 2020
# Tree plot

library(phytools)
library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggrepel)

# Files import
tree <- read.newick("trees/sim_seqs/ntt_sim.tree", node.label = "support")
tree@phylo <- root(tree@phylo, outgroup = "758802_AQT80891_rbcL")
tree

# Taxonomy data
taxa <- data.frame(old = read.csv("../../../preliminar/names/tree_species_names.txt", sep = "\t")[,1],
                   new = read.csv("../../../preliminar/names/tree_species_names.txt", sep = "\t")[,3])
head(taxa)
tnamed <- rename_taxa(tree, taxa, key = old, value = new)
tnamed

# Plot outgroups with different linetype
## Select nodes
# p0 <- ggtree(tree)
# pnodes <- ggtree(tnamed) +
#   geom_tiplab(size = 2, aes(label = p0$data$node)) +
#   geom_nodelab(size = 2, aes(label = p0$data$node), nudge_x = 0.001) 
# pnodes

# ggtree(tnamed) +
#   geom_label_repel(aes(label=node), size = 1)

## Define outgroups and its nodes
# tnamed <- groupOTU(tnamed, .node = c(46,47,48))
root <- rootnode(tnamed@phylo)

to_drop <- c("Nocardia cyriacigeorgica",
             "Nocardia mangyaensis",
             "Nocardia otitidiscaviarum",
             "Nocardia seriolae",
             "Nocardia terpenica",
             "Mycolicibacterium litorale")
tnamed <- drop.tip(tnamed, to_drop)

# Plot tree
tnamed@data$support[which(tnamed@data$support <= 75)] <- NA
tnamed@data$support <- cut(x = tnamed@data$support, c(100, 99.9999999, 90, 85, 75),
                         labels = c("100", "90-99", "85-90", "75-85"),
                         ordered_result = F)

p <- ggtree(tnamed, size = 0.2) +
  geom_tiplab(hjust = -0.0060, align = F, geom = "text", size = 1.25) +
  # xlim(0, 0.9) +
  geom_point2(aes(x = branch, subset = !isTip & node != root,
                  fill = support), shape = 21,lwd = 0.9, stroke = 0) +
  theme_tree(legend.position = "bottom", bgcolor = NA) +
  scale_fill_grey(na.translate = FALSE,
                  name = "BP") +
  scale_linetype(guide = "none")
p

# Branch transformation for better visualisation
p$data[p$data$node %in% c(81), "x"] <- 1.2*mean(p$data$x)
p

pnodes <- ggtree(tnamed) +
  geom_tiplab(size = 2, aes(label = p$data$node)) +
  geom_nodelab(size = 2, aes(label = p$data$node), nudge_x = 0.001)
pnodes

p2 <- rotate(p, 127) %>% rotate(128) %>% rotate(168)
p2
