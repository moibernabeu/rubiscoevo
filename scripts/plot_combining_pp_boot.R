# Moisès Bernabeu
# Started Cocentaina February 13th, 2020
# Bootstrap and posterior probabilities

setwd("~/Desktop/support/plot_ettema/")

library(treeio)
library(ggtree)
library(ggplot2)

tree <- read.newick("bootstrap.tree", node.label = "support")
posterior <- read.newick("posterior_probabilities.tre", node.label = "support")

tree@data$support_pp <- append(posterior@data$support, 
                               rep(NA, times = length(tree@data$support)-length(posterior@data$support)))

p0 <- ggtree(tree)
pnodes <- ggtree(tree) +
  geom_tiplab(size = 2, aes(label = p0$data$node)) +
  geom_nodelab(size = 2, aes(label = p0$data$node), nudge_x = 0.001) 
pnodes

tree <- groupOTU(tree, .node = c(1,2,3))
root <- rootnode(tree@phylo)

bp <- as.numeric(length(tree@data$support[tree@data$support <= 25]))
xl <- c(0, 0.2)

cut_pp <- cut(tree@data$support_pp, c(0.5, 0.75, 0.9, 0.9999999, 1)); cut_pp
cut_bp <- cut(tree@data$support, c(0.5, 0.75, 0.9, 0.9999999, 1)); cut_bp

combined <- c()
for (i in 1:length(tree@data$support)) {
  if (!is.na(cut_pp[i]) | !is.na(cut_bp[i])) { # Here there is a problem!!!!
    if (cut_pp[i] == '(0.9999999,1]' | cut_bp[i] == '(0.9999999,1]' | is.na(cut_bp[i]) | is.na(cut_pp[i])) {
      combined[i] <- 1
    } else if (cut_pp[i] == '(0.9,0.9999999]' | cut_bp[i] == '(0.9,0.9999999]' | is.na(cut_bp[i]) | is.na(cut_pp[i])) {
      combined[i] <- 2
    } else if (cut_pp[i] == '(0.75,0.9]' | cut_bp[i] == '(0.75,0.9]' | is.na(cut_bp[i]) | is.na(cut_pp[i])) {
      combined[i] <- 3
    } else if (cut_pp[i] == '(0.5,0.75]' | cut_bp[i] == '(0.5,0.75]' | is.na(cut_bp[i]) | is.na(cut_pp[i])) {
      combined[i] <- 4
    } 
  } else {combined[i] <- NA}
}
combined

tree@data$combined <- factor(combined, levels = c(1,2,3),
                                labels = c("100/100", "0.9-1/0.9-1","0.75-0.9/0.75-0.9"))
tree@data$combined

p <- ggtree(tree, aes(linetype = group)) +
  geom_tiplab(hjust = -0.060, align = TRUE, geom = "text", size = 2.5) +
  xlim(xl) + 
  geom_point2(aes(x = branch, subset=!isTip & node != root,
                  fill=combined), 
              shape=21) +
  theme_tree(legend.position="bottom", bgcolor = NA) +
  scale_fill_manual(values=c("black", "darkgrey", "grey"),
                    name='BP/PP',
                    breaks = c("100/100", "0.9-1/0.9-1","0.75-0.9/0.75-0.9"),
                    labels=expression("100/100", "0.9-1/0.9-1","0.75-0.9/0.75-0.9")) +
  scale_linetype(guide = "none")
p
