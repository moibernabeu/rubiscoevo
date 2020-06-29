setwd("OneDrive - Universitat de Valencia/RECERCA/PROJECTS/TFG/rubevo/cured/all_data_csv/")

library(ape)
library(mclust)

#----Main theme----
main_theme <- theme(panel.border = element_rect(colour = "darkgrey", fill = NA, size = 0.5),
                    panel.background = element_blank(),
                    panel.grid = element_blank(),
                    plot.background = element_blank(),
                    legend.direction = "horizontal",
                    legend.justification = c("right", "top"),
                    legend.position = c(.99, .99),
                    legend.background = element_blank(),
                    legend.text = element_text(colour = "black"),
                    legend.title = element_text(colour = "black"),
                    legend.key = element_rect(fill = NA, color = NA),
                    axis.ticks = element_line(colour = "darkgrey"),
                    axis.text = element_text(colour = "darkgrey"),
                    axis.title = element_text(colour = "darkgrey"),
                    plot.title = element_text(colour = "black"),
                    plot.subtitle = element_text(colour = "black"),
                    plot.caption = element_text(colour = "black"),
                    plot.tag = element_text(colour = "black", face = "bold"))

#----Plot minor ticks functions----
GeomTicks <- ggproto(
  "GeomTicks", Geom,
  extra_params = "",
  handle_na = function(data, params) {
    data
  },
  
  draw_panel = function(data,
                        panel_scales,
                        coord,
                        base = c(10, 10),
                        sides = c("b", "l"),
                        scaled = TRUE,
                        ticklength = unit(0.1, "cm"),
                        ticks_per_base = base - 1,
                        delog = c(x = TRUE, y = TRUE)) {
    ticks <- list()
    
    for (s in 1:length(sides)) {
      if (grepl("[b|t]", sides[s])) {
        
        xticks <- panel_scales$x.minor
        
        # Make the grobs
        if (grepl("b", sides[s])) {
          ticks$x_b <- with(
            data,
            segmentsGrob(
              x0 = unit(xticks, "npc"),
              x1 = unit(xticks, "npc"),
              y0 = unit(0, "npc"),
              y1 = ticklength,
              gp = gpar(
                col = alpha(colour, alpha),
                lty = linetype,
                lwd = size * .pt
              )
            )
          )
        }
        if (grepl("t", sides[s])) {
          ticks$x_t <- with(
            data,
            segmentsGrob(
              x0 = unit(xticks, "npc"),
              x1 = unit(xticks, "npc"),
              y0 = unit(1, "npc"),
              y1 = unit(1, "npc") - ticklength,
              gp = gpar(
                col = alpha(colour, alpha),
                lty = linetype,
                lwd = size * .pt
              )
            )
          )
        }
      }
      
      
      if (grepl("[l|r]", sides[s])) {
        
        yticks <- panel_scales$y.minor
        
        # Make the grobs
        if (grepl("l", sides[s])) {
          ticks$y_l <- with(
            data,
            segmentsGrob(
              y0 = unit(yticks, "npc"),
              y1 = unit(yticks, "npc"),
              x0 = unit(0, "npc"),
              x1 = ticklength,
              gp = gpar(
                col = alpha(colour, alpha),
                lty = linetype, lwd = size * .pt
              )
            )
          )
        }
        if (grepl("r", sides[s])) {
          ticks$y_r <- with(
            data,
            segmentsGrob(
              y0 = unit(yticks, "npc"),
              y1 = unit(yticks, "npc"),
              x0 = unit(1, "npc"),
              x1 = unit(1, "npc") - ticklength,
              gp = gpar(
                col = alpha(colour, alpha),
                lty = linetype,
                lwd = size * .pt
              )
            )
          )
        }
      }
    }
    gTree(children = do.call("gList", ticks))
  },
  default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = 1)
)

annotation_ticks <- function(sides = "b",
                             scale = "identity",
                             scaled = TRUE,
                             ticklength = unit(0.1, "cm"),
                             colour = "black",
                             size = 0.5,
                             linetype = 1,
                             alpha = 1,
                             color = NULL,
                             ticks_per_base = NULL,
                             ...) {
  if (!is.null(color)) {
    colour <- color
  }
  
  # check for invalid side
  if (grepl("[^btlr]", sides)) {
    stop(gsub("[btlr]", "", sides), " is not a valid side: b,t,l,r are valid")
  }
  
  # split sides to character vector
  sides <- strsplit(sides, "")[[1]]
  
  if (length(sides) != length(scale)) {
    if (length(scale) == 1) {
      scale <- rep(scale, length(sides))
    } else {
      stop("Number of scales does not match the number of sides")
    }
  }
  
  base <- sapply(scale, function(x) switch(x, "identity" = 10, "log10" = 10, "log" = exp(1)),
                 USE.NAMES = FALSE)
  
  if (missing(ticks_per_base)) {
    ticks_per_base <- base - 1
  } else {
    if ((length(sides) != length(ticks_per_base))) {
      if (length(ticks_per_base) == 1) {
        ticks_per_base <- rep(ticks_per_base, length(sides))
      } else {
        stop("Number of ticks_per_base does not match the number of sides")
      }
    }
  }
  
  delog <- scale %in% "identity"
  
  layer(
    data = data.frame(x = NA),
    mapping = NULL,
    stat = StatIdentity,
    geom = GeomTicks,
    position = PositionIdentity,
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(
      base = base,
      sides = sides,
      scaled = scaled,
      ticklength = ticklength,
      colour = colour,
      size = size,
      linetype = linetype,
      alpha = alpha,
      ticks_per_base = ticks_per_base,
      delog = delog,
      ...
    )
  )
}


# Data --------------------------------------------------------------------
files <- list.files(pattern = "csv")

pdf("trees_hierarchy.pdf")
for (i in 1:length(files)) {
  # Arbre
  data <- read.csv(files[i], row.names = 2)
  subdata <- data[-c(1,3:5,7,9,11,13,15)]
  d <- dist(subdata, method = "euclidean")
  plot(hclust(d, method="ward.D"), type = "unrooted",
       cex = 0.4, no.margin = TRUE)
}; dev.off()

pdf("trees.pdf")
for (i in 1:length(files)) {
  # Arbre
  data <- read.csv(files[i], row.names = 2)
  subdata <- data[-c(1,3:5,7,9,11,13,15)]
  d <- dist(subdata, method = "euclidean")
  plot(as.phylo(hclust(d, method="ward.D")), type = "unrooted",
       cex = 0.4, no.margin = TRUE)
}; dev.off()

pdf("cluster.pdf")
for (i in 1:length(files)) {
  # Cluster
  data <- read.csv(files[i], row.names = 2)
  subdata <- data[-c(1,3:5,7,9,11,13,15)]
  name <- paste("clust_", i, sep = "")
  assign(name, Mclust(subdata))
  plot(get(name), what = "classification")
  a <- summary(get(name))
  capture.output(a, file = "clust_analysis.txt", append = T)
}; dev.off()

library(factoextra)
# BIC values used for choosing the number of clusters
fviz_mclust(clust_1, "BIC", palette = "jco")
# Classification: plot showing the clustering

for (i in 1:length(files)) {
  assign(paste("clust_final_plot_", i, sep = ""), fviz_mclust(get(paste("clust_", i, sep = "")), "classification", geom = "point", 
                                                             pointsize = 1.5, palette = "jco") +
           main_theme +
           theme(axis.line =  element_line(size = 0, colour = "grey80")) +
           labs(title = "", subtitle = ""))
  get(paste("clust_final_plot", i, sep = ""))
}

pdf("classif_plots.pdf", width = 6.5, height = 4)
clust_final_plot_8
dev.off()

# Classification uncertainty
fviz_mclust(clust_1, "uncertainty", palette = "jco")

# Individual --------------------------------------------------------------
data <- read.csv("SM_10.csv", row.names = 2)
subdata <- data[-c(1,3:5,7,9,11,13,15)]

# Arbre
d <- dist(subdata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D")
plot(fit, cex=0.3) # display dendogram
plot(as.phylo(fit), type = "unrooted", cex = 0.4, no.margin = TRUE)

# Cluster
fit <- Mclust(subdata)
plot(fit) # plot results
summary(fit)

plot <- recordPlot(clust_1); plot
