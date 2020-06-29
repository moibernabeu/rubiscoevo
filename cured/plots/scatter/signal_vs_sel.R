# Moisès Bernabeu
# Started Cocentaina February 4th, 2020
# Model selection criterion vs. phylogenetic signal

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(readxl)

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

#----Plot main theme----
# Data import
aa_data <- read.csv("aa_models_signal.csv", header = TRUE, sep = ";", row.names = 1)
nt_data <- read.csv("nt_models_signal.csv", header = TRUE, sep = ";", row.names = 1)

aa_data <- read_excel("data.xlsx", sheet = 6)
nt_data <- read_excel("data.xlsx", sheet = 8)

aa_lnl <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(lnL, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(lnL, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(lnL, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(lnL, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(lnL, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("-lnL")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

aa_BIC <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(BIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(BIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(BIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(BIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(BIC, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("BIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

aa_AIC <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(AIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AIC, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("AIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

aa_AICc <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(AICc, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AICc, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AICc, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AICc, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AICc, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("AICc")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

ntt_sim_lnl <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(lnL, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(lnL, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(lnL, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(lnL, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(lnL, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("-lnL")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

nt_BIC <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(BIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(BIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(BIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(BIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(BIC, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("BIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

nt_AIC <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(AIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AIC, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("AIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

nt_AICc <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(AICc, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AICc, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AICc, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AICc, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AICc, I, colour = "I")) +
  main_theme +
  theme(legend.position = "bottom") +
  xlab("AICc")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "darkgrey", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "darkgrey", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

arr <- ggarrange(plotlist = list(aa_BIC, nt_BIC, aa_lnl, nt_lnl, aa_AIC, nt_AIC),
          common.legend = TRUE, ncol = 2, nrow = 3,
          legend = "bottom")
annotate_figure(arr,
                left = text_grob("Phylogenetic signal", color = "darkgrey",
                                 rot = 90, vjust = 1.7, just = "center",
                                 hjust = 0.25))

#-----------Classic-----------
aa_lnl <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(lnL, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(lnL, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(lnL, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(lnL, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(lnL, I, colour = "I")) +
  theme_classic() +
  xlab("-lnL")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

aa_BIC <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(BIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(BIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(BIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(BIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(BIC, I, colour = "I")) +
  theme_classic() +
  xlab("BIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

aa_AIC <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(AIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AIC, I, colour = "I")) +
  theme_classic() +
  xlab("AIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

aa_AICc <- ggplot(aa_data) +
  geom_point(alpha=0.4, aes(AICc, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AICc, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AICc, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AICc, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AICc, I, colour = "I")) +
  theme_classic() +
  xlab("AICc")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

nt_lnl <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(lnL, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(lnL, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(lnL, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(lnL, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(lnL, I, colour = "I")) +
  theme_classic() +
  xlab("-lnL")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

nt_BIC <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(BIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(BIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(BIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(BIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(BIC, I, colour = "I")) +
  theme_classic() +
  xlab("BIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

nt_AIC <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(AIC, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AIC, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AIC, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AIC, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AIC, I, colour = "I")) +
  theme_classic() +
  xlab("AIC")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

nt_AICc <- ggplot(nt_data) +
  geom_point(alpha=0.4, aes(AICc, Cmean, colour = "Cmean")) +
  geom_point(alpha=0.4, aes(AICc, Lambda, colour = "Pagel's Lambda")) +
  geom_point(alpha=0.4, aes(AICc, K, colour = "Bloomberg's K")) +
  geom_point(alpha=0.4, aes(AICc, Kstar, colour = "K star")) +
  geom_point(alpha=0.4, aes(AICc, I, colour = "I")) +
  theme_classic() +
  xlab("AICc")+
  ylab("") +
  labs(colour = "Signal method") +
  theme(legend.position = "bottom", legend.justification = "centre") +
  annotation_ticks(color = "black", sides = c("l"), ticklength = -0.5 * unit(0.1, "cm")) +
  annotation_ticks(color = "black", sides = c("b"), ticklength = -0.5 * unit(0.1, "cm")) +
  coord_cartesian(clip = "off")

arr <- ggarrange(plotlist = list(aa_BIC, nt_BIC, aa_lnl, nt_lnl, aa_AIC, nt_AIC),
                 common.legend = TRUE, ncol = 2, nrow = 3,
                 legend = "bottom")
arr
arr2 <- ggarrange(plotlist = list(aa_lnl, aa_sim_lnl, aat_lnl, aat_sim_lnl,
                                  nt_lnl, nt_sim_lnl, ntt_lnl, ntt_sim_lnl),
                 common.legend = TRUE, ncol = 4, nrow = 2, align = "hv",
                 legend = "bottom")
arr2
annotate_figure(arr2,
                left = text_grob("Phylogenetic signal", color = "darkgrey",
                                 rot = 90, vjust = 1.7, just = "center",
                                 hjust = 0.25))
