library(cowplot)
library(ggpubr)
# install.packages(c("factoextra", "FactoMineR"))

# adjust font family and size in annotate to match that of ggplot 
# https://stackoverflow.com/questions/65057107/annotate-font-is-different-from-the-other-texts-in-the-figure-with-ggplot2
family <- "sans"
fontsize <- 16 / .pt

gg.pentagon <- function(data) {
  colour.crosswalk <- c("positive" = "black", "negative" = "red")
  alpha.crosswalk <- c("non-significant" = 0.1, "marginal-significant" = 0.3, "significant"=0.95)
  
  nudge <- 0.25
  ggplot(data = data, aes(x = x, y = y)) +
    annotate(geom = "text", x = 0, y = (1+nudge), label = "Invariability", size = fontsize, family = "sans") +
    annotate(geom = "text", x = (-0.951-nudge), y = 0.309, angle = 72, label = "Resistance_Dry", size = fontsize, family = "sans", color="#E69F00") +
    annotate(geom = "text", x = (0.951+nudge), y = .309, angle = -72, label = "Resistance_Wet", size = fontsize, family = "sans", color= "#0072B2") +
    annotate(geom = "text", x = -0.6, y = (-0.809-nudge), angle = -38, label = "Recovery_Dry", size = fontsize, family = "sans", color="#E69F00") +
    annotate(geom = "text", x = 0.6, y = (-0.809-nudge), angle = 38, label = "Recovery_Wet", size = fontsize, family = "sans", color= "#0072B2") +
    geom_segment(aes(xend = xend, yend = yend, 
                     colour = as.character(direction), alpha=as.character(sig), linewidth = width),
                 lineend = "round") +
    scale_colour_manual(values = colour.crosswalk) +
    scale_alpha_manual(values = alpha.crosswalk) +
    lims(x = c((-0.951-2*nudge), (0.951+2*nudge)),
         y = c((-0.809-3*nudge), (1+2*nudge))) +
    theme_void() +
    theme(legend.position = "none",
          aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", colour = "white"))
}


gg.triangle <- function(data) {
  colour.crosswalk <- c("positive" = "black", "negative" = "red")
  alpha.crosswalk <- c("non-significant" = 0.1, "marginal-significant" = 0.3, "significant"=0.95)
  
  nudge <- 0.25
  ggplot(data = data, aes(x = x, y = y)) +
    annotate(geom = "text", x = 0, y = (1+nudge), label = "Biomass", size = fontsize, family = "sans") +
    annotate(geom = "text", x = -0.866, y = (-0.5-nudge), angle = -38, label =  "Richness", size = fontsize, family = "sans") +
    annotate(geom = "text", x = 0.866, y = (-0.5-nudge), angle = 38, label = "Composition", size = fontsize, family = "sans") +
    geom_segment(aes(xend = xend, yend = yend, 
                     colour = as.character(direction), alpha=as.character(sig), linewidth = width),
                 lineend = "round") +
    scale_colour_manual(values = colour.crosswalk) +
    scale_alpha_manual(values = alpha.crosswalk) +
    lims(x = c((-0.866-2.2*nudge), (0.866+2.4*nudge)),
         y = c((-0.5-3*nudge), (1+2*nudge))) +
    theme_void() +
    theme(legend.position = "none",
          aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", colour = "white"))
}

