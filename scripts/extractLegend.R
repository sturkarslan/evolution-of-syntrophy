############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Helper function  to extract legend from ggplot object
###############################################################################

library('ggplot2');library('reshape2');library(gridExtra);library(gplots)

extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}