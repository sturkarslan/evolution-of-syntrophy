library('ggplot2');library('reshape2');library(gridExtra);library(gplots)

# Extract legend from ggplot object
extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}