library('ggplot2');library('reshape2');library('pheatmap');library('grid');library('gridExtra')
setwd('/Volumes/Macintosh HD/Users/serdarturkaslan/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates')
# load growth data fro grofit results
growth.file <- "growth_data_UA3_formatted_Nejc.txt"
growth.data <- read.delim(growth.file, header = T, sep="\t", stringsAsFactors = F)
# melt growth data
growth.melted <- melt(growth.data, measure.vars = c("rate", "yield"), id.vars = colnames(growth.data)[3:13])

# ggplot boxpot for rate
growth.rate <- growth.melted[which(growth.melted$variable == "rate"),]
growth.yield <- growth.melted[which(growth.melted$variable == "yield"),]

growth.rate.f <- subset(growth.rate, (pairing_type %in% c("Dv-A_Mm-A","Dv-A_Mm-E","Dv-E_Mm-A","Dv-E_Mm-E","EPD") & clone %in% c("0","3","bulk")))
growth.yield.f <- subset(growth.yield, (pairing_type %in% c("Dv-A_Mm-A","Dv-A_Mm-E","Dv-E_Mm-A","Dv-E_Mm-E","EPD") & clone %in% c("0","3","bulk")))

rate.plot <- ggplot(growth.rate.f, aes(x=pairing_type, y=value, group=pairing_type))
rate.plot <- rate.plot + geom_violin()
rate.plot <- rate.plot + geom_point(aes(col=clone))
#rate.plot <- rate.plot + geom_point(aes(y=value, shape=isolate_1_clone))
rate.plot <- rate.plot + facet_grid(.~factor(EPD, levels=c("An", "3", "9")), scales = "free_x")
rate.plot <- rate.plot + labs(x="pairings", y="Growth rate", title = "Growth Rate Comparison", subtitle = "UE3 Clonal Isolate Pairings, EPDs and Ancestors")
rate.plot <- rate.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
rate.plot

yield.plot <- ggplot(growth.yield.f, aes(x=pairing_type, y=value, group=pairing_type))
yield.plot <- yield.plot + geom_violin()
yield.plot <- yield.plot + geom_point(aes(col=clone))
#yield.plot <- yield.plot + geom_point(aes(y=value, shape=isolate_1_clone))
yield.plot <- yield.plot + facet_grid(.~factor(EPD, levels=c("An", "3", "9")), scales = "free_x")
yield.plot <- yield.plot + labs(x="pairings", y="Growth yield", title = "Growth yield Comparison", subtitle = "UE3 Clonal Isolate Pairings, EPDs and Ancestors")
yield.plot <- yield.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
yield.plot

grid.arrange(rate.plot, yield.plot)

do.ttest <- function(input_data = growth.rate.f, epd = NULL){
  col.names <- unique(input_data$pairing_type) 
  row.names <- unique(input_data$pairing_type) 
  ttest.matrix <- matrix(nrow = length(col.names), ncol = length(col.names), dimnames = list(col.names,col.names))
  
  for(name1 in col.names){
    for(name2 in col.names){
      if(name1 == "Dv-A_Mm-A"){
        set1 <- subset(input_data, subset = (pairing_type == name1 & EPD == "An"))
      }else{
        set1 <- subset(input_data, subset = (pairing_type == name1 & EPD == epd))
      }
      if(name2 == "Dv-A_Mm-A"){
        set2 <- subset(input_data, subset = (pairing_type == name2 & EPD == "An"))
      }else{
        set2 <- subset(input_data, subset = (pairing_type == name2 & EPD == epd))
      }
      if(name1 == name2){
        pvalue <- NA
      } else{
        pvalue <- t.test(set1$value, set2$value)$p.value
      }
      ttest.matrix[name1,name2] <- pvalue
    }
  }
  return(ttest.matrix)
}

epd3.growth_rate.matrix <- do.ttest(input_data = growth.rate.f, epd = 3)
epd3.growth_rate.matrix[upper.tri(epd3.growth_rate.matrix)] <- NA

epd9.growth_rate.matrix <- do.ttest(input_data = growth.rate.f, epd = 9)
epd9.growth_rate.matrix[upper.tri(epd9.growth_rate.matrix)] <- NA

epd3.growth_yield.matrix <- do.ttest(input_data = growth.yield.f, epd = 3)
epd3.growth_yield.matrix[upper.tri(epd3.growth_yield.matrix)] <- NA

epd9.growth_yield.matrix <- do.ttest(input_data = growth.yield.f, epd = 9)
epd9.growth_yield.matrix[upper.tri(epd9.growth_yield.matrix)] <- NA


heatmap.rate3 <- pheatmap(log(epd3.growth_rate.matrix), scale = "none", cluster_rows = F, cluster_cols = F, display_numbers = T)
heatmap.rate9 <- pheatmap(log(epd9.growth_rate.matrix), scale = "none", cluster_rows = F, cluster_cols = F, display_numbers = T)
heatmap.yield3 <- pheatmap(log(epd3.growth_yield.matrix), scale = "none", cluster_rows = F, cluster_cols = F, display_numbers = T)
heatmap.yield9 <- pheatmap(log(epd9.growth_yield.matrix), scale = "none", cluster_rows = F, cluster_cols = F, display_numbers = T)

grid.arrange(heatmap.rate3, heatmap.rate9, heatmap.yield3, heatmap.yield9)
par(mfrow=c(2,2))







