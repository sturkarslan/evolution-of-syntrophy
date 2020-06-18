############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Process growth rate and yield data for various clonal isolate pairings
# and plot violin plots. Also calculates EOB.
###############################################################################

## Load clonal isolate data
setwd("~/Documents/GitHub/evolution-of-syntrophy/Figures/Figure 6")
clonal.data <- read.delim("UE3_ClonalIsolate_Pairings_Growth_Data.txt", header=T, sep="\t", stringsAsFactors=F)

# Excess over Bliss analysis
epds <- c("3", "9")
synergy.table <- data.frame()
for(epd in epds){
  ### Clone: D1, EPD:03
  # D1 evolved x Wt-Mmp
  f.03.D1 <- subset(clonal.data, EPD == epd & pairing_type == "Dv-E_Mm-A" & isolate_1_clone == 'D01')
  
  # Wt-Dv x M1 evolved
  f.03.M1 <- subset(clonal.data, EPD == epd & pairing_type == "Dv-A_Mm-E" & isolate_1_clone == 'M01')
  
  # D1 evolved x M1 evolved
  f.03.D1M1 <- subset(clonal.data, EPD == epd & pairing_type == "Dv-E_Mm-E" & isolate_1_clone == 'D01')
  
  # growth rates: D1 evolved x Wt-Mmp
  f.03.D1.mean.rate <- mean(f.03.D1$rate)
  f.03.D1.sd.rate <- sd(f.03.D1$rate)
  f.03.D1.mean.yield <- mean(f.03.D1$yield)
  f.03.D1.sd.yield <- sd(f.03.D1$yield)
  # growth rates: Wt-Dv x M1 evolved
  f.03.M1.mean.rate <- mean(f.03.M1$rate)
  f.03.M1.sd.rate <- sd(f.03.M1$rate)
  f.03.M1.mean.yield <- mean(f.03.M1$yield)
  f.03.M1.sd.yield <- sd(f.03.M1$yield)
  # growth rates: D1 evolved x M1 evolved
  f.03.D1M1.mean.rate <- mean(f.03.D1M1$rate)
  f.03.D1M1.sd.rate <- sd(f.03.D1M1$rate)
  f.03.D1M1.mean.yield <- mean(f.03.D1M1$yield)
  f.03.D1M1.sd.yield <- sd(f.03.D1M1$yield)
  # calculate fractions
  f.03.D1M1.rate.fraction <-  ((f.03.D1.mean.rate + f.03.M1.mean.rate) - (f.03.D1.mean.rate * f.03.M1.mean.rate))
  f.03.D1M1.rate.fraction.sd <-  (f.03.D1.sd.rate + f.03.M1.sd.rate) - (f.03.D1.sd.rate * f.03.M1.sd.rate)
  f.03.D1M1.yield.fraction <-  ((f.03.D1.mean.yield + f.03.M1.mean.yield) - (f.03.D1.mean.yield * f.03.M1.mean.yield))
  f.03.D1M1.yield.fraction.sd <-  (f.03.D1.sd.yield + f.03.M1.sd.yield) - (f.03.D1.sd.yield * f.03.M1.sd.yield)
  
  # calculate EOB
  EOB.03.D1.rate <- (f.03.D1M1.mean.rate - f.03.D1M1.rate.fraction) * 100
  EOB.03.D1.yield <- (f.03.D1M1.mean.yield - f.03.D1M1.yield.fraction) * 100
  
  synergy.table <- rbind.data.frame(synergy.table, cbind.data.frame(clone="D01", 
                                              epd = epd, 
                                              eob.rate=EOB.03.D1.rate, 
                                              eob.rate.sd =f.03.D1M1.rate.fraction.sd,
                                              eob.yield=EOB.03.D1.yield,
                                              eob.yield.sd=f.03.D1M1.yield.fraction.sd))
}


### EOB SYnergy Plots
synergy.melted <- melt(synergy.table)
synergy.melted$type <- sapply(synergy.melted$variable, function(i) strsplit(as.character(i), split=".", fixed = T)[[1]][2])
## function to have 2 digits for y-axis labels
scaleFUN <- function(x) sprintf("%.2f", x)
## EOB calculations for Growth Rate 
eob.rate.plot <- ggplot(data=synergy.table)
eob.rate.plot <- eob.rate.plot + geom_point(aes(x=epd, y=eob.rate), size=4, color="red")
eob.rate.plot <- eob.rate.plot + geom_errorbar(aes(x=epd, ymin=eob.rate-eob.rate.sd, ymax=eob.rate+eob.rate.sd), width=0.3, color="black")
eob.rate.plot <- eob.rate.plot + geom_hline(yintercept = 0, linetype="dashed")
eob.rate.plot <- eob.rate.plot + facet_grid(.~epd, scales = "free_x",space = "free_x", margins = c(5,5,5,5))
eob.rate.plot <- eob.rate.plot + theme(axis.title.y = element_text(), axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
eob.rate.plot <- eob.rate.plot + scale_y_continuous(labels=scaleFUN, limits = c(-0.1, 1.50), breaks = c(-0.5, 0, 0.5,1,1.50))
eob.rate.plot
## EOB Calculations for Growth Yield
eob.yield.plot <- ggplot(data=synergy.table)
eob.yield.plot <- eob.yield.plot + geom_point(aes(x=epd, y=eob.yield), size=4, color="blue", shape=17)
eob.yield.plot <- eob.yield.plot + geom_errorbar(aes(x=epd, ymin=eob.yield-eob.yield.sd, ymax=eob.yield+eob.yield.sd), width=0.3, color="black")
eob.yield.plot <- eob.yield.plot + geom_hline(yintercept = 0, linetype="dashed")
eob.yield.plot <- eob.yield.plot + facet_grid(.~epd, scales = "free_x", space = "free_x", margins = c(5,5,5,5))
eob.yield.plot <- eob.yield.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
eob.yield.plot <- eob.yield.plot + scale_y_continuous(labels=scaleFUN, limits = c(-8.50,0), breaks=c(-8,-6,-4,-2,0))
eob.yield.plot

grid.arrange(eob.rate.plot, eob.yield.plot)

### Growth Files
# load growth data fro grofit results
growth.file <- "growth_data_UA3_formatted_Nejc-v2.txt"
growth.data <- read.delim(growth.file, header = T, sep="\t", stringsAsFactors = F)
# melt growth data
growth.melted <- melt(growth.data, measure.vars = c("rate", "yield"), id.vars = colnames(growth.data)[3:13])

# ggplot boxpot for rate
growth.rate <- growth.melted[which(growth.melted$variable == "rate"),]
growth.yield <- growth.melted[which(growth.melted$variable == "yield"),]

growth.rate.f <- subset(growth.rate, (pairing_type %in% c("Dv-A_Mm-A","Dv-A_Mm-E","Dv-E_Mm-A","Dv-E_Mm-E","EPD") & clone %in% c("0","1","bulk")))
growth.yield.f <- subset(growth.yield, (pairing_type %in% c("Dv-A_Mm-A","Dv-A_Mm-E","Dv-E_Mm-A","Dv-E_Mm-E","EPD") & clone %in% c("0","1","bulk")))

rate.plot <- ggplot(growth.rate.f, aes(x=pairing_type, y=value, group=pairing_type))
rate.plot <- rate.plot + geom_violin()
rate.plot <- rate.plot + geom_point(aes(col=clone))
#rate.plot <- rate.plot + geom_point(aes(y=value, shape=isolate_1_clone))
rate.plot <- rate.plot + facet_grid(.~factor(EPD, levels=c("An", "3", "9")), scales = "free_x", space = "free_x")
rate.plot <- rate.plot + labs(x="pairings", y="Growth rate")
rate.plot <- rate.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
rate.plot

yield.plot <- ggplot(growth.yield.f, aes(x=pairing_type, y=value, group=pairing_type))
yield.plot <- yield.plot + geom_violin()
yield.plot <- yield.plot + geom_point(aes(col=clone))
#yield.plot <- yield.plot + geom_point(aes(y=value, shape=isolate_1_clone))
yield.plot <- yield.plot + facet_grid(.~factor(EPD, levels=c("An", "3", "9")), scales = "free_x", space = "free_x",)
yield.plot <- yield.plot + labs(x="pairings", y="Growth yield")
yield.plot <- yield.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
yield.plot

### Plot values for only both evolved for all clones

evolved <- subset(growth.melted, pairing_type %in% c("Dv-E_Mm-E", "EPD") )

# ggplot boxpot for rate
evolved.growth.rate <- evolved[which(evolved$variable == "rate"),]
evolved.growth.yield <- evolved[which(evolved$variable == "yield"),]

#evolved.growth.rate.f <- subset(evolved.growth.rate, (pairing_type %in% c("Dv-A_Mm-A","Dv-A_Mm-E","Dv-E_Mm-A","Dv-E_Mm-E","EPD") & clone %in% c("0","1","bulk")))
#evolved.growth.yield.f <- subset(evolved.growth.yield, (pairing_type %in% c("Dv-A_Mm-A","Dv-A_Mm-E","Dv-E_Mm-A","Dv-E_Mm-E","EPD") & clone %in% c("0","1","bulk")))

evolved.rate.plot <- ggplot(evolved.growth.rate, aes(x=clone, y=value, group=clone))
evolved.rate.plot <- evolved.rate.plot + geom_violin()
evolved.rate.plot <- evolved.rate.plot + geom_point(aes(col=clone))
#evolved.rate.plot <- evolved.rate.plot + geom_point(aes(y=value, shape=isolate_1_clone))
evolved.rate.plot <- evolved.rate.plot + facet_grid(.~factor(EPD, levels=c("An", "3", "9")), scales = "free_x", space = "free_x")
evolved.rate.plot <- evolved.rate.plot + labs(x="pairings", y="Growth rate")
evolved.rate.plot <- evolved.rate.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
#evolved.rate.plot <- evolved.rate.plot + scale_y_continuous(labels=function(x) sprintf("%.2f", x))
evolved.rate.plot

evolved.yield.plot <- ggplot(evolved.growth.yield, aes(x=clone, y=value, group=clone))
evolved.yield.plot <- evolved.yield.plot + geom_violin()
evolved.yield.plot <- evolved.yield.plot + geom_point(aes(col=clone))
#evolved.yield.plot <- evolved.yield.plot + geom_point(aes(y=value, shape=isolate_1_clone))
evolved.yield.plot <- evolved.yield.plot + facet_grid(.~factor(EPD, levels=c("An", "3", "9")), scales = "free_x", space = "free_x",)
evolved.yield.plot <- evolved.yield.plot + labs(x="pairings", y="Growth yield")
evolved.yield.plot <- evolved.yield.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
evolved.yield.plot <- evolved.yield.plot + scale_y_continuous(labels=function(x) sprintf("%.3f", x))

evolved.yield.plot



eg <- read.delim("~/Downloads/early_evolved_lines.csv", sep=",", header=T)
ue3 <- subset(eg, evolution_line %in% c("KSEA", "UA3"))
eg.rate.plot <- ggplot(data=ue3, aes(x=transfer, y=rate, group=transfer))
eg.rate.plot <- eg.rate.plot + geom_point()
eg.rate.plot <- eg.rate.plot + geom_boxplot()
eg.rate.plot

