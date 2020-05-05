## Load clonal isolate data
setwd("~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/")
### Growth Files
# load growth data fro grofit results
growth.file <- "growth_data_UE3_formatted_Kris.txt"
growth.data <- read.delim(growth.file, header = T, sep="\t", stringsAsFactors = F)
# melt growth data
growth.melted <- melt(growth.data, measure.vars = c("rate_Kris", "OD_Kris" ), id.vars = c("pair","pairing_type","inter_vs_intra","dvh","D.line","D.generation","D.epd","D.colony","mmp","M.line","M.generation","M.epd","M.colony","rate_Kris_sd","OD_Kris_sd"))

# ggplot boxpot for rate
growth.rate <- growth.melted[which(growth.melted$variable == "rate_Kris"),]
growth.yield <- growth.melted[which(growth.melted$variable == "OD_Kris"),]

growth.rate.f <- subset(growth.rate, (pairing_type %in% c("Dan_Mev","Dev_Man","Dev_Mev")))
growth.yield.f <- subset(growth.yield, (pairing_type %in% c("Dan_Mev","Dev_Man","Dev_Mev")))

rate.plot <- ggplot(growth.rate.f, aes(x=pairing_type, y=value, group=pairing_type))
rate.plot <- rate.plot + geom_violin()
rate.plot <- rate.plot + geom_point(aes(col=D.colony))
#rate.plot <- rate.plot + geom_point(aes(y=value, shape=isolate_1_clone))
rate.plot <- rate.plot + facet_grid(.~factor(inter_vs_intra), scales = "free_x", space = "free_x")
rate.plot <- rate.plot + labs(x="pairings", y="Growth rate")
rate.plot <- rate.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
rate.plot

yield.plot <- ggplot(growth.yield.f, aes(x=pairing_type, y=value, group=pairing_type))
yield.plot <- yield.plot + geom_violin()
yield.plot <- yield.plot + geom_point(aes(col=dvh_strain))
#yield.plot <- yield.plot + geom_point(aes(y=value, shape=isolate_1_clone))
yield.plot <- yield.plot + facet_grid(.~factor(line), scales = "free_x", space = "free_x",)
yield.plot <- yield.plot + labs(x="pairings", y="Growth yield")
yield.plot <- yield.plot + theme(axis.text.x = element_text(angle = 90, size = 13), strip.text.x = element_text(size = 14, colour = "red"))
yield.plot
