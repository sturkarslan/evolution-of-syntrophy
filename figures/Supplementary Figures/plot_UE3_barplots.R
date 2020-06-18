############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Plot barplots for UE3 line of Dv and Mm
# For barplots in Figure 5 and Supplementary Figures 2-3
###############################################################################
library('ggplot2'); library('reshape2');library(gridExtra);library(gplots);library('pheatmap')
source("~/Documents/Github/evolution-of-syntrophy/scripts/extractLegend.R")

# load early generation mutations data
mut.data.dv <- read.delim("~/Documents/GitHub/evolution-of-syntrophy/data/dvh_mutations_allsamples_attributes_3112020.txt", header=T, sep="\t", stringsAsFactors=FALSE)
mut.data.mm <- read.delim("~/Documents/GitHub/evolution-of-syntrophy/data/mmp_mutations_allsamples_attributes_3112020.txt", header=T, sep="\t", stringsAsFactors=FALSE)
ci.names <- read.delim("~/Documents/GitHub/evolution-of-syntrophy/data/clonal_isolate_pairs_IDs.txt", header=T, sep="\t")

# get DV Clones
ci.names.dv <- ci.names[grep("D", ci.names$isolate),]
mut.data.dv$sample.name <- mut.data.dv$sample
mut.data.dv <- mut.data.dv[grep("CI_00", mut.data.dv$sample, invert = T),]

# get Mm Clones
ci.names.mm <- ci.names[grep("M", ci.names$isolate),]
mut.data.mmsample.name <- mut.data.mm$sample
mut.data.mm <- mut.data.mm[grep("(CI_00|CI_36|CI_37)", mut.data.mm$sample, invert = T),]
mut.data.mm <- mut.data.mm[grep("WT", mut.data.mm$sample, invert = T),]

### replace geneid for intergenic mutations with IG for DV and translate clonal isolate names
for(row in 1:length(mut.data.dv$variant_id)){
  variant <- mut.data.dv[row,"variant_id"]
  locus_name <- strsplit(as.character(variant), split = "-", fixed=T)[[1]][2]
  position <- strsplit(as.character(variant), split = "-", fixed=T)[[1]][3]
  #copy sample names into a new column
  sample <- mut.data.dv[row,"sample"]
  sample.type <- strsplit(as.character(sample), split = "_", fixed=T)[[1]][1]
  sample.no <- strsplit(as.character(sample), split = "_", fixed=T)[[1]][2]
  
  if(sample.type == "CI"){
    sample.name <-  paste(as.character(ci.names.dv[which(ci.names.dv$pair == sample.no),"isolate"]), sep = "")
  } else {
    sample.name <- sample
  }
  
  if(locus_name == "IG"){
    locus_name <- paste("IG", as.character(position), sep="_")
  }
  mut.data.dv[row,"locus_name"] <- locus_name
  mut.data.dv[row,"sample.name"] <- sample.name
}

### replace geneid for intergenic mutations with IG for Mm and translate clonal isolate names
for(row in 1:length(mut.data.mm$variant_id)){
  variant <- mut.data.mm[row,"variant_id"]
  locus_name <- strsplit(as.character(variant), split = "-", fixed=T)[[1]][2]
  position <- strsplit(as.character(variant), split = "-", fixed=T)[[1]][3]
  #copy sample names into a new column
  sample <- mut.data.mm[row,"sample"]
  sample.type <- strsplit(as.character(sample), split = "_", fixed=T)[[1]][1]
  sample.no <- strsplit(as.character(sample), split = "_", fixed=T)[[1]][2]
  
  if(sample.type == "CI" & !(sample.no %in% c("36","37"))){
    sample.name <-  paste(as.character(ci.names.mm[which(ci.names.mm$pair == sample.no),"isolate"]), sep = "")
  } else {
    sample.name <- sample
  }
  
  if(locus_name == "IG"){
    locus_name <- paste("IG", as.character(position), sep="_")
  }
  mut.data.mm[row,"locus_name"] <- locus_name
  mut.data.mm[row,"sample.name"] <- sample.name
}

###### format early generations mutations data for dvh
mut.data.eg.dv <- mut.data.dv[grep(c("(AN_Dv-Ancestor-1|_UE3|CI(_2_|_6_|_8_|_20_|_21_|_1(1_|3_|4_)|_37_))"), mut.data.dv$sample),]  
mut.data.eg.dv$transfer <- sapply(mut.data.eg.dv$sample, function(i) strsplit(as.character(i), split = "-")[[1]][2])
mut.data.eg.dv$line <- sapply(mut.data.eg.dv$sample, function(i) strsplit(as.character(i), split = "-")[[1]][1])
mut.data.eg.dv$freq2 <- as.numeric(as.character(sub("%", "",mut.data.eg.dv$freq)))
mut.data.eg.dv[is.na(mut.data.eg.dv$transfer),"transfer"] <- 152
## select UA3 line (change it if you want to filter for the line)
mut.data.ua3.dv <- mut.data.eg.dv[grep("", mut.data.eg.dv$sample),]
# order based on transfer
mut.data.ua3.dv <- mut.data.ua3.dv[order(as.numeric(as.character(mut.data.ua3.dv$transfer))),]
# add line info
mut.data.ua3.dv$eline <- sapply(mut.data.ua3.dv$line, function(i) strsplit(as.character(i), split = "_", fixed = T)[[1]][2])

# get ancestral mutations
ancestor.mutations <- mut.data.dv[which(mut.data.dv$sample == "AN_Dv-Ancestor-1"),"variant_id"]
for(variant in unique(mut.data.ua3.dv$variant_id)){
  if(variant %in% ancestor.mutations){
    mut.data.ua3.dv[which(mut.data.ua3.dv$variant_id == variant),"ancestor"] <- "ancestral"
  } else {
    mut.data.ua3.dv[which(mut.data.ua3.dv$variant_id == variant),"ancestor"] <- "non-ancestral"
  }
}

###### format early generations mutations data for mmp
mut.data.eg.mm <- mut.data.mm[grep(c("(AN_Coculture-Ancestor|_UE3|CI(_3_|_7_|_20_|_27_|_29_|_3(1_|2_|4_)|_10_))"), mut.data.mm$sample),]  
mut.data.eg.mm$transfer <- sapply(mut.data.eg.mm$sample, function(i) strsplit(as.character(i), split = "-")[[1]][2])
mut.data.eg.mm$line <- sapply(mut.data.eg.mm$sample, function(i) strsplit(as.character(i), split = "-")[[1]][1])
mut.data.eg.mm$freq2 <- as.numeric(as.character(sub("%", "",mut.data.eg.mm$freq)))
mut.data.eg.mm[is.na(mut.data.eg.mm$transfer),"transfer"] <- 152
## select UA3 line (change it if you want to filter for the line)
mut.data.ua3.mm <- mut.data.eg.mm[grep("", mut.data.eg.mm$sample),]
# order based on transfer
mut.data.ua3.mm <- mut.data.ua3.mm[order(as.numeric(as.character(mut.data.ua3.mm$transfer))),]
# add line info
mut.data.ua3.mm$eline <- sapply(mut.data.ua3.mm$line, function(i) strsplit(as.character(i), split = "_", fixed = T)[[1]][2])

# get ancestral mutations
ancestor.mutations.mm <- mut.data.mm[which(mut.data.mm$sample == "AN_Coculture-Ancestor"),"variant_id"]
for(variant in unique(mut.data.ua3.mm$variant_id)){
  if(variant %in% ancestor.mutations.mm){
    mut.data.ua3.mm[which(mut.data.ua3.mm$variant_id == variant),"ancestor"] <- "ancestral"
  } else {
    mut.data.ua3.mm[which(mut.data.ua3.mm$variant_id == variant),"ancestor"] <- "non-ancestral"
  }
}

#slect non-ancestral mutations only dvh
#mutations.1000.dv <- mut.data.ua3.dv[which(mut.data.ua3.dv$ancestor != "ancestral"),]
mutations.1000.dv <- mut.data.ua3.dv
#slect non-ancestral mutations only Mm
#mutations.1000.mm <- mut.data.ua3.mm[which(mut.data.ua3.mm$ancestor != "ancestral"),]
## if you keep ancestral mutations for plotting
mutations.1000.mm <- mut.data.ua3.mm


##### barplot for dv
my.samples.dv <-
  subset(
    mutations.1000.dv,
    sample.name %in% c(
      "AN_Dv-Ancestor-1",
      "TG_UE3",
      "EP_UE3_03",
      "UA3.152.03.D01",
      "UA3.152.03.D02",
      "UA3.152.03.D03",
      "EP_UE3_10",
      # "UA3.152.10.D01",
      "UA3.152.10.D02",
      "UA3.152.10.D03",
      "EP_UE3_09",
      "UA3.152.09.D01",
      "UA3.152.09.D02",
      "UA3.152.09.D03"
    )
  )

my.samples.dv$variant_id <- sub("Chr-", "", my.samples.dv$variant_id)
my.samples.dv$variant_id <- sub("pDV-", "", my.samples.dv$variant_id)

sample.order <- c(
  "AN_Dv-Ancestor-1",
  "TG_UE3",
  "EP_UE3_03",
  "UA3.152.03.D01",
  "UA3.152.03.D02",
  "UA3.152.03.D03",
  "EP_UE3_10",
  #"UA3.152.10.D01",
  "UA3.152.10.D02",
  "UA3.152.10.D03",
  "EP_UE3_09",
  "UA3.152.09.D01",
  "UA3.152.09.D02",
  "UA3.152.09.D03"
)


library(dplyr)
variant.order <- my.samples.dv %>%
  group_by(variant_id) %>%
  tally(sort = T) %>%
  select(variant_id) %>%
  pull()

dv.bar <- ggplot(my.samples.dv, aes(x=factor(variant_id, levels=variant.order), y=freq2, fill=variant_id))
dv.bar <- dv.bar + geom_bar(stat="identity")
dv.bar <- dv.bar + coord_flip()
dv.bar <- dv.bar + facet_grid(.~factor(sample.name, levels=sample.order))
dv.bar <- dv.bar + theme(legend.position = "none")
dv.bar <- dv.bar + labs(y="Frequency", x="Mutations")
dv.bar

ggsave(dv.bar, filename = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/UE3_dvh_frequency_barplot.pdf",device = "pdf")

##### barplot for mm
my.samples.mm <-
  subset(
    mutations.1000.mm,
    sample.name %in% c(
      "AN_Coculture-Ancestor",
      "TG_UE3",
      "EP_UE3_03",
      "UA3.152.03.M01",
      "UA3.152.03.M02",
      "UA3.152.03.M03",
      "EP_UE3_10",
      "UA3.152.10.M01",
      "UA3.152.10.M02",
      "UA3.152.10.M03",
      "EP_UE3_09",
      "UA3.152.09.M01",
      "UA3.152.09.M02",
      "UA3.152.09.M03"
    )
  )

my.samples.mm$variant_id <- sub("BX950229-", "", my.samples.mm$variant_id)

sample.order <- c(
  "AN_Coculture-Ancestor",
  "TG_UE3",
  "EP_UE3_03",
  "UA3.152.03.M01",
  "UA3.152.03.M02",
  "UA3.152.03.M03",
  "EP_UE3_10",
  "UA3.152.10.M01",
  "UA3.152.10.M02",
  "UA3.152.10.M03",
  "EP_UE3_09",
  "UA3.152.09.M01",
  "UA3.152.09.M02",
  "UA3.152.09.M03"
)


library(dplyr)
variant.order <- my.samples.mm %>%
  group_by(variant_id) %>%
  tally(sort = T) %>%
  select(variant_id) %>%
  pull()

mm.bar <- ggplot(my.samples.mm, aes(x=factor(variant_id, levels=variant.order), y=freq2, fill=variant_id))
mm.bar <- mm.bar + geom_bar(stat="identity", width = 0.5)
mm.bar <- mm.bar + coord_flip()
mm.bar <- mm.bar + facet_grid(.~factor(sample.name, levels=sample.order))
mm.bar <- mm.bar + theme(legend.position = "none")
mm.bar <- mm.bar + labs(y="Frequency", x="Mutations")
mm.bar

ggsave(mm.bar, filename = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/UE3_mmp_frequency_barplot.pdf",device = "pdf")


