############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Plot heatmap for HR2 line of Dv and Mm
# also creates SiFit files for lineage building
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
mut.data.eg.dv <- mut.data.dv[grep(c("(_HR2|CI(_27_|_28_|_19_|_36_|_16_|_5_|_22_|_33_|_24_))"), mut.data.dv$sample),]  
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
ancestor.mutations <- mut.data.dv[which(mut.data.dv$sample == "AN_Coculture-Ancestor"),"variant_id"]
for(variant in unique(mut.data.ua3.dv$variant_id)){
  if(variant %in% ancestor.mutations){
    mut.data.ua3.dv[which(mut.data.ua3.dv$variant_id == variant),"ancestor"] <- "ancestral"
  } else {
    mut.data.ua3.dv[which(mut.data.ua3.dv$variant_id == variant),"ancestor"] <- "non-ancestral"
  }
}

###### format early generations mutations data for mmp
mut.data.eg.mm <- mut.data.mm[grep(c("(_HR2|CI(_19_|_14_|_35_|_17_|_11_|_24_|_6_|_8_|_16_))"), mut.data.mm$sample),]  
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
mutations.1000.dv <- mut.data.ua3.dv[which(mut.data.ua3.dv$ancestor != "ancestral"),]
#slect non-ancestral mutations only Mm
mutations.1000.mm <- mut.data.ua3.mm[which(mut.data.ua3.mm$ancestor != "ancestral"),]


plot.line.pheatmap <-function(org=c("dvh","mmp")){
  library(RColorBrewer)
  library(viridis)
  
  dv.order <- c("TG_HR2",
                "EP_HR2_01","HR2.152.01.D01","HR2.152.01.D02","HR2.152.01.D03",
                "EP_HR2_05","HR2.152.05.D01","HR2.152.05.D02","HR2.152.05.D03",
                "EP_HR2_10","HR2.152.10.D01","HR2.152.10.D02","HR2.152.10.D03"
  )
  
  mm.order <- c("TG_HR2",
                "EP_HR2_01","HR2.152.01.M01","HR2.152.01.M02","HR2.152.01.M03",
                "EP_HR2_05","HR2.152.05.M01","HR2.152.05.M02","HR2.152.05.M03",
                "EP_HR2_10","HR2.152.10.M01","HR2.152.10.M02","HR2.152.10.M03"
  )
  
  ## Select samples
  if(org == "dvh"){
    sample.order <- dv.order
    mutations.file <- mutations.1000.dv
    mutations.file$variant_id <- sub("Chr-", "", mutations.file$variant_id)
    mutations.file$variant_id <- sub("pDV-", "", mutations.file$variant_id)
  }
  if(org == "mmp"){
    sample.order <- mm.order
    mutations.file <- mutations.1000.mm
    mutations.file$variant_id <- sub("BX950229-", "", mutations.file$variant_id)
  }
  
  ## select samples for freq heatmap
  my.samples.freq <- subset(mutations.file, sample.name %in% sample.order)
  
  
  ## create matrix for freq
  my.matrix.freq <- matrix(
    nrow = length(unique(my.samples.freq$variant_id)), 
    ncol = length(unique(my.samples.freq$sample.name)),
    dimnames = list(unique(my.samples.freq$variant_id), 
                    unique(my.samples.freq$sample.name))
  )
  
  ## fill the matrix
  for(variant in unique(my.samples.freq$variant_id)){
    for(sample in my.samples.freq$sample.name){
      my.data <- my.samples.freq[which(my.samples.freq$variant_id == variant & my.samples.freq$sample.name == sample),"freq2"]
      if(length(my.data) != 0){
        my.matrix.freq[variant,sample] <- my.data
      }else{
        my.matrix.freq[variant,sample] <- 0
      }
    }
  }
  
  
  ## custom sort matrix
  my.matrix.freq <- my.matrix.freq[,order(factor(colnames(my.matrix.freq), levels=sample.order))]
  
  ## annotations for the sample type
  #annotations <- data.frame(Type = c("Early-Gen","Early-Gen","Early-Gen","Early-Gen","1K","EPD","Clonal","Clonal","Clonal","EPD","Clonal","Clonal","Clonal","EPD","Clonal","Clonal","Clonal"))
  col.annot <- unique(data.frame(my.samples.freq[,c("sample.name","experiment")]))
  row.names(col.annot) <- col.annot$sample.name
  col.annot <- col.annot[-1]
  ## plot heatmap
  pheatmap(clustering_distance_rows = "correlation",my.matrix.freq,cluster_cols = F,cluster_rows = T, annotation_legend = F, treeheight_row = 0, scale = "none", color = viridis(30),border_color = "#333333", annotation_col = col.annot, drop_levels = F)
}

## multiple plotting
mm <- list(plot.line.pheatmap(org = "dvh")[[4]])
mm[[2]] <- plot.line.pheatmap(org = "mmp")[[4]]
z <- do.call(grid.arrange,mm)
plot(z)


#### SiFit files


#### dv matrix for SiFit
sifit_files_dv <- function() {
  my.samples <- subset(
    mutations.1000.dv,
    sample.name %in% c(
      "TG_HR2",
      "EP_HR2_05",
      "HR2.152.05.D01",
      "HR2.152.05.D02",
      "HR2.152.05.D03",
      "EP_HR2_10",
      "HR2.152.10.D01",
      "HR2.152.10.D02",
      "HR2.152.10.D03",
      "EP_HR2_01",
      "HR2.152.01.D01",
      "HR2.152.01.D02",
      "HR2.152.01.D03"
    )
  )
  
  # create matrix
  my.matrix <- matrix(
    nrow = length(unique(my.samples$variant_id)),
    ncol = length(unique(my.samples$sample.name)),
    dimnames = list(
      unique(my.samples$variant_id),
      unique(my.samples$sample.name)
    )
  )
  
  # fill the matrix
  for (variant in unique(my.samples$variant_id)) {
    for (sample in my.samples$sample.name) {
      my.data <-
        my.samples[which(my.samples$variant_id == variant &
                           my.samples$sample.name == sample), "freq2"]
      if (length(my.data) != 0) {
        my.matrix[variant, sample] <- 1
      } else{
        my.matrix[variant, sample] <- 0
      }
    }
  }
  
  # write matrix files
  #mutation names
  write.table(
    row.names(my.matrix),
    file = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/siFit/hr2_dv_mutation_names.txt",
    sep = " ",
    row.names = F,
    col.names = F,
    quote = F
  )
  # sample names
  writeLines(colnames(my.matrix), con = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/siFit/hr2_dv_sample_names.txt", sep =
               " ")
  # mutation matrix
  my.matrix <- cbind(id = seq(1, dim(my.matrix)[1]), my.matrix)
  write.table(
    my.matrix,
    file = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/siFit/hr2_dv_mutation_matrix.txt",
    sep = " ",
    row.names = F,
    col.names = F,
    quote = F
  )
}


#### mm matrix for SiFit
sifit_files_mm <- function() {
  my.samples.mm <- subset(
    mutations.1000.mm,
    sample.name %in% c(
      "TG_HR2",
      "EP_HR2_05",
      "HR2.152.05.M01",
      "HR2.152.05.M02",
      "HR2.152.05.M03",
      "EP_HR2_10",
      "HR2.152.10.M01",
      "HR2.152.10.M02",
      "HR2.152.10.M03",
      "EP_HR2_01",
      "HR2.152.01.M01",
      "HR2.152.01.M02",
      "HR2.152.01.M03"
    )
  )
  
  # create matrix
  my.matrix.mm <- matrix(
    nrow = length(unique(my.samples.mm$variant_id)),
    ncol = length(unique(my.samples.mm$sample.name)),
    dimnames = list(
      unique(my.samples.mm$variant_id),
      unique(my.samples.mm$sample.name)
    )
  )
  
  # fill the matrix
  for (variant in unique(my.samples.mm$variant_id)) {
    for (sample in my.samples.mm$sample.name) {
      my.data <-
        my.samples.mm[which(my.samples.mm$variant_id == variant &
                              my.samples.mm$sample.name == sample), "freq2"]
      if (length(my.data) != 0) {
        my.matrix.mm[variant, sample] <- 1
      } else{
        my.matrix.mm[variant, sample] <- 0
      }
    }
  }
  
  # write matrix files
  #mutation names
  write.table(
    row.names(my.matrix.mm),
    file = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/siFit/hr2_mm_mutation_names.txt",
    sep = " ",
    row.names = F,
    col.names = F,
    quote = F
  )
  # sample names
  writeLines(colnames(my.matrix.mm), con = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/siFit/hr2_mm_sample_names.txt", sep =
               " ")
  # mutation matrix
  my.matrix.mm <- cbind(id = seq(1, dim(my.matrix.mm)[1]), my.matrix.mm)
  write.table(
    my.matrix.mm,
    file = "~/Documents/GitHub/evolution-of-syntrophy/Clonal-Isolates/siFit/hr2_mm_mutation_matrix.txt",
    sep = " ",
    row.names = F,
    col.names = F,
    quote = F
  )
}



