############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Create plots for visualizing Early Generation mutations for Dv and Mm
###############################################################################

library('ggplot2');library('reshape2');library(gridExtra);library(gplots);library(pheatmap);library('tictoc')
setwd("/Volumes/Macintosh HD/Users/serdarturkaslan/Documents/GitHub/evolution-of-syntrophy/Early-Gen")
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
source("../scripts/extractLegend.R")
source("../scripts/load_mutation_data.R")

### Function to load mutation Data and modify gene names ###
get.modify.data <- function(mutation.file=NULL){
  ## load mutation data
  mutation.data <- read.delim(mutation.file, header=T, sep="\t")
  
  ## replace geneid for intergenic mutations with IG
  cat("Processing the data file and modifying gene IDs...\n")
  pb <- txtProgressBar(min = 1, max = length(mutation.data$variant_id), style = 3)
  for(row in 1:length(mutation.data$variant_id)){
    variant <- mutation.data[row,"variant_id"]
    locus_name <- strsplit(as.character(variant), split = "-", fixed=T)[[1]][2]
    position <- strsplit(as.character(variant), split = "-", fixed=T)[[1]][3]
    callers = strsplit(as.character(mutation.data[row,"predictor"]), split = ":")[[1]]
    frequencies = as.numeric(sub("%","",strsplit(as.character(mutation.data[row,"freq_predictor"]),split = ":")[[1]]))
    #cat("callers:", callers, "frequencies:", frequencies,"\n")
    #if(length(mutation.data[row,"freq2"]) == 0){
    if("samtools" %in% callers){
      #cat("samtools\n")
      final_freq = frequencies[grep("samtools",callers)]
    } else if ("varscan" %in% callers & length(callers) > 1){
      final_freq = frequencies[grep("varscan",callers)]
      #cat("varscan\n")
    } else {
      #cat("nathing\n")
      final_freq = 0
    }
    mutation.data[row,"final_freq"] <- final_freq
    if(locus_name == "IG"){
      locus_name <- paste("IG", as.character(position), sep="_")
    }
    mutation.data[row,"locus_name"] <- locus_name
    setTxtProgressBar(pb, row)
  }
  close(pb)
  
  # Only include early geenration, 1K and Ancestor data
  mutation.data.filtered <- mutation.data[mutation.data$sample != "EG_HA3_118.New" & (mutation.data$experiment == "1000-gen" | mutation.data$experiment == "Early-gen" | mutation.data$sample == "AN_Coculture-Ancestor"),]
  
  # Format transfer and line names
  mutation.data.filtered$line <- sapply(mutation.data.filtered$sample, function(i) strsplit(as.character(i), split = "_")[[1]][2])
  mutation.data.filtered$transfer <- sapply(mutation.data.filtered$sample, function(i) strsplit(as.character(i), split = "_")[[1]][3])
  mutation.data.filtered$transfer <- sub("118.New", 119, mutation.data.filtered$transfer)
  mutation.data.filtered$transfer <- sub("118.Early", 118, mutation.data.filtered$transfer)
  mutation.data.filtered[which(mutation.data.filtered$line == "Coculture-Ancestor"),"transfer"] <- 0
  mutation.data.filtered[which(mutation.data.filtered$experiment == "1000-gen"),"transfer"] <- 152
  # remove 1000-gen lines without early gen data
  mutation.data.filtered <- mutation.data.filtered[!(mutation.data.filtered$line %in% c("HE2", "HR1", "UA2", "UE2")),]
  # order based on transfer
  mutation.data.filtered <- mutation.data.filtered[order(as.numeric(as.character(mutation.data.filtered$transfer))),]
  # get ancestral mutations and add ancestral annotation to mutation data frame
  ancestor.mutations <- mutation.data[which(mutation.data$sample == "AN_Coculture-Ancestor"),"variant_id"]
  for(variant in unique(mutation.data.filtered$variant_id)){
    if(variant %in% ancestor.mutations){
      mutation.data.filtered[which(mutation.data.filtered$variant_id == variant),"ancestor"] <- "ancestral"
    } else {
      mutation.data.filtered[which(mutation.data.filtered$variant_id == variant),"ancestor"] <- "non-ancestral"
    }
  }
  
  #select non-ancestral mutations only
  mutation.data.noancestor <- mutation.data.filtered[which(mutation.data.filtered$ancestor != "ancestral"),]
  # some statistics
  unique.ancestral.variants <- unique(mutation.data.filtered[which(mutation.data.filtered$ancestor == "ancestral"),"variant_id"])
  unique.ancestral.locus <-  unique(mutation.data.filtered[which(mutation.data.filtered$ancestor == "ancestral"),"locus_name"])
  unique.nonancestral.variants <- unique(mutation.data.filtered[which(mutation.data.filtered$ancestor == "non-ancestral"),"variant_id"])
  unique.nonancestral.locus <- unique(mutation.data.filtered[which(mutation.data.filtered$ancestor == "non-ancestral"),"locus_name"])
  
  cat("Statistics\n" )
  cat("Ancestral variants: ", length(unique.ancestral.variants), "Ancestral locus: ", length(unique.ancestral.locus), "\n" )
  cat("Early-Gen variants: ", length(unique.nonancestral.variants), "Early-Gen locus: ", length(unique.nonancestral.locus), "\n" )
  
  return(list(mutation.data=mutation.data, ancestral.variants=unique.ancestral.variants, mutation.data.noancestor=mutation.data.noancestor))
}

## Function to create an ordering based on # of mutations ##
order.by.mutations <- function(input.data=mutation.data.df$mutation.data.noancestor){
  order.table <- data.frame()
  for(locus in unique(input.data$locus_name)){
    count <- length(input.data[which(input.data$locus_name == locus),"variant_id"])
    order.table <- rbind(order.table, cbind(locus, count))
    order.table <- order.table[order(as.numeric(as.character(order.table$count)),decreasing = T),]
  }
  ordered.locus.counts <- as.vector(order.table$locus)
  return(ordered.locus.counts)
}

## create ordering from GScores
  order.by.gscore <- function(gscores.data=NULL,
                            input.data= mutation.data.df$mutation.data.noancestor ){
  gscores <- read.delim(gscores.data, sep="\t")
  gscores$new.locus <- sapply(gscores$locus, function(i) sub("_","", i))
  my.locuses <- unique(input.data$locus_name)
  gscores.selected = data.frame()
  # get intergenic and other locus that were not in the gscore list but in the early generation list
  missing.locus <- setdiff(my.locuses, gscores$new.locus)
  for(locus in my.locuses){
    gscores.selected <- rbind(gscores.selected, gscores[which(gscores$new.locus == locus),])
    locus.order <- gscores.selected$new.locus
  }
  ordered.locus.gscore <- append(locus.order, missing.locus)
  gscores.locus.only <- gscores[which(gscores$observed > 1),"new.locus"]
  
  return(list(ordered.locus.gscore=ordered.locus.gscore, gscores.locus.only=gscores.locus.only))
}

## order by ranking across lines ##
order.by.rank <- function(input.data= mutation.data.df$mutation.data.noancestor){
  mutation.data <- input.data
  locus.list <- unique(mutation.data$locus_name)
  sample.list <- unique(as.character(mutation.data$sample))
  line.list <- unique(mutation.data$line)
  # create a count matrix for # of mutations for a given locus in each line
  count.matrix <- matrix(nrow=length(locus.list),
                         ncol=length(line.list),
                         dimnames = list(c(locus.list),
                                         c(line.list)))
  for(loci in locus.list){
    for(line in line.list){
      line.count <- length(mutation.data[which(mutation.data$locus_name == loci &
                                                 mutation.data$line == line),"transfer"])
      if(line.count > 0){
        count.matrix[loci,line] <- line.count
      } else {
        count.matrix[loci,line] <- 0
      }
    }
  }
  
  ## create a rank order matrix mutations based on the count matrix
  rank.matrix <- matrix(nrow=length(locus.list),
                        ncol=length(line.list),
                        dimnames = list(c(1:length(locus.list)),
                                        c(line.list)))
  for(line in line.list){
    tmp1 <- order(count.matrix[,line],decreasing = T)
    rank.matrix[,line] <- tmp1
  }
  
  ## create a data frame for ranks for each locus
  rank.order.df <- data.frame()
  for(row in 1:length(rank.matrix[,1])){
    locus = row.names(count.matrix)[row]
    rank.sum <- sum(rank.matrix[row,])
    rank.order.df <- rbind(rank.order.df, cbind(locus, rank.sum))
    #rank.order.df <- rank.order.df[sort(rank.order.df$rank.sum),]
  }
  return(as.character(rank.order.df$locus))
  
}


## order by line count ##
order.by.line <- function(input.data= mutation.data.df$mutation.data.noancestor, filter.single=T){
  mutation.data <- input.data
  line.count.df <- data.frame()
  locus.list <- unique(mutation.data$locus_name)
  for(loci in locus.list){
    line.count <- length(unique(mutation.data[which(mutation.data$locus_name == loci),"line"]))
    line.count.df <- rbind(line.count.df, cbind(locus=loci, line.count))
  }
  if(filter.single == T){
    line.count.df <- subset(line.count.df,subset = as.numeric(as.character(line.count)) > 1)
  } else {
    line.count.df <- line.count.df
  }
  line.count.df <- line.count.df[order(as.numeric(as.character(line.count.df$line.count)),decreasing = T),]
  return(as.character(line.count.df$locus))
}

## create ordering from heatmap
heatmap.ordering <- function(input.data=mutation.data.df$mutation.data.noancestor){
  heatmap.order.matrix <- matrix(nrow=length(unique(as.character(input.data$locus_name))), ncol = length(unique(as.character(input.data$transfer))),dimnames = list(c(unique(as.character(input.data$locus_name))), c(unique(as.character(input.data$transfer)))))
  for(row in unique(as.character(input.data$locus_name))){
    for(column in unique(as.character(input.data$transfer))){
      tmp1 <- length(input.data[which(as.character(input.data$locus_name) == row & as.character(input.data$transfer) == column),"line"])
      if(length(tmp1 > 0)){
        heatmap.order.matrix[row,column] = tmp1
      } else{
        heatmap.order.matrix[row,column] = 0
      }
    }
  }
  my.heatmap <- pheatmap(heatmap.order.matrix, cluster_cols = T, cluster_rows = T)
  heatmap.order <- my.heatmap$tree_row$labels[my.heatmap$tree_row$order]
  heatmap.order <- heatmap.order[c(1:7, 53:70, 8:52)]
  return(heatmap.order)
}

### heatmap with mutation data with no ancestor mutations
mutation.heatmap <- function(input.data=mutation.data.df$mutation.data.noancestor,
                             order=ordered.locus.heatmap, gene.list=NULL){
  
  if(!is.null(gene.list)){
    input.data <- subset(input.data, subset = input.data$locus_name %in% gene.list)
  } else {
    input.data <- input.data
  }
  
  
  plot.heatmap <- ggplot(data=input.data, 
                         aes(x=factor(transfer, levels = c( "15","45","76","118","152")),
                             y=factor(locus_name, levels=rev(order)),
                             fill = final_freq))
  plot.heatmap <- plot.heatmap + geom_raster(vjust=0, hjust=0)
  plot.heatmap <- plot.heatmap + facet_grid(.~line, scales = "fixed", drop = F )
  
  plot.heatmap <- plot.heatmap + scale_fill_distiller(type="seq", palette="YlGnBu", direction = 1, na.value = "gray90") 
  plot.heatmap <- plot.heatmap + scale_x_discrete(expand = c(0, 0.1)) + scale_y_discrete(expand = c(0, 0.1))
  plot.heatmap <- plot.heatmap + labs(x="Early Generations", y="Mutations")
  legend.plot.heatmap <- extractLegend(plot.heatmap)
  plot.heatmap <- plot.heatmap + theme(legend.position = "none", 
                                         axis.text.x = element_text(vjust=-0.5, angle = 90),
                                         axis.text.y = element_text(vjust = 1, hjust=0, margin=margin(t=2, r=0)),
                                         axis.ticks.x = element_blank(),
                                         axis.ticks.y = element_blank())
  
  
  print(plot.heatmap)
  return(list(legend.heatmap=legend.plot.heatmap, plot=plot.heatmap))
}


## Function to count mutation effects for each line
count.mutation.effects <- function(input.data=mutation.data.df$mutation.data.noancestor){
  # count and add mutation numbers based on impact
  mutation.effects= data.frame()
  for(line in unique(input.data$line)){
    for(transfer in unique(input.data$transfer)){
      mutations = input.data[which(input.data$line == line & input.data$transfer == transfer),"effect"]
      low = length(grep("LOW", mutations))
      high = length(grep("HIGH", mutations))
      modifier = length(grep("MODIFIER", mutations))
      moderate = length(grep("MODERATE", mutations))
      total = length(mutations)
      mutation.effects <- rbind(mutation.effects, cbind(
        line = line, low=as.numeric(as.character(low)), high=as.numeric(as.character(high)), moderate=as.numeric(as.character(moderate)), modifier=as.numeric(as.character(modifier)) ,total=as.numeric(as.character(total)), transfer=as.character(transfer)
      ))
    }
  }
  return(mutation.effects)
}

## Function to plot mutation effects barplots
mutation.barplot <- function(input.data=mutation.effects){
  egbar.dv <- melt(input.data, id.vars = c("line","transfer"), measure.vars = c("low", "modifier", "moderate", "high"))
  #remove data with 0 values
  egbar.dv <- egbar.dv[which(egbar.dv$value != 0),]
  
  dv.bar <- ggplot(egbar.dv, aes(x = factor(transfer, levels = c("15","45","76","118","152"))))
  dv.bar <- dv.bar + geom_bar(stat="identity", aes(y=as.numeric(as.character(value)), group=transfer, fill=variable))
  dv.bar <- dv.bar + theme_gray() + scale_fill_brewer(type = "div", palette = "RdYlGn", direction = -1)
  dv.bar <- dv.bar + facet_grid(.~ factor(line, levels = c("HA2","HA3","HE3","HR2","HS3","UA3","UE3","UR1","US1")),scales = "fixed")
  dv.bar <- dv.bar + theme(axis.text.x =element_blank(),
                           legend.direction = "horizontal",
                           legend.title = element_text(size=6), legend.key.size = unit(1,"cm"),
                           legend.text = element_text(size = 6)) +
    guides(fill = guide_legend(title.position = "left",
                               title.hjust = 0.5, 
                               label.position = "bottom", 
                               label.hjust = 0.5))
  legend.barplot <- extractLegend(dv.bar)
  dv.bar <- dv.bar + labs(x="", y="Mutations") + theme(legend.position = "none")
  print(dv.bar)
  return(list(legend.barplot=legend.barplot, plot=dv.bar))
}


## Find how many G-score mutations per line 
find.gscore.mutations <- function(gscores.data= NULL,
                            input.data= mutation.data.df$mutation.data.noancestor ){
  gscores <- read.delim(gscores.data, sep="\t")
  gscores <- subset(gscores, subset = observed > 1)
  gscores$new.locus <- sapply(gscores$locus, function(i) sub("_","", i))
  my.locuses <- unique(input.data$locus_name)
  
  my.list <- list()
  for(myline in unique(input.data$line)){
    my.list[[myline]] <- dim(unique(subset(input.data, subset = locus_name %in% gscores$new.locus & input.data$line == myline, select = locus_name)))[[1]][1]
  }
  
  }
  
#### Line specific heatmap panels
gene.list <- c("DVU2210", "DVU0150", "DVU0796", "DVU1833", "DVU1083", "DVU2396", "DVU2244", "DVU2609", "DVU0845", "DVU0672", "DVU0436", "DVU1788")
gene.list <- c("MMP1186", "MMP0378", "MMP0705", "MMP0986", "MMP1170", "MMP0033", "MMP0234", "MMP0939", "MMP0694", "MMP0146", "MMP0643", "IG_1439720", "MMP0026", "MMP0466", "MMP1557", "MMP0952", "MMP1077", "MMP1479", "MMP1180", "IG_646987", "MMP1363", "MMP1720")

### heatmap with mutation data with no ancestor mutations
mutation.heatmap.linespecific <- function(input.data=mutation.data.df$mutation.data.noancestor,
                             gene.list=gene.list){
  
  if(!is.null(gene.list)){
    input.data <- subset(input.data, subset = input.data$locus_name %in% gene.list)
  } else {
    input.data <- input.data
  }
  
  
  plot.heatmap <- ggplot(data=input.data, 
                         aes(x=factor(transfer, levels = c( "15","45","76","118","152")),
                             y=factor(locus_name, levels=rev(gene.list)),
                             #y=locus_name,
                             fill = final_freq))
  plot.heatmap <- plot.heatmap + geom_raster(vjust=0, hjust=0)
  plot.heatmap <- plot.heatmap + facet_grid(.~line, scales = "fixed", drop = F )
  
  plot.heatmap <- plot.heatmap + scale_fill_distiller(type="seq", palette="YlGnBu", direction = 1, na.value = "gray90") 
  plot.heatmap <- plot.heatmap + scale_x_discrete(expand = c(0, 0.1)) + scale_y_discrete(expand = c(0, 0.1))
  plot.heatmap <- plot.heatmap + labs(x="Early Generations", y="Mutations")
  legend.plot.heatmap <- extractLegend(plot.heatmap)
  plot.heatmap <- plot.heatmap + theme(legend.position = "none", 
                                       axis.text.x = element_text(vjust=-0.5, angle = 90),
                                       axis.text.y = element_text(vjust = 1, hjust=0, margin=margin(t=2, r=0)),
                                       axis.ticks.x = element_blank(),
                                       axis.ticks.y = element_blank())
  
  
  print(plot.heatmap)
  return(list(legend.heatmap=legend.plot.heatmap, plot=plot.heatmap))
}

#### Dv heatmap 
## Load and process mutations data ##
mutation.data.df <- get.modify.data(mutation.file = "~/Documents/GitHub/evolution-of-syntrophy/data/dvh_mutations_allsamples_attributes_3112020.txt")

## Order based on number of mutations ##
ordered.locus.counts <- order.by.mutations()

## Order based on Gscores ##
ordered.locus.gscore <- order.by.gscore(gscores.data ="../Figures/Figure 2/simulated_Gscores_Dv.txt")

## Order based on heatmap clustering ##
ordered.locus.heatmap <- heatmap.ordering()

## Order based on ranking ##
ordered.locus.ranking <- order.by.rank()

## Order based on lines ##
ordered.locus.line <- order.by.line(filter.single = T)

## Plot heatmap based selected ordering ##
ordered.locus.heatmap <- mutation.heatmap(order=ordered.locus.line, gene.list = ordered.locus.line)

## Count mutation effects ##
mutation.effects <- count.mutation.effects()

## Plot mutation barplots ##
ordered.locus.barplot <- mutation.barplot(input.data = mutation.effects)

## Plot all braphs combined ##
grid.arrange(
  arrangeGrob(ordered.locus.barplot$plot, ordered.locus.heatmap$plot, ncol=1, nrow=2,widths=c(60),heights=c(10,50),padding = c(0,0,0,50)), 
  arrangeGrob(ordered.locus.barplot$legend.barplot, ordered.locus.heatmap$legend.heatmap, ncol=1,nrow=3), 
  ncol=2, nrow=1)



#### Mm heatmap 
## Load and process mutations data ##
mutation.data.df <- get.modify.data(mutation.file = "~/Documents/GitHub/evolution-of-syntrophy/data/mmp_mutations_allsamples_attributes_3112020.txt")

## Order based on number of mutations ##
ordered.locus.counts <- order.by.mutations()

## Order based on Gscores ##
ordered.locus.gscore <- order.by.gscore(gscores.data ="../Figures/Figure 2/simulated_Gscores_Mm.txt")

## Order based on heatmap clustering ##
ordered.locus.heatmap <- heatmap.ordering()

## Order based on ranking ##
ordered.locus.ranking <- order.by.rank()

## Order based on lines ##
ordered.locus.line <- order.by.line(filter.single = T)

## Plot heatmap based selected ordering ##
ordered.locus.heatmap <- mutation.heatmap(order=ordered.locus.line, gene.list = ordered.locus.line)

## Count mutation effects ##
mutation.effects <- count.mutation.effects()

## Plot mutation barplots ##
ordered.locus.barplot <- mutation.barplot(input.data = mutation.effects)

## Plot all braphs combined ##
grid.arrange(
  arrangeGrob(ordered.locus.barplot$plot, ordered.locus.heatmap$plot, ncol=1, nrow=2,widths=c(60),heights=c(10,50),padding = c(0,0,0,50)), 
  arrangeGrob(ordered.locus.barplot$legend.barplot, ordered.locus.heatmap$legend.heatmap, ncol=1,nrow=3), 
  ncol=2, nrow=1)


### heatmap with mutation data with no ancestor mutations for HS3 sweep.
HS3.sweep.heatmap <- function(plot.line="HS3", gene.list=gene.list){
  ## gene list that had selective sweep
  gene.list <- c("DVU0799","DVU0597","DVU0797","DVU0013","DVU0001","MMP1511","MMP1362","MMP1255","MMP1611",
                 "DVU2394","DVU2451","DVU1862","IG_184033","MMP0466","MMP1557","MMP0952", "MMP1077","MMP1479",
                 "DVU2609","DVU0845","DVU0672")
  ## Get DV mutations specific to Early gene and 1000-gen for HS3
  dv.mutations <- load.mutation.data(org = "dvh",sample_type = "all")
  dv.mutations <- dv.mutations$mutation.data.noancestor
  dv.mutations <- subset(dv.mutations, line==plot.line & experiment %in% c("Early-gen", "1000-gen"))
  ## Get Mm mutations specific to Early gene and 1000-gen for HS3
  mm.mutations <- load.mutation.data(org = "mmp",sample_type = "all")
  mm.mutations <- mm.mutations$mutation.data.noancestor
  mm.mutations <- subset(mm.mutations, line==plot.line & experiment %in% c("Early-gen", "1000-gen"))
  ## Combine mutations
  all.mutations <- rbind(dv.mutations, mm.mutations)
  ## Filter data for selected genes
  input.data <- subset(all.mutations, subset = all.mutations$locus_name %in% gene.list)
  ## Plot heatmap
  plot.heatmap <- ggplot(data=input.data, 
                         aes(x=factor(transfer, levels = c( "15","45","76","118","152")),
                             y=factor(locus_name, levels=rev(gene.list)),
                             #y=locus_name,
                             fill = final_freq))
  plot.heatmap <- plot.heatmap + geom_raster(vjust=0, hjust=0)
  plot.heatmap <- plot.heatmap + facet_grid(.~line, scales = "fixed", drop = F )
  
  plot.heatmap <- plot.heatmap + scale_fill_distiller(type="seq", palette="YlGnBu", direction = 1, na.value = "gray90") 
  plot.heatmap <- plot.heatmap + scale_x_discrete(expand = c(0, 0.1)) + scale_y_discrete(expand = c(0, 0.1))
  plot.heatmap <- plot.heatmap + labs(x="Early Generations", y="Mutations")
  legend.plot.heatmap <- extractLegend(plot.heatmap)
  plot.heatmap <- plot.heatmap + theme(legend.position = "none", 
                                       axis.text.x = element_text(vjust=-0.5, angle = 90),
                                       axis.text.y = element_text(vjust = 1, hjust=0, margin=margin(t=2, r=0)),
                                       axis.ticks.x = element_blank(),
                                       axis.ticks.y = element_blank())
  
  
  print(plot.heatmap)
  return(list(legend.heatmap=legend.plot.heatmap, plot=plot.heatmap))
}

hs3.heatmap <- HS3.sweep.heatmap(plot.line = "HS3")





