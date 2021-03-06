## calculate G-score based on the https://www.nature.com/articles/nature18959#s1
# Li: length of protein-coding gene
# Ni: the number of independent nonsynonymous mutations observed in that gene across all clones
# Ltot: summed length of all protein-coding genes in the ancestral genome
# Ntot: summed total number of nonsynonymous mutations in these genes
# Ei: expected number of non-synonymous mutations for a given gene
# Ei = Ntot(Li/Ltot)
#  Goodness of fit score: Gi = 2Niloge(Ni/Ei)
library('ggplot2'); library('reshape2');library(gridExtra);library(gplots)
source("~/Google Drive File Stream/My Drive/R_Scripts/extractLegend.R")

# load mutations data
#DvH
mut.data.dv <- read.delim("/Volumes/omics4tb/sturkarslan/syntrophy_raw_sequences/dvh_mutations_allsamples_attributes_3112020.txt", header=T, sep="\t",)
mut.data.dv$organism <- "Desulfovibrio"
length(mut.data.dv$variant_id) # 3911
length(unique(mut.data.dv$variant_id)) # 278

# # remove TG_HS3 data
# mut.data.dv <- mut.data.dv[which(mut.data.dv$sample != "TG_HS3"),]
# length(mut.data.dv$variant_id) # 3446
# length(unique(mut.data.dv$variant_id)) # 319

# get 1000-gen data
mutations.1000 <- mut.data.dv[which(mut.data.dv$experiment == "1000-gen"),]
length(mutations.1000$variant_id) # 557
length(unique(mutations.1000$variant_id)) # 159

# get ancestor mutations
ancestor.mutations <- mut.data.dv[which(mut.data.dv$sample == "AN_Coculture-Ancestor"),]
ancestor.mutated.genes <- unique(ancestor.mutations$gene_id)
ancestor.mutation.ids <- unique(ancestor.mutations$variant_id)
length(ancestor.mutations$variant_id) # 33
length(ancestor.mutated.genes) #16

# Filter mutation data for non-ancestral mutations
allmutation.ids <- unique(mutations.1000$variant_id)
nonancestor.mutation.ids <- setdiff(allmutation.ids, ancestor.mutation.ids)
length(nonancestor.mutation.ids) # 126

mutations.1000.na <- data.frame()
for(mutation.id in nonancestor.mutation.ids){
  mutations.1000.na <- rbind(mutations.1000.na, mutations.1000[which(mutations.1000$variant_id == mutation.id),])
}
length(mutations.1000.na$variant_id) # 138

# mutations per line
mutations.line <- data.frame()
for(line in unique(mutations.1000.na$sample)){
  count <- length(mutations.1000.na[which(mutations.1000.na$sample == line),"variant_id"])
  mutations.line <- rbind(mutations.line, cbind(sample=line, mutation.count=count))
}
# sample mutation.count
# TG_HA2             13
# TG_US1             21
# TG_HA3             17
# TG_HR2             17
# TG_HR1             23
# TG_UA3             16
# TG_UE3             11
# TG_HE2             19
# TG_UE2             12
# TG_HE3             15
# TG_UA2             13
# TG_UR1             13
# TOTAL             190
# mean              15.83
mean.line.mutations <- mean(as.numeric(as.character(mutations.line$mutation.count)))

# load dvh genome features file
dvh.genome <- read.delim("~/Google Drive File Stream/My Drive/Manuscripts/Syntrophy-SingleCell/GenomeFiles/GCF_000195755.1_ASM19575v1_feature_table.txt",
                         header = T, sep="\t", stringsAsFactors=F)
# collect only genes
dvh.genes <- dvh.genome[which(dvh.genome$X..feature == "gene"),]
dvh.cds <- dvh.genome[which(dvh.genome$X..feature == "CDS" | dvh.genome$X..feature == "tRNA" | dvh.genome$X..feature == "rRNA" | dvh.genome$X..feature == "misc_RNA"),]


##### all mutations except low and modifier forDv
# calculate total coding length (Ltot)
Ltot <- sum(dvh.genes$feature_interval_length) #3313981
frameshift.mutations <- mutations.1000.na[grep("FRAME_SHIFT", mutations.1000.na$mutation),]
length(frameshift.mutations$variant_id) #45
syncoding.mutations <- mutations.1000.na[which(mutations.1000.na$mutation == "SYNONYMOUS_CODING"),]
length((syncoding.mutations$variant_id)) # 4

# get all mutations except synonymous and intergenic
all.mutations <- mutations.1000.na[which(mutations.1000.na$effect == "HIGH" | mutations.1000.na$effect == "MODERATE"),]
Ntot <- length(all.mutations$variant_id) # 124

all.table <- data.frame()
for(gene in dvh.genes$GeneID){
  locus <- dvh.genes[which(dvh.genes$GeneID == gene), "locus_tag"]
  locus <- sub("DVUA", "DVKA", locus)
  locus <- sub("DVU", "DVU_", locus)
  locus <- sub("DVKA", "DVUA", locus)
  name <- dvh.cds[which(dvh.cds$GeneID == gene), "name"]
  symbol <- dvh.genes[which(dvh.genes$GeneID == gene), "symbol"]
  
  Li <- dvh.genes[which(dvh.genes$GeneID == gene),"feature_interval_length"]
  Ni <- length(all.mutations[which(all.mutations$gene_id == locus),"variant_id"])
  Nframeshift <- length(frameshift.mutations[which(frameshift.mutations$gene_id == locus),"variant_id"])
  Nsynonymous <- length(syncoding.mutations[which(syncoding.mutations$gene_id == locus),"variant_id"])
  # calculate expected number of mutations
  Ei <- Ntot*(Li/Ltot)
  # calculate GScore
  Gscore <- 2*Ni*log1p(Ni/Ei)
  # # add ancestor annotation
  # if(locus %in% ancestor.mutated.genes){
  #   ancestor = "ancestral"
  # }else{
  #   ancestor = "non-ancestral"
  # }
  all.table <- rbind(all.table, cbind(locus=locus,
                                      symbol= if(length(symbol !=0)) {paste(symbol)}else{paste("")},
                                      name=if(length(name !=0)) {paste(name)}else{paste("")},
                                      Length=Li,
                                      observed=Ni,
                                      expected=Ei,
                                      Gscore=Gscore,
                                      Nframeshift=Nframeshift,
                                      Nsynonymous=Nsynonymous
                                      ))
  cat(gene, locus, Li, Ni, Ei, "\n")
}

# sort table based on Gscores
all.table <- all.table[order(as.numeric(as.character(all.table$Gscore)), decreasing = T ),]
sum.gscores <- sum(as.numeric(as.character(all.table$Gscore))) # 1092.617

# write results into a table
write.table(all.table, file="Gscore_high_moderate_dvh.txt", sep="\t")

#### mutations for individual genes
ind.data.dv.tmp1 <- mutations.1000.na[,c("variant_id", "sample","position", "effect", "gene_id", "freq")]
# filter low impact mutations
ind.data.dv.tmp1 <- ind.data.dv.tmp1[which(ind.data.dv.tmp1$effect != "LOW"),]

ind.data.dv <- data.frame()
for(row in 1:length(ind.data.dv.tmp1$variant_id)){
  variant <- ind.data.dv.tmp1[row,"variant_id"]
  sample <-  ind.data.dv.tmp1[row,"sample"]
  if(variant %in% ancestor.mutations){
    ancestor <- "ancestral"
  }else{
    ancestor <- "non-ancestral"
  }
  ind.data.dv <- rbind(ind.data.dv, cbind(ind.data.dv.tmp1[row,], ancestor = ancestor))
}
ind.data.dv.noanc <- ind.data.dv[which(ind.data.dv$ancestor !="ancestral"),]

## prepare plotting for dotplot of individual genes
# get the G-score list
#granks <- read.delim("~/Google Drive File Stream/My Drive/Manuscripts/Syntrophy-SingleCell/Data/Gscore_high_moderate_dvh.txt", sep="\t", header=T)
# filter gransks for ancetral mutations
#granks <- granks[which(granks$ancestor != "ancestral"),]
# plot GScores
gsplot.data <- all.table[which(as.numeric(as.character(all.table$observed)) > 0),]
gsplot.data$gene <- sapply(gsplot.data$locus, function(i) sub("_","", i))
gsplot.data <- subset(gsplot.data, subset = as.numeric(as.character(observed)) > 1)
gsplot.data1 <- gsplot.data[1:12,]
gsplot.data2 <- gsplot.data[13:24,]
#gsplot.data2 <- gsplot.data[c(13:23,24),]

# # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # #
# SET 1
# use top level 15 genes for the plotting of the first set
top15 = paste("(",paste(gsplot.data1$gene, sep = "", collapse = "|"), ")",sep = "", collapse = "")
set1 <- ind.data.dv[grep(top15, ind.data.dv$variant_id),]
set1 <- set1[which(set1$ancestor == "non-ancestral"),]
set1$gene <- sapply(set1$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])
set1$value <- as.numeric(as.character(sub("%", "",set1$freq)))
set1 <- set1[!is.na(set1$value),]
facet.genes1 <- paste("'",gsplot.data1$gene, "'",sep = "", collapse = ",")

## Dotplot for mutation coordinares across lines
pset1 <- ggplot(data=set1, aes(factor(sample)))
pset1 <- pset1 + geom_point(shape=21, color="black", aes(y=position, fill=factor(effect), size=as.numeric(value)), alpha=0.8) + theme_gray()
pset1 <- pset1 + theme(legend.direction = "horizontal",
                       legend.title = element_text(size=6), 
                       legend.key.size = unit(1,"cm"),
                       legend.text = element_text(size = 6))
pset1 <- pset1 + scale_y_continuous(expand = c(0.5,0.5))
pset1 <- pset1 + facet_grid(factor(gene, levels=gsplot.data1$gene)~., scales = "free")
#pset1 <- pset1 + facet_grid(factor(gene, levels=c('DVU0799','DVU2776','DVU1283','DVU1260','DVU2394','DVU1295','DVU0846','DVU2451','DVU1862','DVU2894','DVU0001','DVU1092','DVU0597'))~., scales = "free")
pset1 <- pset1 + scale_fill_manual(values=c("#d62027", "#f9ac61"))
pset1 <- pset1 + guides(fill = guide_legend(title.position = "left",
                                            title.hjust = 0.5, 
                                            label.position = "bottom", 
                                            label.hjust = 0.5))
# extract legend to plot later by using extract legend function
legend.pset1 <-extractLegend(pset1)
pset1 <- pset1 + theme(legend.position = "none", axis.text.x = element_text(size = 8)) + labs(x="Evolution Line", y="Mutations")

# # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # #
# SET2
# use top level 15-30 genes for the plotting of the second set
top30 = paste("(",paste(gsplot.data2$gene, sep = "", collapse = "|"), ")",sep = "", collapse = "")
set12 <- ind.data.dv[grep(top30, ind.data.dv$variant_id),]
set12 <- set12[which(set12$ancestor == "non-ancestral"),]
set12$gene <- sapply(set12$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])
set12$value <- as.numeric(as.character(sub("%", "",set12$freq)))
set12 <- set12[!is.na(set12$value),]
facet.genes2 <- paste("'",gsplot.data2$gene, "'",sep = "", collapse = ",")

## Dotplot for mutation coordinates across lines
pset12 <- ggplot(data=set12, aes(factor(sample)))
pset12 <- pset12 + geom_point(shape=21, color="black", aes(y=position, fill=factor(effect), size=as.numeric(value)),alpha=0.8) + theme_gray()
pset12 <- pset12 + theme(legend.direction = "horizontal",
                       legend.title = element_text(size=6), 
                       legend.key.size = unit(1,"cm"),
                       legend.text = element_text(size = 6))
pset12 <- pset12 + scale_y_continuous(expand = c(0.5,0.5))
#pset12 <- pset12 + facet_grid(factor(gene, levels=c('DVU0013','DVU0847','DVU1214','DVU2395','DVU0797','DVU0150','DVU2305','DVU3227','DVU0436','DVU0876','DVU2210','DVU1833'))~., scales = "free")
pset12 <- pset12 + facet_grid(factor(gene, levels=gsplot.data2$gene)~., scales = "free")
pset12 <- pset12 + scale_fill_manual(values=c("#d62027", "#f9ac61"))
pset12 <- pset12 + guides(fill = guide_legend(title.position = "left",
                                            title.hjust = 0.5, 
                                            label.position = "bottom", 
                                            label.hjust = 0.5))
pset12 <- pset12 + theme(legend.position = "none", axis.text.x = element_text(size = 8)) + labs(x="Evolution Line", y="Mutations")


# SET 1
# Plotting Text values next to dotplot
# prepare data for plotting of text values of name, gscore and # of mutations for individual genes 1st set
set2 <- ind.data.dv
set2 <- set2[!is.na(set2$freq),]
#set2$coord <- sapply(set2$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][3])
set2$gene<- sapply(set2$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])

# get counts of mutations
line.counts <- data.frame()
for(mgene in unique(set2$gene)){
  lines <- set2[which(set2$gene == mgene),"sample"]
  count <- length(lines)
  line.counts <- rbind(line.counts, cbind(gene=mgene, linecount=count, lines=paste(lines, sep = "", collapse = ":")))
}

## 1set of gene plots
top15.1 = paste("(",paste(gsplot.data1$locus, sep = "", collapse = "|"), ")",sep = "", collapse = "")
pset2.data <- line.counts[grep(top15.1, line.counts$gene),]
pset2 <- ggplot(data=gsplot.data1, aes(x=1))
#pset2 <- pset2 + geom_bar(aes(y=as.numeric(as.character(linecount)), fill=gene),stat = "identity", width = 0.1) + theme_gray()
pset2 <- pset2 + lims(y=c(0,16), x=c(0,3))
#gsplot.data1$gene <- gsplot.data1$locus
pset2 <- pset2 + geom_label(data=gsplot.data1, aes(x=2.5, y=14, label=round(as.numeric(as.character(Gscore)), digits = 1), color="red"),size=3)
pset2 <- pset2 + geom_label(data=gsplot.data1, aes(x=2.5, y=3, label=observed), color="blue", size=3)
pset2 <- pset2 + geom_text(data=gsplot.data1, aes(x=0.5, y=10, label=name), color="black", size=2, hjust="center",vjust="center")
facet.genes1 <- paste("'",gsplot.data1$locus, "'",sep = "", collapse = ",")

# extract legend to plot later
legend.pset2 <- extractLegend(pset2)
pset2 <- pset2 + facet_grid(factor(gsplot.data1$gene, levels=gsplot.data1$gene)~.,scales = "free")  + coord_flip()
#pset2 <- pset2 + facet_grid(factor(gene, levels=c('DVU_0799','DVU_2776','DVU_1283','DVU_1260','DVU_2394','DVU_1295','DVU_0846','DVU_2451','DVU_1862','DVU_2894','DVU_0001','DVU_1092','DVU_0597'))~.,)  + coord_flip()
pset2 <- pset2 + guides(fill = guide_legend(title.position = "top",
                                            title.hjust = 0.5, 
                                            label.position = "bottom", 
                                            label.hjust = 0.5))
pset2 <- pset2 + scale_fill_manual(values = rep("gray50", 12))
pset2 <- pset2 + theme(axis.text.y = element_blank(),
                       axis.title.y = element_blank(),
                       axis.title.x = element_text(),
                       panel.grid = element_blank(),
                       axis.ticks = element_blank(),
                       strip.background = element_blank(),
                       strip.text.x = element_blank(),
                       legend.position = "none") + labs(x="Evolution Line", y="Mutations")
pset2 <- pset2 + theme(legend.position = "none", axis.text.x = element_text(size = 8)) + labs(x="Evolution Line", y="Mutations")


# SET 2
# Plotting Text values next to dotplot
# prepare data for plotting of text values of name, gscore and # of mutations for individual genes 1st set
set21 <- ind.data.dv
set21 <- set21[!is.na(set21$freq),]
#set2$coord <- sapply(set2$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][3])
set21$gene<- sapply(set21$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])

# get counts of mutations
line.counts <- data.frame()
for(mgene in unique(set21$gene)){
  lines <- set21[which(set21$gene == mgene),"sample"]
  count <- length(lines)
  line.counts <- rbind(line.counts, cbind(gene=mgene, linecount=count, lines=paste(lines, sep = "", collapse = ":")))
}

## 2nd set of gene plots
top30.1 = paste("(",paste(gsplot.data2$locus, sep = "", collapse = "|"), ")",sep = "", collapse = "")
pset21.data <- line.counts[grep(top30.1, line.counts$gene),]
pset21 <- ggplot(data=gsplot.data2, aes(x=1))
#pset21 <- pset21 + geom_bar(aes(y=as.numeric(as.character(linecount)), fill=gene),stat = "identity", width = 0.1) + theme_gray()
pset21 <- pset21 + lims(y=c(0,16), x=c(0,3))
#gsplot.data2$gene <- gsplot.data2$locus
pset21 <- pset21 + geom_label(data=gsplot.data2, aes(x=2.5, y=14, label=round(as.numeric(as.character(Gscore)), digits = 1), color="red"),size=3)
pset21 <- pset21 + geom_label(data=gsplot.data2, aes(x=2.5, y=3, label=observed), color="blue", size=3)
pset21 <- pset21 + geom_text(data=gsplot.data2, aes(x=0.5, y=10, label=name), color="black", size=2, hjust="center",vjust="center")

facet.genes2 <- paste("'",gsplot.data2$locus, "'",sep = "", collapse = ",")

# extract legend
legend.pset21 <- extractLegend(pset21)
pset21 <- pset21 + facet_grid(factor(gsplot.data2$gene, levels=gsplot.data2$gene)~.)  + coord_flip()
#pset21 <- pset21 + facet_grid(factor(gene, levels=c('DVU_0013','DVU_0847','DVU_1214','DVU_2395','DVU_0797','DVU_0150','DVU_2305','DVU_3227','DVU_0436','DVU_0876','DVU_2210','DVU_1833'))~.,)  + coord_flip()
pset21 <- pset21 + guides(fill = guide_legend(title.position = "top",
                                            title.hjust = 0.5, 
                                            label.position = "bottom", 
                                            label.hjust = 0.5))

pset21 <- pset21 + scale_fill_manual(values = rep("gray50", 12))
pset21 <- pset21 + theme(axis.text.y = element_blank(),
                       axis.title.y = element_blank(),
                       axis.title.x = element_text(),
                       panel.grid = element_blank(),
                       axis.ticks = element_blank(),
                       strip.background = element_blank(),
                       strip.text.x = element_blank(),
                       legend.position = "none") + labs(x="Evolution Line", y="Mutations")
pset21 <- pset21 + theme(legend.position = "none", axis.text.x = element_text(size = 8)) + labs(x="Evolution Line", y="Mutations")

## Put all plots together
grid.arrange(
  arrangeGrob(pset1, pset2, ncol=2, nrow=1,widths=c(70,30),heights=c(60)), 
  arrangeGrob(pset12,pset21, ncol=2,nrow=1, widths=c(70,30), heights=c(60)), 
  arrangeGrob(legend.pset1, legend.pset2, legend.pset21,ncol=1,nrow=3,widths=c(10)), 
  ncol=3, nrow=1, widths=c(40,40,20))




