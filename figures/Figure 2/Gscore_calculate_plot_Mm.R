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
mut.data.mm <- read.delim("/Volumes/omics4tb/sturkarslan/syntrophy_raw_sequences/mmp_mutations_allsamples_attributes_3112020.txt", header=T, sep="\t")
mut.data.mm$organism <- "Methanococcus"
length(mut.data.mm$variant_id) # 941
length(unique(mut.data.mm$variant_id)) # 142

# remove TG_HS3 data
# mut.data.mm <- mut.data.mm[which(mut.data.mm$sample != "TG_HS3"),]
# length(mut.data.mm$variant_id) # 935
# length(unique(mut.data.mm$variant_id)) # 141

# get 1000-gen data
mutations.1000.mm <- mut.data.mm[which(mut.data.mm$experiment == "1000-gen"),]
length(mutations.1000.mm$variant_id) # 181
length(unique(mutations.1000.mm$variant_id)) # 103

# get ancestor mutations
ancestor.mutations.mm <- mut.data.mm[which(mut.data.mm$sample == "AN_Coculture-Ancestor"),]
ancestor.mutated.genes.mm <- unique(ancestor.mutations.mm$gene_id)
ancestor.mutation.ids.mm <- unique(ancestor.mutations.mm$variant_id)
length(ancestor.mutations.mm$variant_id) # 5
length(ancestor.mutated.genes.mm) #5

# Filter mutation data for non-ancestral mutations
allmutation.ids.mm <- unique(mutations.1000.mm$variant_id)
nonancestor.mutation.ids.mm <- setdiff(allmutation.ids.mm, ancestor.mutation.ids.mm)
mutations.1000.mm.na <- data.frame()
for(mutation.id in nonancestor.mutation.ids.mm){
  mutations.1000.mm.na <- rbind(mutations.1000.mm.na, mutations.1000.mm[which(mutations.1000.mm$variant_id == mutation.id),])
}
length(mutations.1000.mm.na$variant_id) # 121

# mutations per line
mutations.line.mm <- data.frame()
for(line in unique(mutations.1000.mm.na$sample)){
  count <- length(mutations.1000.mm.na[which(mutations.1000.mm.na$sample == line),"variant_id"])
  mutations.line.mm <- rbind(mutations.line.mm, cbind(sample=line, mutation.count=count))
}
# sample mutation.count
# TG_HA2              9
# TG_HE2             11
# TG_HR2             12
# TG_UE2              8
# TG_US1             11
# TG_HA3             10
# TG_HE3              9
# TG_HR1              8
# TG_UA2              9
# TG_UA3              8
# TG_UE3              9
# TG_UR1             17
# Total              121
# mean              10.08
mean.line.mutations.mm <- mean(as.numeric(as.character(mutations.line.mm$mutation.count)))

# load dvh genome features file
mmp.genome <- read.delim("~/Google Drive File Stream/My Drive/Manuscripts/Syntrophy-SingleCell/GenomeFiles/GCA_000011585.1_ASM1158v1_feature_table.txt",
                         header = T, sep="\t", stringsAsFactors=F)
# collect only genes
mmp.genes <- mmp.genome[which(mmp.genome$X..feature == "gene"),]
mmp.cds <- mmp.genome[which(mmp.genome$X..feature == "CDS" | mmp.genome$X..feature == "tRNA" | mmp.genome$X..feature == "rRNA" | mmp.genome$X..feature == "misc_RNA"),]

# # get ancestor mutations
# ancestor.mutations <- mut.data.mm[which(mut.data.mm$sample == "AN_Coculture-Ancestor"),]
# ancestor.mutated.genes <- unique(ancestor.mutations$gene_id)
# 
# ## collect 1000 gen data for DvH
# tg.data.mm <- mut.data.mm[grep(c("(TG_)"), mut.data.mm$sample),]
# tg.mut.data.mm = data.frame()
# for(tg in unique(tg.data.mm$sample)){
#   # get list of mutatons for given clone
#   tg.mutations <- tg.data.mm[which(tg.data.mm$sample == tg & tg.data.mm$organism == "Methanococcus"), "effect"]
#   low = length(grep("LOW", tg.mutations))
#   high = length(grep("HIGH", tg.mutations))
#   modifier = length(grep("MODIFIER", tg.mutations))
#   moderate = length(grep("MODERATE", tg.mutations))
#   total = length(tg.mutations)
#   tg.mut.data.mm <- rbind(tg.mut.data.mm, 
#                           cbind(tg, 
#                                 low=as.numeric(as.character(low)), 
#                                 high=as.numeric(as.character(high)),
#                                 modifier=as.numeric(as.character(modifier)),
#                                 moderate=as.numeric(as.character(moderate)),
#                                 total=as.numeric(as.character(total)),
#                                 epd="NA",
#                                 cloneid="NA",
#                                 organism="Methanococcus"))
# }


##### all mutations except low and modifier for Mmp
# calculate total coding length (Ltot)
Ltot <- sum(mmp.genes$feature_interval_length) #1494525
frameshift.mutations <- mutations.1000.mm[grep("FRAME_SHIFT", mutations.1000.mm$mutation),]
length(frameshift.mutations$variant_id) #58
syncoding.mutations <- mutations.1000.mm[which(mutations.1000.mm$mutation == "SYNONYMOUS_CODING"),]
length((syncoding.mutations$variant_id)) # 9

all.mutations.mm <- mutations.1000.mm.na[which(mutations.1000.mm.na$effect == "HIGH" | mutations.1000.mm.na$effect == "MODERATE"),]
Ntot <- length(all.mutations.mm$variant_id) #103

all.table.mm <- data.frame()
for(gene in mmp.genes$locus_tag){
  locus <- gene
  name <- mmp.cds[which(mmp.cds$locus_tag == gene), "name"]
  symbol <- mmp.genes[which(mmp.genes$locus_tag == gene), "symbol"]
  
  Li <- mmp.genes[which(mmp.genes$locus_tag == gene),"feature_interval_length"]
  Ni <- length(all.mutations.mm[which(all.mutations.mm$gene_id == locus),"variant_id"])
  Nframeshift <- length(frameshift.mutations[which(frameshift.mutations$gene_id == locus),"variant_id"])
  Nsynonymous <- length(syncoding.mutations[which(syncoding.mutations$gene_id == locus),"variant_id"])
  # calculate expected number of mutations
  Ei <- Ntot*(Li/Ltot)
  # calculate GScore
  Gscore <- 2*Ni*log1p(Ni/Ei)
  # add ancestor annotation
  if(locus %in% ancestor.mutated.genes.mm){
    ancestor = "ancestral"
  }else{
    ancestor = "non-ancestral"
  }
  all.table.mm <- rbind(all.table.mm, cbind(locus=locus,
                                            symbol= if(length(symbol !=0)) {paste(symbol)}else{paste("")},
                                            name=if(length(name !=0)) {paste(name)}else{paste("")},
                                            Length=Li,
                                            observed=Ni,
                                            expected=Ei,
                                            Gscore=Gscore,
                                            Nframeshift=Nframeshift,
                                            Nsynonymous=Nsynonymous,
                                            ancestor = ancestor))
  cat(gene, locus, Li, Ni, Ei, "\n")
}
all.table.mm <- all.table.mm[order(as.numeric(as.character(all.table.mm$Gscore)), decreasing = T ),]
# filter ancestral mutations
all.table.mm <- all.table.mm[which(all.table.mm$ancestor != "ancestral"),]
sum.gscores <- sum(as.numeric(as.character(all.table.mm$Gscore))) # 839.5635

# write results into a table
write.table(all.table.mm, file="~/Google Drive File Stream/My Drive/Manuscripts/Syntrophy-SingleCell/Data/Gscore_frameshift_moderate_mmp.txt", sep="\t", row.names = F)


######
#### mutations for individual genes
ind.data.mm.tmp1 <- mutations.1000.mm[,c("variant_id", "sample","position", "effect", "gene_id", "freq")]
# filter low impact mutations
ind.data.mm.tmp1 <- ind.data.mm.tmp1[which(ind.data.mm.tmp1$effect != "LOW"),]

ind.data.mm <- data.frame()
for(row in 1:length(ind.data.mm.tmp1$variant_id)){
  variant <- ind.data.mm.tmp1[row,"variant_id"]
  sample <-  ind.data.mm.tmp1[row,"sample"]
  if(variant %in% ancestor.mutations.mm){
    ancestor <- "ancestral"
  }else{
    ancestor <- "non-ancestral"
  }
  ind.data.mm <- rbind(ind.data.mm, cbind(ind.data.mm.tmp1[row,], ancestor = ancestor))
}
ind.data.mm.noanc <- ind.data.mm[which(ind.data.mm$ancestor !="ancestral"),]

## prepare plotting for dotplot of individual genes
# get the G-score list
#granks <- read.delim("~/Google Drive File Stream/My Drive/Manuscripts/Syntrophy-SingleCell/Data/Gscore_high_moderate_dvh.txt", sep="\t", header=T)
# filter gransks for ancetral mutations
#granks <- granks[which(granks$ancestor != "ancestral"),]
# plot GScores
gsplot.data <- all.table.mm[which(as.numeric(as.character(all.table.mm$observed)) > 0),]
gsplot.data$gene <- sapply(gsplot.data$locus, function(i) sub("_","", i))
gsplot.data1 <- gsplot.data[1:12,]
gsplot.data2 <- gsplot.data[c(13:24),] #also add the MMP1718 to make sure to draw all lines. otherwise it display only 9 lines

# # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # #
# SET 1
# use top level 15 genes for the plotting of the first set
top15 = paste("(",paste(gsplot.data1$gene, sep = "", collapse = "|"), ")",sep = "", collapse = "")
set1 <- ind.data.mm[grep(top15, ind.data.mm$variant_id),]
set1 <- set1[which(set1$ancestor == "non-ancestral"),]
set1$gene <- sapply(set1$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])
set1$value <- as.numeric(as.character(sub("%", "",set1$freq)))
set1 <- set1[!is.na(set1$value),]
facet.genes1 <- paste("'",gsplot.data1$gene, "'",sep = "", collapse = ",")

## Dotplot for mutation coordinares across lines
pset1 <- ggplot(data=set1, aes(factor(sample)))
pset1 <- pset1 + geom_point(shape=21, color="black",aes(y=position, fill=factor(effect), size=as.numeric(value)), alpha=0.8) + theme_gray()
pset1 <- pset1 + theme(legend.direction = "horizontal",
                       legend.title = element_text(size=6), 
                       legend.key.size = unit(1,"cm"),
                       legend.text = element_text(size = 6))
pset1 <- pset1 + scale_y_continuous(expand = c(0.5,0.5))
pset1 <- pset1 + facet_grid(factor(gene, levels=gsplot.data1$gene)~., scales = "free")
#pset1 <- pset1 + facet_grid(factor(gene, levels=c('MMP1718','MMP1511','MMP0419','MMP1227','MMP0335','MMP0111','MMP1362','MMP1303','MMP0209','MMP1611','MMP0033','MMP0166','MMP1361'))~., scales = "free")
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
set12 <- ind.data.mm[grep(top30, ind.data.mm$variant_id),]
set12 <- set12[which(set12$ancestor == "non-ancestral"),]
set12$gene <- sapply(set12$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])
set12$value <- as.numeric(as.character(sub("%", "",set12$freq)))
set12 <- set12[!is.na(set12$value),]
facet.genes2 <- paste("'",gsplot.data2$gene, "'",sep = "", collapse = ",")

## Dotplot for mutation coordinates across lines
pset12 <- ggplot(data=set12, aes(factor(sample)))
pset12 <- pset12 + geom_point(shape=21, color="black", aes(y=position, fill=factor(effect), size=as.numeric(value)), alpha=0.8) + theme_gray()
pset12 <- pset12 + theme(legend.direction = "horizontal",
                         legend.title = element_text(size=6), 
                         legend.key.size = unit(1,"cm"),
                         legend.text = element_text(size = 6))
pset12 <- pset12 + scale_y_continuous(expand = c(0.5,0.5))
pset12 <- pset12 + facet_grid(factor(gene, levels=gsplot.data2$gene)~., scales = "free")
#pset12 <- pset12 + facet_grid(factor(gene, levels=c('MMP0327','MMP0694','MMP1363','MMP0257','MMP0359','MMP0939','MMP1170','MMP0251','MMP1498','MMP1232','MMP0295','MMP0167','MMP1718'))~., scales = "free")
pset12 <- pset12 + scale_fill_manual(values=c("#d62027", "#f9ac61"))
pset12 <- pset12 + guides(fill = guide_legend(title.position = "left",
                                              title.hjust = 0.5, 
                                              label.position = "bottom", 
                                              label.hjust = 0.5))
pset12 <- pset12 + theme(legend.position = "none", axis.text.x = element_text(size = 8)) + labs(x="Evolution Line", y="Mutations")

# SET 1
# Plotting Text values next to dotplot
# prepare data for plotting of text values of name, gscore and # of mutations for individual genes 1st set
set2 <- ind.data.mm
set2 <- set2[!is.na(set2$freq),]
#set2$coord <- sapply(set2$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][3])
set2$gene<- sapply(set2$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])

# get counts of mutations
line.counts <- data.frame()
for(mgene in unique(ind.data.mm$gene)){
  lines <- set2[which(ind.data.mm$gene == mgene),"sample"]
  count <- length(lines)
  line.counts <- rbind(line.counts, cbind(gene=mgene, linecount=count, lines=paste(lines, sep = "", collapse = ":")))
}

## 1set of gene plots
top15.1 = paste("(",paste(gsplot.data1$gene, sep = "", collapse = "|"), ")",sep = "", collapse = "")
pset2.data <- line.counts[grep(top15.1, line.counts$gene),]
pset2 <- ggplot(data=gsplot.data1, aes(x=1))
#pset2 <- pset2 + geom_bar(aes(y=as.numeric(as.character(linecount)), fill=gene),stat = "identity", width = 0.1) + theme_gray()
#gsplot.data1$gene <- gsplot.data1$locus
pset2 <- pset2 + lims(y=c(0,16), x=c(0,3))
pset2 <- pset2 + geom_label(data=gsplot.data1, aes(x=2.5, y=14, label=round(as.numeric(as.character(Gscore)), digits = 1), color="red"),size=3)
pset2 <- pset2 + geom_label(data=gsplot.data1, aes(x=2.5, y=3, label=observed), color="blue", size=3)
pset2 <- pset2 + geom_text(data=gsplot.data1, aes(x=0.5, y=10, label=name), color="black", size=2, hjust="center",vjust="center")
facet.genes1 <- paste("'",gsplot.data1$locus, "'",sep = "", collapse = ",")

# extract legend to plot later
legend.pset2 <- extractLegend(pset2)
pset2 <- pset2 + facet_grid(factor(gsplot.data1$gene, levels=gsplot.data1$locus)~., scales = "free")  + coord_flip()
#pset2 <- pset2 + facet_grid(factor(gene, levels=c('MMP1718','MMP1511','MMP0419','MMP1227','MMP0335','MMP0111','MMP1362','MMP1303','MMP0209','MMP1611','MMP0033','MMP0166','MMP1361'))~.,)  + coord_flip()
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
set21 <- ind.data.mm
set21 <- set21[!is.na(set21$freq),]
#set2$coord <- sapply(set2$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][3])
set21$gene<- sapply(set21$variant_id, function(i) strsplit(as.character(i), split = "-", fixed = T)[[1]][2])

# get counts of mutations
line.counts <- data.frame()
for(mgene in unique(ind.data.mm$gene)){
  lines <- set2[which(ind.data.mm$gene == mgene),"sample"]
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
pset21 <- pset21 + facet_grid(factor(gsplot.data2$gene, levels=gsplot.data2$locus)~.,)  + coord_flip()
#pset21 <- pset21 + facet_grid(factor(gene, levels=c('MMP0327','MMP0694','MMP1363','MMP0257','MMP0359','MMP0939','MMP1170','MMP0251','MMP1498','MMP1232','MMP0295','MMP0167','MMP0643'))~.,)  + coord_flip()
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