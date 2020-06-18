############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Simulates similar number of mutations in Desulfovibrio genome to calculate G-scores
###############################################################################

library('progress')
# load dvh genome features file
dvh.genome <- read.delim("~/Documents/GitHub/evolution-of-syntrophy/data/GCF_000195755.1_ASM19575v1_feature_table.txt",
                         header = T, sep="\t", stringsAsFactors=F)
# collect only genes
dvh.genes <- dvh.genome[which(dvh.genome$X..feature == "gene"),]

# create position list for each gene
coordinate.list <- list()
all.coordinates <- vector()
for(gene in dvh.genes$locus_tag){
  start <- dvh.genes[which(dvh.genes$locus_tag == gene),"start"]
  end <- dvh.genes[which(dvh.genes$locus_tag == gene),"end"]
  initial.pos <- min(start, end)
  final.pos <-max(start, end)
  # if plasmid, add coordinates to the end of the lst chromosome position
  if(dvh.genes[which(dvh.genes$locus_tag == gene),"seq_type"] == "chromosome"){
    interval <- seq(from=initial.pos, to=final.pos, by = 1)
  } else{
    interval <- seq(from= (3570714 + initial.pos), to= (3570714 + final.pos), by = 1)
  }
  coordinate.list[[gene]] <- interval
  all.coordinates <- c(all.coordinates, interval)
}

### randomly select coordinates from mutated genome
# total # of mutations
total.genes <- length(dvh.genes$GeneID)
Ntot <- 166
# total length of coding genome
Ltot <- sum(dvh.genes$feature_interval_length)
gscore.table <- data.frame()
gscores = vector()
simulation.target = 1000
pb <- progress_bar$new(
  format = " :what Remaining [:bar] :percent eta: :eta (:spin)",
  total = simulation.target, clear = FALSE, width= 100)

#pb <- txtProgressBar(min = 1, max = simulation.target, style = 3 ) # progress bar
#cat("Completed %:", "\n")
for(simulation.no in 1:simulation.target){
  pb$tick(tokens = list(what = paste("Sim # ", simulation.no, sep="")))
  Sys.sleep(1 / simulation.target)
  #cat("Running simulation # ", simulation.no, "\n")
  #progress <- setTxtProgressBar(pb, simulation.no)
  #cat(paste("Running simulation # ", simulation.no, progress, "\n", sep = " "))

  # Random position selection
  random.positions <- sample(all.coordinates, Ntot)
  hit.genes <- vector()
  for(position in random.positions){
    hit.indice <- which(sapply(1:length(coordinate.list), function(x) any(coordinate.list[[x]] == position)))
    hit.gene <- names(coordinate.list[hit.indice])
    #hit.gene <- names(coordinate.list[grep(position, coordinate.list)])
    hit.genes <- append(hit.genes, hit.gene)
  }
  
  # calculate Gscore for each gene
  #cat("Calculating Gscores for all genes in the genome..\n")
  gscore.sum = 0
  gene.no = 0
  
  for(gene in dvh.genes$locus_tag){
    Li <- dvh.genes[which(dvh.genes$locus_tag == gene),"feature_interval_length"]
    Ni <- length(grep(gene, hit.genes))
    # calculate expected number of mutations
    Ei <- Ntot*(Li/Ltot)
    # calculate GScore
    Gscore <- 2*Ni*log1p(Ni/Ei)
    
    #cat(gene, Li, Ni, Ei, Gscore, "\n")
    gscore.sum = gscore.sum + Gscore
    #gscore.table <- rbind(gscore.table, 
                          #cbind(iteration= simulation.no, locus=gene, Length=Li, observed=Ni, expected=Ei, Gscore=Gscore))
  }
  #cat("Appending Gscore from this iteration.. \n")
  gscores <- append(gscores, gscore.sum)
}

write.table(gscores, file="~/Google Drive File Stream/My Drive/Manuscripts/Syntrophy-SingleCell/Data/simulated_Gscores_DvH.txt", sep="\t")
mean.gscores <- mean(gscores) # 980.2369
sd.gscores <- sd(gscores) #18.10171
observed.gscores <- read.delim("~/Google Drive File Stream/My Drive/Manuscripts/Syntrophy-SingleCell/Data/Gscore_high_moderate_dvh.txt", sep="\t", header=T)
sum.observed <- sum(observed.gscores[,"Gscore"]) #1092.617
zscore <- (sum.observed - mean.gscores)/sd.gscores #6.208239






