#### Calculate mutation rate per nucleotide per generation for each line
# mu = # of mutations/(number of nuclotides * number of generations) 
setwd("~/Documents/GitHub/evolution-of-syntrophy/Gscores")
source('load_1K_mutation_data.R')

mm.genome.size <- 1661137
dv.genome.size <- 3773159

## get mutations for Dv
dv.mutations <- load_1K_mutation_data(include.ancestor = F, organism = "Desulfovibrio")
dv.mutations <- dv.mutations$mutations.1K.na

## calculate mu for each line for Dv
mus.dv <- data.frame()
for(line in unique(dv.mutations$sample)){
  count <- dim(dv.mutations[dv.mutations$sample == line,])[[1]]
  mu <- count/(dv.genome.size*1000)
  mus.dv <- rbind(mus.dv, cbind(line, count, mu))
}
mus.dv.mean <- mean(as.vector(as.numeric(as.character(mus.dv$mu)))) # 2.813394e-09
mus.dv.sd <- sd(as.vector(as.numeric(as.character(mus.dv$mu)))) #9.588716e-10


## get mutations for Mm
mm.mutations <- load_1K_mutation_data(include.ancestor = F, organism = "Methanococcus")
mm.mutations <- mm.mutations$mutations.1K.na

## calculate mu for each line
mus.dv <- data.frame()
for(line in unique(mm.mutations$sample)){
  count <- dim(mm.mutations[mm.mutations$sample == line,])[[1]]
  mu <- count/(mm.genome.size*1000)
  mus.dv <- rbind(mus.dv, cbind(line, count, mu))
}
mus.mm.mean <- mean(as.vector(as.numeric(as.character(mus.dv$mu)))) # 5.325361e-09
mus.mm.sd <- sd(as.vector(as.numeric(as.character(mus.dv$mu)))) # 9.783198e-10


dv.mutations[dv.mutations$gene_id == "DVU_1295",c("sample","predictor","freq_predictor","freq2")]

mm.mutations[mm.mutations$gene_id == "MMP1718",c("sample", "mutation","effect")]

