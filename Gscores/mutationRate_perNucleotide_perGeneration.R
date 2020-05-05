#### Calculate mutation rate per nucleotide per generation for each line
# mu = # of mutations/(number of nuclotides * number of generations) 
setwd("~/Documents/GitHub/evolution-of-syntrophy/Gscores")
source('load_1K_mutation_data.R')
#source('~/Google Drive File Stream/My Drive/R_Scripts/load_mutation_data.R')

mm.genome.size <- 1661137
dv.genome.size <- 3773159

## get mutations for Dv
dv.mutations <- load_1K_mutation_data(include.ancestor = F, organism = "Desulfovibrio")
dv.mutations <- dv.mutations$mutations.1K.na

# dv.early.mutations <- load.mutation.data(org = "dvh",sample_type = "all", remove_missing_EG = T)
# dv.early.mutations <- dv.early.mutations$mutation.data.noancestor
# 
# ## calculate mu for each early generation line for Dv
# emus.dv <- data.frame()
# for(sample in unique(dv.early.mutations$sample)){
#   count <- dim(dv.early.mutations[dv.early.mutations$sample == sample,])[[1]]
#   experiment <- unique(dv.early.mutations[dv.early.mutations$sample == sample, "experiment"])
#   generation <- as.numeric(unique(dv.early.mutations[dv.early.mutations$sample == sample, "transfer"]))*6
#   line <- unique(dv.early.mutations[dv.early.mutations$sample == sample, "line"])
#   mu <- count/(dv.genome.size*generation)
#   emus.dv <- rbind(emus.dv, cbind(sample, line, generation, count, mu))
# }

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




###### Unused Per early gen calculations


# function to get mutation counts and rates for each early gen
get.mutation.counts <- function(org="dvh"){
  
  mutation.counts <- data.frame()
  if(org == "dvh"){
    mutations.file <- dv.early.mutations
  }
  if(org == "mmp"){
    mutations.file <- mmp.early.mutations
  }
  
  for(line in unique(mutations.file$line)){
    mut.15 <-mutations.file[which(mutations.file$line == line & mutations.file$transfer == 15),"variant_id"]
    mut.45 <-mutations.file[which(mutations.file$line == line & mutations.file$transfer == 45),"variant_id"]
    mut.76 <-mutations.file[which(mutations.file$line == line & mutations.file$transfer == 76),"variant_id"]
    mut.118 <-mutations.file[which(mutations.file$line == line & mutations.file$transfer == 118),"variant_id"]
    mut.152 <-mutations.file[which(mutations.file$line == line & mutations.file$transfer == 152),"variant_id"]
    first.union  <- union(mut.15, mut.45)
    second.union <- union(mut.76, first.union)
    third.union <- union(mut.118, second.union)
    total.mutations <- length(third.union)
    
    count.15 <- round(length(mut.15)/(15*6),digits = 2)
    count.45 <- round(length(setdiff(mut.45, mut.15))/(30*6), digits = 2)
    count.76 <- round(length(setdiff(mut.76, first.union))/(31*6),digits = 2)
    count.118 <- round(length(setdiff(mut.118, second.union))/(42*6),digits = 2)
    count.152 <- round(length(setdiff(mut.152, third.union))/(34*6),digits = 2)
    mutation.counts <- rbind(mutation.counts, 
                             cbind(line = line, mut.15 = count.15, mut.45=count.45, mut.76=count.76, mut.118=count.118, mut.152=count.152, total=total.mutations))
    
  }
  return(mutation.counts)
}

# collect mutation counts per generation
dvh.mutation.counts <- get.mutation.counts(org = "dvh")
mmp.mutation.counts <- get.mutation.counts(org = "mmp")


## plotting for mutation counts for dvh
plot.data <- melt(dvh.mutation.counts, measure.vars = c("mut.15", "mut.45", "mut.76", "mut.118", "mut.152"))
plot.data$transfer <- sapply(plot.data$variable, function(i) strsplit(as.character(i), split = "mut.")[[1]][2])
plot.data$org <- "Desulfovibrio"

rateplot <- ggplot(plot.data, aes(x=factor(transfer, levels=c("15","45","76","118","152")), fill=line, color=line, group=line))
#rateplot <- rateplot + geom_bar(aes(y=value), stat="identity",position = "stack")
rateplot <- rateplot + geom_line(aes(y=value))
rateplot <- rateplot + geom_point(aes(y=value))
rateplot <- rateplot + facet_grid(line~.)
rateplot <- rateplot + labs(x="Transfer", y="Mutations/Generation")
rateplot

