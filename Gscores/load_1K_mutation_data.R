## This function loads 1000 generation mutation data for DvH or Mmp
# Use include.ancestor=F if you dont want ancestral mutations

load_1K_mutation_data <-function(include.ancestor=F, organism="Desulfovibrio"){
  if(organism == "Desulfovibrio"){
    mut.data <- read.delim("/Volumes/omics4tb/sturkarslan/syntrophy_raw_sequences/dvh_mutations_allsamples_attributes_5152019.txt", header=T, sep="\t",)
  }
  if(organism == "Methanococcus"){
    mut.data <- read.delim("/Volumes/omics4tb/sturkarslan/syntrophy_raw_sequences/mmp_mutations_allsamples_attributes_5152019.txt", header=T, sep="\t",)
  }
  # get 1000-gen data
  mutations.1000 <- mut.data[which(mut.data$experiment == "1000-gen"),]
  # get ancestor mutations
  ancestor.mutations <- mut.data[which(mut.data$sample == "AN_Coculture-Ancestor"),]
  ancestor.mutated.genes <- unique(ancestor.mutations$gene_id)
  ancestor.mutation.ids <- unique(ancestor.mutations$variant_id)
  # Filter mutation data for non-ancestral mutations
  allmutation.ids <- unique(mutations.1000$variant_id)
  nonancestor.mutation.ids <- setdiff(allmutation.ids, ancestor.mutation.ids)
  # only get non-ancestral mutations
  mutations.1000.na <- data.frame()
  for(mutation.id in nonancestor.mutation.ids){
    mutations.1000.na <- rbind(mutations.1000.na, mutations.1000[which(mutations.1000$variant_id == mutation.id),])
  }
  return(list(ancestor.mutated.genes = ancestor.mutated.genes,
                ancestor.mutation.ids = ancestor.mutation.ids,
                mutations.1K.na = mutations.1000.na,
              all.mutations = mut.data))
}
