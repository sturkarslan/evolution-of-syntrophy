############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Function to load all the mutation Data,  modify gene names and other 
# modifications to prepare for plotting
###############################################################################

load.mutation.data <- function(org="dvh", sample_type=c("1000-gen", "all", "Early-gen"), remove_missing_EG=FALSE){
  # load early generation mutations data
  if(org == "dvh"){
    file <- "~/Documents/GitHub/evolution-of-syntrophy/data/dvh_mutations_allsamples_attributes_3112020.txt"
  }
  if(org == "mmp"){
    file <- "~/Documents/GitHub/evolution-of-syntrophy/data/mmp_mutations_allsamples_attributes_3112020.txt"
  }
  
  ## load mutation data
  mutation.data <- read.delim(file, header=T, sep="\t")
  
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
  
  # filter for different samples
  if(sample_type=="all"){
    cat("Getting all the data without any filtering\n" )
    # Only include early geenration, 1K and Ancestor data
    mutation.data.filtered <- mutation.data[mutation.data$sample != "EG_HA3_118.New" & (mutation.data$experiment == "1000-gen" | mutation.data$experiment == "Early-gen" | mutation.data$experiment == "EPD" | mutation.data$experiment == "Clonal-isolates" | mutation.data$sample == "AN_Coculture-Ancestor"),]
  } 
  if(sample_type == "1000-gen"){
    cat("Getting only 1000-gen data\n" )
    # Only include 1000-gen
    mutation.data.filtered <- mutation.data[mutation.data$sample != "EG_HA3_118.New" & (mutation.data$experiment == "1000-gen"),]
  }
  if(sample_type=="Early-gen"){
    cat("Getting only Early-gen data\n" )
    # Only include early generations
    mutation.data.filtered <- mutation.data[mutation.data$sample != "EG_HA3_118.New" & (mutation.data$experiment == "Early-gen"),]
  }
  if(sample_type=="ancestor"){
    cat("Getting only ancestral data\n" )
    # Only include ancestral mutations
    mutation.data.filtered <- mutation.data[mutation.data$sample != "EG_HA3_118.New" & (mutation.data$sample == "AN_Coculture-Ancestor"),]
  } 
  
  # Format transfer and line names
  mutation.data.filtered$line <- sapply(mutation.data.filtered$sample, function(i) strsplit(as.character(i), split = "_")[[1]][2])
  mutation.data.filtered$transfer <- sapply(mutation.data.filtered$sample, function(i) strsplit(as.character(i), split = "_")[[1]][3])
  mutation.data.filtered$transfer <- sub("118.New", 119, mutation.data.filtered$transfer)
  mutation.data.filtered$transfer <- sub("118.Early", 118, mutation.data.filtered$transfer)
  mutation.data.filtered[which(mutation.data.filtered$line == "Coculture-Ancestor"),"transfer"] <- 0
  mutation.data.filtered[which(mutation.data.filtered$experiment == "1000-gen"),"transfer"] <- 152
  
  if(remove_missing_EG == TRUE){
    # remove 1000-gen lines without early gen data
    mutation.data.filtered <- mutation.data.filtered[!(mutation.data.filtered$line %in% c("HE2", "HR1", "UA2", "UE2")),]
  }
  
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
  cat(paste(sample_type, " variants: ", sep = ""), length(unique.nonancestral.variants), paste(sample_type, " locus: ", sep = ""), length(unique.nonancestral.locus), "\n" )
  return(list(mutation.data=mutation.data.filtered, ancestral.variants=unique.ancestral.variants, mutation.data.noancestor=mutation.data.noancestor))
}
