# Script takes in haplotype tallies from a given gene in given patient samples.
# It calculates the haplotype frequency at each timepoint.
# Meant to be run from top level of Github repository.

library(dplyr)

args<-commandArgs(TRUE)
patient<-args[1]
gene<-args[2]

dir<-"analysis/figures/Haplotypes/"

# Read in haplotype frequency data.
Data <- read.table(paste(dir,patient,"/",patient,"-",gene,"-summaries.data", sep=""),
                   header=FALSE, stringsAsFactors = FALSE)
colnames(Data) <- c("Count","Haplotype",
                    "Sample","Patient","Timepoint","Site","Aliquot","Replicate")

# Replace "NA" values, which represent an uncalled base followed by an A,
# with the string "NA" instead of the "not available" indicator.
Data[is.na(Data)] <- "NA"

# Exclude haplotypes originating from low-quality timepoints.
source("analysis/figures/FilterTimepoints.R")
Data <- FilterAllowedTimepoints(Data,
                                read.table("analysis/figures/ReplicateVariability/HighQualitySamples.data",
                                           header=TRUE, stringsAsFactors = FALSE))

# Exclude haplotypes that are not fully called,
# i.e. that contain N's.
Data <- Data %>% filter(!grepl("N",Haplotype))

# Convert nucleotide haplotypes to character states.
# Read in file of sites, ancestral, and derived states.
Sites <- read.table(paste(dir,patient,"/",patient,"-",gene,"-sites.data", sep=""),
                    header=FALSE, stringsAsFactors = FALSE)
colnames(Sites) <- c("Base","AA","Anc","Der")

# Write a function that converts a base haplotype (i.e. "ATCA")
# to a character haplotype (i.e. "2020")
# using the Sites table above.
HaplotypeBaseToState <- function(haplotype){
  
  # Write a function to convert a string to a vector of characters.
  StringToCharVector <- function(string){
    return(sapply(seq(1:nchar(haplotype)), function(x) substr(haplotype,x,x)))
  }
  
  # Write a function that takes in an index (position in haplotype)
  # and the base at that position
  # and converts that base to 0 if it matches the ancestral state,
  # 1 if it matches the derived state,
  # and 2 otherwise.
  # This uses the Sites dataframe just imported above.
  BaseToState <- function(index, base){
    return(ifelse(base==Sites[index,"Anc"],0,
                  ifelse(base==Sites[index,"Der"],1,2)))
  }
  
  return(paste(mapply(BaseToState, 
                      seq(1:nchar(haplotype)), 
                      StringToCharVector(haplotype)),
         sep="",collapse=""))
}

# Convert base haplotypes to character haplotypes
# using the function defined above.
Data$Haplotype <- sapply(Data$Haplotype, HaplotypeBaseToState)

# Exclude haplotypes that contain states that are neither derived
# nor ancestral.
Data <- Data %>% filter(!grepl("2",Haplotype))

# Calculate the frequency of remaining haplotypes.
Data <- Data %>% group_by(Timepoint, Replicate) %>%
  mutate(Freq=Count/sum(Count)) %>%
  dplyr::arrange(Haplotype) %>% 
  dplyr::select(Timepoint, Replicate, Haplotype, Freq) %>%
  ungroup()

# Create a function to add rows corresponding to haplotypes
# that are not detected at a given timepoint.
# These empty rows, which are akin to pseudocounts,
# are necessary for plotting purposes.
AddEmptyRows <- function(dataframe, timepoint, replicate, haplotype, value){
  if(nrow(dataframe %>% 
          filter(Timepoint==timepoint, Replicate==replicate,
                 Haplotype==haplotype))==0){
    ToAdd <- data.frame(timepoint, replicate, haplotype, value)
    colnames(ToAdd) <- c("Timepoint","Replicate","Haplotype","Freq")
    bind_rows(dataframe, ToAdd)
  }
  else{
    dataframe
  }
}

# Iterate through the dataframe and add rows at timepoints
# where haplotypes are missing.
Timepoints <- sort(unique(Data$Timepoint))
Replicates <- sort(unique(Data$Replicate))
Haplotypes <- sort(unique(Data$Haplotype))
for(i in seq(1, length(Timepoints))){
  for(j in seq(1, length(Replicates))){
    for(k in seq(1, length(Haplotypes))){
      Data <- AddEmptyRows(Data,
                           Timepoints[i],
                           Replicates[j],
                           Haplotypes[k], 0.0001)
    }
  }
}

# Export haplotype frequencies.
write.table(Data %>% arrange(Timepoint, Replicate, Haplotype),
            paste(dir,patient,"/",patient,"-",gene,"-frequencies.data", sep=""),
            quote=FALSE, row.names=FALSE)
