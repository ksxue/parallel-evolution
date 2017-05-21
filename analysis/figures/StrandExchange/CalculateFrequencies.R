# Script takes in haplotype tallies from a given gene in given patient samples.
# It calculates the haplotype frequency at each timepoint.
# Meant to be run from top level of Github repository.

library(dplyr)
library(tidyr)

dir<-"analysis/figures/StrandExchange/"

# Read in haplotype frequency data.
Data <- read.table(paste(dir,"MIX1A.hapsummary", sep=""),
                   header=FALSE, stringsAsFactors = FALSE)
colnames(Data) <- c("Count","Haplotype","Replicate")

# Replace "NA" values, which represent an uncalled base followed by an A,
# with the string "NA" instead of the "not available" indicator.
Data[is.na(Data)] <- "NA"

# Convert nucleotide haplotypes to character states.
# Read in file of sites, ancestral, and derived states.
Sites <- read.table(paste(dir,"MIX1A-4-HA-sites.data", sep=""),
                    header=FALSE, stringsAsFactors = FALSE)
colnames(Sites) <- c("Base","Anc","Der")

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
Data <- Data %>% group_by(Replicate) %>%
  mutate(Freq=Count/sum(Count)) %>%
  dplyr::arrange(Haplotype) %>% 
  dplyr::select(Replicate, Haplotype, Count, Freq) %>%
  ungroup()

# Convert each haplotype to a string for further analysis.
Data$Haplotype <- as.character(Data$Haplotype)


# For each replicate, start at the beginning of each haplotype.
# Proceed along the haplotype, counting how many have experienced
# recombination events. For instance, the haplotype 0011 has
# undergone recombination by the end of the third site.
CalculateRecombination <- function(df, numsites){
  if(numsites==1){
    df <- df %>% mutate(Sum=as.integer(substr(Haplotype,numsites,numsites))) %>%
      mutate_(.dots=setNames("'NR'",
                             paste0("Site",as.character(Sites$Base[numsites]))))
  } else{
    df <- CalculateRecombination(df, numsites-1)
    df <- df %>% mutate(Sum=Sum + as.integer(substr(Haplotype,numsites,numsites))) %>%
      mutate_(.dots=setNames(
        paste0("ifelse(Site",as.character(Sites$Base[numsites-1]),"=='R','R',",
               "ifelse(Sum==",as.character(numsites),",'NR',",
               "ifelse(Sum==0,'NR','R')))"),
        paste0("Site",as.character(Sites$Base[numsites]))))
  }
  return(df)
}
Data <- CalculateRecombination(Data, nrow(Sites))

# Gather the data into tidy form.
DataTidy <- Data %>% gather(Site, State, 
                            -Haplotype, -Replicate, -Count, -Freq, -Sum) %>%
  mutate(Site=as.integer(substr(Site,5,length(Site))))

# Summarize the proportion of recombinant sequences by sequence position
# for each library replicate.
DataTidy <- merge(DataTidy %>% group_by(Replicate, Site) %>%
                    filter(State=="NR") %>% 
                    summarize(FreqRecombinant=1-sum(Freq)),
                  DataTidy %>% group_by(Replicate, Site) %>%
                    filter(State=="R") %>%
                    summarize(MaxRecombHaplotype=max(Freq)),
                  by=c("Replicate","Site"), all=TRUE) %>%
  mutate(MaxRecombHaplotype=ifelse(is.na(MaxRecombHaplotype),0,MaxRecombHaplotype))

# Export summary of recombination frequencies along the length of the haplotype.
write.table(DataTidy, paste0(dir, "CumulativeRecombination.data"),
            row.names=FALSE, quote=FALSE)
