# Script calls variants from each patient sample based on frequency and coverage criteria.
# Meant to be run from top level of Github repository.

library(dplyr)

# Write a function to read a concatenated summary file for a single patient,
# call sites as variable based on a set of frequency and coverage thresholds,
# and calculate the correlation between replicates.
# Note that the same sites may be plotted from different timepoints.
CallVariants <- function(patient, subtype, minfreq, mincoverage){
  
  
  # Read in data.
  Data <- read.table(gzfile(paste("nobackup/SCCA/",patient,"-annotated.summary.gz",
                                  sep="")), header=FALSE, stringsAsFactors = FALSE)
  colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                      "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                      "Sample","Patient","Timepoint","Site","Aliquot","Replicate")
  
  # Calculate base frequencies at each position.
  Data <- Data %>% group_by(Sample, Chr, Pos) %>%
    mutate(Consensus=Base[which.max(Count)], Freq=Count/sum(Count),
           Coverage=sum(Count), Subtype=subtype)
  
  # Add HA codon numbering as well as a text label indicating gene and codon number.
  Data$Codon <- as.numeric(Data$Codon)
  Data <- Data %>% ungroup() %>%
    mutate(CodonHA=ifelse(Chr=="4-HA",
                          ifelse(Subtype=="H3", Codon-16, Codon-13),Codon))
  
  # Ensure that positions in multiple genes are not counted multiple times.
  Data <- Data %>% 
    dplyr::select(Patient, Timepoint, Site, Aliquot, Chr, Pos, Base, GenomePos, Replicate,
                  Count, Coverage, Freq) %>%
    distinct()
  
  # Determine the initial consensus at the first sequenced timepoint.
  InitialConsensus <- Data %>% filter(Timepoint==min(Timepoint), Replicate==1) %>%
    group_by(GenomePos) %>% filter(Freq==max(Freq)) %>%
    dplyr::select(GenomePos, Chr, Base)
  InitialConsensus$Base <- as.character(InitialConsensus$Base)
  
  # Label each genome position in the main dataframe with the initial consensus base.
  Data <- Data %>%
    mutate(InitBase=InitialConsensus$Base[GenomePos])
  
  # Call variants that reach a certain frequency and coverage threshold.
  Data <- Data %>%
    mutate(Variant=ifelse((Base!=InitBase & Freq>minfreq & Coverage>mincoverage), 1, 0))
  
  # Split and then merge information for the two replicates.
  Data <- Data %>%
    dplyr::select(Patient, Timepoint, Site, Aliquot, Chr, Pos, Base, GenomePos, Replicate,
                  Count, Coverage, Freq, InitBase, Variant)
  Data <- merge(Data %>% filter(Replicate==1), Data %>% filter(Replicate==2),
                by=c("Patient","Timepoint","Site","Aliquot","Chr","Pos",
                     "Base","GenomePos","InitBase"))
  
  # Extract all sites that are called as variants in BOTH replicates.
  Data <- Data %>% filter(Variant.x==1 & Variant.y==1)
    
  return(Data)
}


write.table(rbind(CallVariants("A","H3",0.01,200),
                  CallVariants("C","H3",0.01,200),
                  CallVariants("D","H3",0.01,200),
                  CallVariants("E","H3",0.01,200)), 
            "analysis/figures/ReplicateVariability/Variants.data",
            quote=FALSE, row.names=FALSE)