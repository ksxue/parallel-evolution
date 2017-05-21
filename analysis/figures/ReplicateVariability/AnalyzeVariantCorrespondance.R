# Script calls variants from each patient sample based on frequency and coverage criteria.
# In particular, it analyzes the correspondance in calls between replicates
# for samples that have been designated high-quality.
# Meant to be run from top level of Github repository.

library(dplyr)

args<-commandArgs(TRUE)
minfreq <- as.double(args[1])
mincoverage <- as.integer(args[2])


# Write a function to read a concatenated summary file for a single patient,
# call sites as variable based on a set of frequency and coverage thresholds
# that are used for calling longitudinal variants in the patient data,
# and calculate the percentage of variants called in one replicate
# that are also called in the other replicate.
# Note that this script excludes low-quality samples based on
# earlier criteria.
CallVariants <- function(patient, subtype, minfreq, mincoverage){
  
  # Read in data.
  Data <- read.table(gzfile(paste("nobackup/SCCA/",patient,"-annotated.summary.gz",
                                  sep="")), header=FALSE, stringsAsFactors = FALSE)
  colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                      "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                      "Sample","Patient","Timepoint","Site","Aliquot","Replicate")
  
  # Import function to exclude haplotypes 
  # originating from low-quality timepoints.
  source("analysis/figures/FilterTimepoints.R")
  Timepoints <- read.table("analysis/figures/ReplicateVariability/HighQualitySamples.data",
                           header=TRUE, stringsAsFactors = FALSE)
  Data <- FilterAllowedTimepoints(Data, Timepoints)
  
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
  
  # For each site that is called as a variant in one replicate,
  # determine whether it is also called as a variant in the other replicate.
  VariantCorrespondance <- 
    (Data %>% mutate(VariantCalls=Variant.x+Variant.y) %>%
          filter(VariantCalls!=0) %>% 
       arrange(Patient, Timepoint, Chr, Pos, Base))
  
  return(VariantCorrespondance)
}


write.table(rbind(CallVariants("A","H3",minfreq,mincoverage),
                  CallVariants("C","H3",minfreq,mincoverage),
                  CallVariants("D","H3",minfreq,mincoverage),
                  CallVariants("E","H3",minfreq,mincoverage)), 
            "analysis/figures/ReplicateVariability/VariantCorrespondance.data",
            quote=FALSE, row.names=FALSE)
