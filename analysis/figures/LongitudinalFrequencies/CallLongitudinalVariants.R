# Script takes in BAM summary files for all timepoints from a single patient
# and extracts sites that vary longitudinally, based on given thresholds.
# Meant to be run from top level of Github repository.

library(dplyr)

args<-commandArgs(TRUE)
patient<-args[1]
subtype<-args[2]
minfreq <- as.double(args[3])
mincoverage <- as.integer(args[4])
mintimepoints <- as.integer(args[5])

# Read in data.
Data <- read.table(gzfile(paste("nobackup/SCCA/",patient,"-annotated.summary.gz",
                                sep="")), stringsAsFactors=FALSE, header=FALSE)
colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                    "Sample","Patient","Timepoint","Site","Aliquot","Replicate")

# Filter out samples from low-quality timepoints.
# Exclude haplotypes originating from low-quality timepoints.
Timepoints <- unique((read.table("analysis/figures/ReplicateVariability/HighQualitySamples.data",
                         header=TRUE, stringsAsFactors = FALSE) %>%
                 filter(Patient==patient))$Timepoint)
Data <- Data %>% filter(Timepoint %in% Timepoints)

# Calculate base frequencies at each position.
# Filter out bases that are in intergenic regions.
Data <- Data %>% group_by(Sample, Gene, Pos) %>%
  mutate(Consensus=Base[which.max(Count)], Freq=Count/sum(Count),
         Coverage=sum(Count), Subtype=subtype) %>%
  filter(Gene!="none")

# Add HA codon numbering.
# This means offsetting the codon numbering by 16 for H3 and 13 for H1.
# All other codon numbering remains the same.
Data$Codon <- as.numeric(Data$Codon)
Data <- Data %>% ungroup() %>%
  mutate(Codon=ifelse(Chr=="4-HA",
                        ifelse(Subtype=="H3", Codon-16, Codon-13),Codon))

# Determine the initial consensus base at the first sequenced timepoint.
# Exclude positions from intergenic sequence,
# whose annotations can be misleading for segments in which there are
# multiple genes, i.e. 7-M and 8-NS.
InitialConsensus <- Data %>%
  filter(Timepoint==min(Timepoint), Replicate==1) %>%
  group_by(GenomePos, Gene) %>% filter(Freq==max(Freq)) %>%
  dplyr::select(GenomePos, Chr, Gene, Base, AltAA)
InitialConsensus$Base <- as.character(InitialConsensus$Base)
InitialConsensus$AltAA <- as.character(InitialConsensus$AltAA)
InitialConsensus <- InitialConsensus %>%
  mutate(GenePos=paste(GenomePos,Gene, sep="-")) %>%
  arrange(GenePos) %>%
  filter(Gene!="none")

# Label each genome position in the main dataframe with the initial consensus
# base and amino acid.
# Filter out intergenic positions.
if((nrow(Data)/nrow(InitialConsensus))%%1!=0){
  print("Invalid variant data.")
}
Data <- Data %>%
  mutate(GenePos=paste(GenomePos,Gene,sep="-")) %>%
  arrange(GenePos, Patient, Timepoint, Replicate) %>%
  mutate(Index=floor((row_number()-1)/(nrow(Data)/nrow(InitialConsensus))+1),
         InitBase=InitialConsensus$Base[Index],
         InitAA=InitialConsensus$AltAA[Index])

# Verify that indices are assigned correctly for each gene-position combination.
# These are the indices used to match a position to its initial consensus base.
SortedIndices <- Data %>% dplyr::select(Index,GenePos) %>% distinct() %>%
  arrange(Index) %>% dplyr::select(GenePos)
if(nrow(SortedIndices)!=nrow(InitialConsensus)){
  print("Invalid variant data.")
}
CheckIndices <- function(x){
  if(SortedIndices$GenePos[x]!=InitialConsensus$GenePos[x]){
    return(1)
  }else{
    return(0)
  }
}
if(sum(sapply(seq(1,nrow(SortedIndices)),CheckIndices))!=0){
  print("Invalid variant data.")
}

# Remove unnecessary columns.
Data <- Data %>% 
  dplyr::select(-RefBase,-RefAA,-Syn,-FourfoldSyn,-Subtype,-Index)

# Determine whether each base individually satisfies the variant criteria.
Data <- Data %>%
  mutate(Variant=ifelse(Freq>minfreq & Coverage>mincoverage,1,0))

# Split and then merge the dataframe to account for the two replicates.
DataMerged <- merge(Data %>% filter(Replicate==1),
              Data %>% filter(Replicate==2),
              by=c("Patient","Timepoint","Site","Aliquot","Chr","Pos",
                   "Base","GenomePos","Gene","InitBase","InitAA","AltAA")) %>%
  arrange(Patient, Timepoint, GenomePos, Gene)

# Extract a list of sites that are called as variants in both replicates
# of the minimum required number of timepoints.
VariableSites <- DataMerged %>% 
  filter(Variant.x==1, Variant.y==1, Base!=InitBase) %>%
  group_by(GenomePos, Gene, Base, AltAA) %>% summarize(NumSamples=n()) %>%
  filter(NumSamples>=mintimepoints) %>%
  dplyr::select(GenomePos, Gene, Base, AltAA) %>%
  filter(Gene!="none") %>%
  mutate(GenePos=paste(GenomePos,Gene,sep="-"))

# For each patient, extract the full longitudinal frequency data
# for sites that have been annotated as variable.
Data <- Data %>% filter(Gene!="none", GenePos %in% VariableSites$GenePos) %>%
  dplyr::select(-Variant)

# Determine the derived base and amino acid at each site.
# This is the listing of variable sites.
Data$DerBase <- sapply(Data$GenePos,
                       function(x) (VariableSites$Base)[which(VariableSites$GenePos==x)])
Data$DerAA <- sapply(Data$GenePos,
                     function(x) (VariableSites$AltAA)[which(VariableSites$GenePos==x)])

# Create a short variant name for each site.
Data <- Data %>%
  mutate(Variant=paste(Gene,"-",InitAA,
                       formatC(Codon,width=3,format="d",flag="0"),
                       DerAA,"-",GenomePos,sep=""))

# Export the variant sites and their frequencies through time.
write.table(Data, paste("analysis/figures/LongitudinalFrequencies/",
                        patient,"-variants.data",sep=""),
            row.names=FALSE, quote=FALSE)
