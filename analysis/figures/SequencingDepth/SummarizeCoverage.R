library(dplyr)

# Write a function to calculate the mode of a vector.
mode <- function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Write a function to take in the patient summary files across all timepoints
# and calculate coverage in 50-bp bins across the genome.
CalculateCoverage <- function(patient){
  
  # Read in concatenated summary of all samples from given patient.
  Data <- read.table(
    paste("nobackup/SCCA/",patient,"-annotated.summary.gz", sep=""),
    header=FALSE, stringsAsFactors = FALSE)

  colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                      "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                      "Sample","Patient","Timepoint","Site","Aliquot","Replicate")
  
  # Exclude BAL samples from this analysis.
  # Study only patients with H3N2 influenza.
  Data <- Data %>% filter(Site=="NW", Patient %in% c("A","C","D","E"))
  
  # Ensure that genome positions located in multiple genes are not counted
  # multiple times.
  Data <- Data %>% dplyr::select(Sample, Replicate, Chr, Pos, GenomePos, Base, Count) %>%
    distinct()
  
  # Identify unique biological samples, and calculate coverage.
  Data <- Data %>% mutate(BioSample=gsub("..$","",Sample))
  Data <- Data %>% group_by(BioSample, Replicate, Chr, Pos, GenomePos) %>%
    summarize(Coverage=sum(Count))
  
  # Summarize coverage in 50-bp bins along the genome.
  DataSummary <- Data %>% ungroup() %>%
    mutate(BinnedPos=floor(GenomePos/50)*50) %>% 
    group_by(BioSample, Replicate, Chr, BinnedPos) %>%
    summarize(AvgCoverage=mean(Coverage), AvgChr=mode(mode(Chr)))
  
  return(DataSummary)
}


write.table(rbind(CalculateCoverage("A"),
                  CalculateCoverage("C"),
                  CalculateCoverage("D"),
                  CalculateCoverage("E")),
            "analysis/figures/SequencingDepth/SampleCoverage.data", 
            quote=FALSE, row.names=FALSE)
