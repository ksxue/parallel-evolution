library(dplyr)

# Write a function to calculate the length of each influenza gene
# in base pairs and codons.
# This helps with some plots and statistical tests.
# And it is useful to do only once.

# Read in the summary file for an arbitrary patient.
Data <- read.table(
  paste("nobackup/SCCA/A-annotated.summary.gz", sep=""),
  header=FALSE, stringsAsFactors = FALSE)

colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                    "Sample","Patient","Timepoint","Site","Aliquot","Replicate")

# For each gene segment, count the total number of base pairs.
# Note that this depends slightly on the particular reference genome
# and how much flanking region was included.
write.table(Data %>% group_by(Chr) %>% summarize(NumBases=max(Pos),
                                                 Midpoint=as.integer(mean(GenomePos)),
                                                 GenomeMin=min(GenomePos),
                                                 GenomeMax=max(GenomePos)),
            "analysis/figures/SequencingDepth/GenomeSegmentSummary.data",
            quote=FALSE, row.names=FALSE)

# For each gene, count the total number of base pairs in the gene
# and the number of codons.
write.table(Data %>% filter(Gene!="none") %>%
              group_by(Gene) %>% 
              summarize(NumBases=length(unique(Pos)), NumCodons=length(unique(Codon))),
            
            "analysis/figures/SequencingDepth/GenomeGeneSummary.data",
            quote=FALSE, row.names=FALSE)
