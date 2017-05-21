library(dplyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
MinCoverage <- as.integer(args[1])
MaxBins <- as.integer(args[2])

# Read in coverage summary data for each sample.
Data <- read.table("analysis/figures/SequencingDepth/SampleCoverage.data",
                   header=TRUE, stringsAsFactors = FALSE)

# Count the number of bins in each sample with coverage below the minimum.
Data <- Data %>% group_by(BioSample) %>%
  filter(AvgCoverage<MinCoverage) %>%
  summarize(NumLowCovBins=n())

# Identify samples for which there are more than the allowed number
# of regions of low coverage.
Data <- Data %>% filter(NumLowCovBins>MaxBins)

# Parse the sample names to recover the patient and timepoint.
Data <- Data %>% mutate(Patient=substr(BioSample,1,1),
                        Timepoint=as.integer(substr(BioSample,2,3))) %>%
  dplyr::select(Patient, Timepoint)

# Output the list of low-coverage samples.
write.table(Data, 
            "analysis/figures/SequencingDepth/LowCoverageSamples.data",
            row.names=FALSE, quote=FALSE)