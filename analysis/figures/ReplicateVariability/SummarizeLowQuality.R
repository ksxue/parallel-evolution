library(dplyr)

# Read in sets of patients and timepoints called as low-quality
# based on either coverage or variant replicability.
LowQ <- rbind(
  read.table("analysis/figures/ReplicateVariability/LowReplicabilitySamples.data",
             header=TRUE, stringsAsFactors = FALSE),
  read.table("analysis/figures/SequencingDepth/LowCoverageSamples.data",
             header=TRUE, stringsAsFactors = FALSE))

# Combine these two datasets and export a single set of low-quality timepoints.
write.table(LowQ %>% distinct() %>% arrange(Patient,Timepoint),
            "analysis/figures/ReplicateVariability/LowQualitySamples.data",
            row.names=FALSE, quote=FALSE)

# Read in the list of all patients and timepoints
# and output a list of high-quality samples.
Samples <- read.table("analysis/figures/PatientTimelines/Samples.txt",
                      header=TRUE, stringsAsFactors = FALSE) %>%
  mutate(PatientTimepoint=paste(Patient,Timepoint,sep="-"))

LowQ <- LowQ %>%
  mutate(PatientTimepoint=paste(Patient,Timepoint,sep="-"))

write.table(Samples %>% 
              filter(!(PatientTimepoint %in% LowQ$PatientTimepoint)) %>%
              dplyr::select(Patient, Timepoint),
            "analysis/figures/ReplicateVariability/HighQualitySamples.data",
            row.names=FALSE, quote=FALSE)