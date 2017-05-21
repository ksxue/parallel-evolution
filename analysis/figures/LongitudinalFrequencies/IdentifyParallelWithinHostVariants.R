library(dplyr)
# Import longitudinal variants.
Data <- read.table(
  "analysis/figures/LongitudinalFrequencies/LongitudinalVariants.data",
  header=TRUE, stringsAsFactors = FALSE)

# Identify variants that are present in more than one patient.
write.table(Data %>% group_by(Patient, Gene, Codon) %>%
  distinct() %>% ungroup() %>% group_by(Gene, Codon) %>%
  summarize(NumPatients=n()) %>% filter(NumPatients > 1) %>%
  dplyr::select(Gene, Codon),
  "analysis/figures/LongitudinalFrequencies/ParallelWithin.data", 
  row.names=FALSE, quote=FALSE)