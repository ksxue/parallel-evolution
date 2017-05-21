# Script analyzes variability between replicates 
# and sets exclusion criteria for low-quality timepoints.
# Meant to be run from top level of Github repository.

library(dplyr)
library(ggplot2)
library(cowplot)

# Read in data.
Data <- read.table("analysis/figures/ReplicateVariability/Variants.data", 
                   header=TRUE, stringsAsFactors = FALSE)

# Calculate the difference in allele frequency between the two replicates.
Data <- Data %>% mutate(XYDistance=abs(Freq.y-Freq.x))

# Write a function to calculate variability between variant frequencies at each timepoint
# in each patient.
CalculateVariability <- function(patient){
  
   DataPatient <- Data %>% filter(Patient==patient)
   DataPatientSummary <- DataPatient %>% 
     group_by(Patient, Timepoint) %>% 
     summarize(AvgXYDistance=mean(XYDistance),
               MaxXYDistance=max(XYDistance))
   
   return(DataPatientSummary)
}

DataSummary <- rbind(CalculateVariability("A"),
      CalculateVariability("C"),
      CalculateVariability("D"),
      CalculateVariability("E"))

# Export the list of variability metrics for each timepoint.
write.table(DataSummary, 
            "analysis/figures/ReplicateVariability/VariabilityMetrics.data",
            quote=FALSE, row.names=FALSE)

# Export the list of timepoints for each patient
# for which the metric of variability lies ABOVE a certain threshold.
# That is, the samples are below a certain quality threshold.
write.table(DataSummary %>% filter(AvgXYDistance>=0.05) %>%
              dplyr::select(Patient, Timepoint), 
            "analysis/figures/ReplicateVariability/LowReplicabilitySamples.data",
            quote=FALSE, row.names=FALSE)