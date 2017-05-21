library(dplyr)

# Write a function to filter a dataframe containing data from many patients and timepoints
# and exclude all data deriving from timepoints that are NOT on the allowed list.
# Both the data and timepoints arguments must contain 
# columns labeled "Patient" and "Timepoint".
FilterAllowedTimepoints <- function(data, timepoints){
  
  data <- data %>% mutate(PatientTimepoint=paste0(Patient,Timepoint))
  timepoints <- timepoints %>% mutate(PatientTimepoint=paste0(Patient,Timepoint))
  
  data <- data %>% filter(PatientTimepoint %in% timepoints$PatientTimepoint)
  
  return(data %>% dplyr::select(-PatientTimepoint))
}
