#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Given an alignment summary file in the Year|Position|Base|Count format
# produced by bin/AnalyzeAASiteFrequencies-1.0.py,
# extracts sites at which a non-consensus base ever exceeds a given
# frequency in a given number of years
# and also extracts the bases that exceed those thresholds.
# Also takes in a minimum year threshold, set inclusively.

library(dplyr)
library(tools)

# Verify that exactly two arguments are given.
if(length(args)!=3){
  stop("Arguments: summary file,
       minimum variant frequency, 
       minimum date threshold", 
       call.=FALSE)
}

File <- args[1]
FreqThreshold <- as.numeric(args[2])
DateThreshold <- as.numeric(args[3])

# Import data.
Data <- read.table(File, header=TRUE,
                   stringsAsFactors=FALSE)

# Clean up data formatting.
Data <- Data %>% filter(Year!="0000" & Year!="unkn")
Data$Year <- as.numeric(Data$Year)
Data$Position <- as.integer(Data$Position)

# If the gene being analyzed is HA,
# then convert position numbering from sequential to HA codon numbering.
# Note that positions outputted by the python script are zero-indexed,
# so if the gene is not HA, then offset to one-indexed numbering.
Data <- Data %>%
  mutate(Position=ifelse(Gene=="4-HA",Position-15, Position+1))


# Bin sequences by year and summarize the base counts per year.
Data <- Data %>% mutate(Year=trunc(Year)) %>%
  group_by(Gene, Year, Position, Base) %>%
  summarize(Count=sum(Count))

# Remove ambiguous amino acids.
Data <- Data %>% filter(Base!="X")

# Calculate the frequencies of each amino acid at each position during each year.
Data <- Data %>% group_by(Gene, Year, Position) %>%
  mutate(Freq=Count/sum(Count))

# Rename the M2-exon2 and NEP-exon2 sequences 
# and combine them with data from exon 1 for both genes.
# The resulting data frame should have year-position-aminoacid-frequency data
# for the complete 7-M2 and 8-NEP genes.
# Note that the first exon of 7-M2 includes nine amino acids,
# and the first exon of 8-NEP includes ten amino acids.
Data <- Data %>% ungroup() %>%
  mutate(Gene=ifelse(Gene=="7-M2-exon2","7-M2",Gene),
         Gene=ifelse(Gene=="8-NEP-exon2","8-NEP",Gene),
         Position=ifelse(Gene=="7-M2",Position+10,Position),
         Position=ifelse(Gene=="8-NEP",Position+9,Position))
Data <- rbind(Data,
              Data %>% filter(Gene=="7-M1", Position<=9) %>% mutate(Gene="7-M2"),
              Data %>% filter(Gene=="8-NEP", Position<=10) %>% mutate(Gene="8-NEP")) %>%
  mutate(GeneSite=paste(Gene,"-",
                        formatC(Position,width=3,format="d",flag="0"), sep="")) %>%
  arrange(Gene,Year,Position,Base)

# Identify bases at each site that reach an intermediate frequency at their sites
# in at least one year at or after the date threshold.
DataSites <- Data %>% 
  filter(Freq > FreqThreshold & Freq < 1-FreqThreshold &
           Year >= DateThreshold) %>%
  ungroup() %>% group_by(Gene, Position, Base, GeneSite) %>%
  summarize(NumYears=n(), AvgFreq=mean(Freq))

# Also identify locations in the flu genome where the consensus base
# changes between years.
# These sites are also variable, though it is important to note that
# if the change happens sufficiently quickly, then there may
# note be a detectable variant in the interim.
DataFlipped <- (Data %>% filter(Year >= DateThreshold) %>%
  group_by(Gene, Position, Year) %>% filter(Count==max(Count)) %>%
  ungroup() %>% group_by(Gene, Position) %>%
  filter(length(unique(Base))>1))

# Print number of variable sites.
VariableSites <- unique(c(DataSites$GeneSite, DataFlipped$GeneSite))

# Output global frequency trajectories at all sites.
write.table(Data %>% filter(Year>=DateThreshold) %>%
              dplyr::select(-GeneSite),
            paste(file_path_sans_ext(File),"-frequencies-all.data",sep=""),
            quote=FALSE, row.names=FALSE)

# Output global frequency trajectories at sites called as variable.
write.table(Data %>% filter(GeneSite %in% VariableSites, Year>=DateThreshold) %>%
              dplyr::select(-GeneSite),
            paste(file_path_sans_ext(File),"-frequencies.data",sep=""),
            quote=FALSE, row.names=FALSE)

# Export sites and bases called as variable in the designated timespan.
write.table(Data %>% filter(GeneSite %in% VariableSites) %>% 
              dplyr::select(Gene, Position, Base) %>% distinct(),
            paste(file_path_sans_ext(File),"-sites.data",sep=""),
            quote=FALSE, row.names=FALSE)
