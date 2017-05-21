#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
SequenceDistances=args[1]

library(dplyr)
library(cowplot)
library(tools)

# Read in data of sequence distances.
Data <- read.table(SequenceDistances, 
                   stringsAsFactors = FALSE)
colnames(Data) <- c("SequenceName","Year","NumAADifferences")

# Round the date of each sequence to the proper year, 
# i.e. 2000.9 would be 2000.
Data <- Data %>% mutate(YearShort=floor(Year))

# For each year, calculate the interquartile range of the
# distances from the reference sequence.
# The minimum allowed number of differences is
# the mean minus five times the IQR, and the maximum is
# the mean plus five times the IQR.
# The "allowed" field is 1 if the sequence falls in this
# allowed range and 0 otherwise.
Data <- Data %>% group_by(YearShort) %>%
  mutate(Min=median(NumAADifferences)-5*max(IQR(NumAADifferences),1), 
         Max=median(NumAADifferences)+5*max(IQR(NumAADifferences),1),
         Allowed=ifelse(NumAADifferences<=Max &
                        NumAADifferences>=Min,1,0),
         Allowed=factor(Allowed, levels=c(1,0)))

# Plot the distribution of amino-acid differences per year
p <- ggplot(Data) +
  geom_point(aes(x=Year, y=NumAADifferences, 
                 color=factor(Allowed)), alpha=0.3) +
  scale_color_manual(values=c("black","firebrick3")) +
  guides(color=FALSE)
save_plot(paste0(file_path_sans_ext(SequenceDistances),"-exclusions.pdf"),p)

# Output a list of sequences that are not allowed,
# based on the distance exclusions.
write.table(Data %>% ungroup() %>% filter(Allowed==0) %>% 
                dplyr::select(SequenceName),
            paste0(file_path_sans_ext(SequenceDistances),"-exclusions.data"),
            row.names = FALSE, quote = FALSE, col.names=FALSE)
