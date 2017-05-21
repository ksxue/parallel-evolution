#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
YearThreshold <- as.integer(args[1])
NumSequences <- as.integer(args[2])

library(dplyr)

# Given an alignment summary file in the Year|Position|Base|Count format
# produced by bin/AnalyzeAASiteFrequencies-1.0.py,
# extracts all sites that ever show variation.
# Also takes in a minimum year threshold, set inclusively.

# Read in data on global frequencies of mutation in each gene.
Data <- read.table("analysis/figures/GlobalFrequencies/H3N2-GISAID.data",
                   header=TRUE, stringsAsFactors = FALSE)

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

# Remove ambiguous amino acids.
Data <- Data %>% filter(Base!="X")

# Retain only sequences since 2000.
Data <- Data %>% filter(Year>=YearThreshold)

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
              Data %>% filter(Gene=="8-NEP", Position<=10) %>% mutate(Gene="8-NEP"))


# Combine all sequences from all years and determine the count of
# amino acid identities at each position.
Data <- Data %>% group_by(Gene, Position, Base) %>%
  summarize(Count=sum(Count))

# Identify the proportion of variant sites in each gene, 
# at which only a single amino acid identity is ever observed.
# NOTE: take in a parameter that sets the number of sequences required
# for a variant to be called as a variant.
Variants <- Data %>% filter(Count>=NumSequences) %>%
  group_by(Gene, Position) %>%
  summarize(NumIdentities=length(unique(Base))) %>%
  ungroup() %>% group_by(Gene) %>%
  summarize(VariableSites=sum(NumIdentities>1),
            VariableProportion=sum(NumIdentities>1)/n())

# Output proportion of variable sites.
write.table(Variants,
            "analysis/figures/Permutations/GlobalVariableSites.data",
            row.names=FALSE, quote=FALSE)
  