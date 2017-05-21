#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Given a set of intra-host and inter-host sites,
# extracts the sites that are shared between the two call sets.

library(dplyr)
library(tools)

FileWithin <- args[1]
FileBetween <- args[2]

Within <- read.table(FileWithin, header=TRUE, stringsAsFactors = FALSE)
Between <- read.table(FileBetween, header=TRUE, stringsAsFactors = FALSE)

# Extract the set of sites and set of mutations that occur globally.
BetweenSites <- unique((Between %>% 
  mutate(GeneSite=paste(Gene,"-",
                        formatC(Position,width=3,format="d",flag="0"), 
                        sep="")))$GeneSite)
BetweenMutations <- unique((Between %>% 
  mutate(GeneMut=paste(Gene,"-",
                        formatC(Position,width=3,format="d",flag="0"), 
                        "-",Base,sep="")))$GeneMut)

# Filter the within-host data for variants that occur between hosts.
ParallelSites <- Within %>% 
  mutate(GeneSite=paste(Gene,"-",
                        formatC(Codon,width=3,format="d",flag="0"), 
                        sep="")) %>%
  filter(GeneSite %in% BetweenSites) %>%
  dplyr::select(Gene, Codon) %>% distinct() %>%
  arrange(Gene, Codon)
ParallelMuts <- Within %>% 
  mutate(GeneMut=paste(Gene,"-",
                        formatC(Codon,width=3,format="d",flag="0"), 
                        "-",DerAA,sep="")) %>%
  filter(GeneMut %in% BetweenMutations) %>%
  dplyr::select(Gene, Codon, DerAA) %>% distinct() %>%
  arrange(Gene, Codon, DerAA)

# Export parallel sites.
write.table(ParallelSites, 
            "analysis/figures/GlobalFrequencies/ParallelSites.data",
            quote=FALSE, row.names=FALSE)

write.table(ParallelMuts, 
            "analysis/figures/GlobalFrequencies/ParallelMuts.data",
            quote=FALSE, row.names=FALSE)
