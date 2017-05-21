#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Takes the complete data on intra-host variants
# and summarizes it into a set of sites for plotting on a crystal structure.

library(dplyr)

File <- args[1]

# Read in data.
Data <- read.table(File, header=TRUE, stringsAsFactors = FALSE)

# Identify all sites that are variable within-host.
DataSites <- Data %>% 
  filter(Chr=="4-HA" | Chr=="6-NA") %>%
  mutate(Syn=ifelse(DerAA!=InitAA,"NS","S")) %>%
  group_by(Chr, Codon, Syn) %>%
  summarize(NumPatients=length(unique(Patient))) %>% 
  arrange(Chr, Codon, Syn, NumPatients)

# Export the list of sites.
# Export all nonsynonymous, synonymous, and parallel sites in HA.
write.table(DataSites %>% ungroup() %>%
              filter(Chr=="4-HA", Syn=="NS") %>% dplyr::select(Codon),
            "analysis/figures/ProteinStructures/4-HA-NS.data",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(DataSites %>% ungroup() %>%
              filter(Chr=="4-HA", Syn=="S") %>% dplyr::select(Codon),
            "analysis/figures/ProteinStructures/4-HA-S.data",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(DataSites %>% ungroup() %>%
              filter(Chr=="4-HA", NumPatients>1) %>% dplyr::select(Codon),
            "analysis/figures/ProteinStructures/4-HA-parallel.data",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# Export all nonsynonymous, synonymous, and parallel sites in HA.
write.table(DataSites %>% ungroup() %>%
              filter(Chr=="6-NA", Syn=="NS") %>% dplyr::select(Codon),
            "analysis/figures/ProteinStructures/6-NA-NS.data",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(DataSites %>% ungroup() %>%
              filter(Chr=="6-NA", Syn=="S") %>% dplyr::select(Codon),
            "analysis/figures/ProteinStructures/6-NA-S.data",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(DataSites %>% ungroup() %>%
              filter(Chr=="6-NA", NumPatients>1) %>% dplyr::select(Codon),
            "analysis/figures/ProteinStructures/6-NA-parallel.data",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
