# Script to plot intrahost variable sites on protein structures.
# This script is meant to be run from the top level of the Github repository.

# Current directory.
dir="analysis/figures/ProteinStructures"

# Summarize the set of variable sites.
Rscript ${dir}/SummarizeVariableSites.R analysis/figures/LongitudinalFrequencies/LongitudinalVariants.data

# Calculate the distances of each residue in HA and NA
# from the receptor-binding pocket.
python ${dir}/CalculateDistancesHA.py
python ${dir}/CalculateDistancesNA.py

# Permute distances of each residue from the active site.
# Determine the likelihood that the observed distribution of distances
# is due to chance.
Rscript ${dir}/PermuteDistances.R 10000