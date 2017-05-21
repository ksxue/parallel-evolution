# Script to identify sites that change longitudinally across the sequenced samples.
# This script is meant to be run from the top level of the Github repository.

# Current directory.
dir="analysis/figures/Permutations"

# Use the global frequency data in the GlobalFrequencies folder
# to identify sites that do not vary from 2000 to the present.
Rscript ${dir}/GlobalVariableSites.R 2000 2

# Perform a permutation test to determine how much parallelism is expected
# given the number of variable sites observed for each patient
# if each site is equally likely to be variable.
# The parameter gives the number of permutations run for each gene.
Rscript ${dir}/ParallelWithinPermutationTest.R 100000

# Perform a permutation test to determine how much parallelism is expected
# given the number of variable sites observed for each patient
# and the number of variable sites in the global population,
# for HA, NA, and the other flu genes.
# The parameter gives the number of permutations run for each gene.
Rscript ${dir}/ParallelScalesPermutationTest.R 10000