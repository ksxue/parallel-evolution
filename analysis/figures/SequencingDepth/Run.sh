# Script to calculate sequencing coverage.
# This script is meant to be run from the top level of the Github repository.

dir=analysis/figures/SequencingDepth

# Run R script that calculates coverage in 50bp bins along the genome for each sample
# and then exports the numbers as a file.
Rscript ${dir}/SummarizeCoverage.R

# Run R script that identifies samples with large regions of no coverage.
# Script runs as IdentifyLowCoverage.R <mincoverage> <maxbins>
# It outputs a list of patient-timepoint combinations that have
# more than <maxbins> bins with less than <mincoverage> average coverage.
# <mincoverage> is currently set to 200, the minimum coverage required to call variants.
# <maxbins> is currently set to 16, two bins for each gene, to account for the
# decreased coverage at the very ends of genes caused by Nextera.
Rscript ${dir}/IdentifyLowCoverage.R 200 16

# Calculate the number of base pairs and codons for each gene segment and gene.
Rscript ${dir}/SummarizeGeneLengths.R