# Script to calculate sequencing coverage.
# This script is meant to be run from the top level of the Github repository.

dir=analysis/figures/ReplicateVariability

# Run R script that calls variants in all sequenced samples
# based on frequency and coverage.
# Exports variants that meet the criteria in both library replicates.
Rscript ${dir}/CallVariants.R

# Run R script that analyzes the variability between replicates
# for each timepoint.
# Variability metric is difference between replicate frequencies.
# Exports variability metrics for each sequenced timepoint.
# Exports list of timepoints for each patient for which the variability
# between replicates is below the specified threshold.
# The current threshold requires the average variant replicability
# to be below 0.1.
Rscript ${dir}/AnalyzeReplicability.R

# Run R script that combines the set of timepoints called as low-quality
# based on either low coverage or low variant replicability.
Rscript ${dir}/SummarizeLowQuality.R

# Analyze correspondance between variant calls at the
# frequency and coverage thresholds used to call
# longitudinal variants.
Rscript ${dir}/AnalyzeVariantCorrespondance.R 0.05 200