# Script to identify sites that change longitudinally across the sequenced samples.
# This script is meant to be run from the top level of the Github repository.

# Current directory.
dir="analysis/figures/LongitudinalFrequencies"

# Script to extract and longitudinally varying sites from each patient of interest.
# Script is run as script <patient> <subtype> <minfreq> <mincoverage> <mintimepoints>
CallLongitudinalVariants="analysis/figures/LongitudinalFrequencies/CallLongitudinalVariants.R"

# For each patient of interest, run the script to summarize haplotypes at HA sites of interest.
Patients=( "A" "C" "D" "E" )
for patient in "${Patients[@]}"
do
  Rscript ${CallLongitudinalVariants} ${patient} H3 0.05 200 1
done


# Concatenate the longitudinal variants for each patient.
# Do this such that headers are kept only at the top of the file.
# Remove intermediate files.
head -1 ${dir}/${Patients[0]}-variants.data > ${dir}/LongitudinalVariants.data
tail -n +2 -q ${dir}/*-variants.data >> ${dir}/LongitudinalVariants.data
rm -f ${dir}/*-variants.data

# Take in the longitudinal variants and export the list of sites
# at which variants are called in multiple independent patients.
Rscript ${dir}/IdentifyParallelWithinHostVariants.R
