========================
LongitudinalFrequencies
========================

This directory calls longitudinal variants for each patient that reach a given frequency threshold relative to the consensus base at the initial sequenced timepoint and outputs the resulting list of variants.

**DATA**

This analysis uses annotated BAM summary files that have been concatenated across all timepoints for each patient.

**ANALYSIS**

The script CallLongitudinalVariants.R takes in the BAM summary file. It identifies all variants that reach a frequency of at least 0.05 in both sequencing replicates of at least one timepoint. It excludes data from low-quality samples, then calculates metrics like coverage and variant frequency at each site in the genome for all samples. It also determines the initial consensus base at each position in the genome at the first sequenced timepoint, and all variants are called relative to this initial consensus. Each site is called separately in each sample and replicate as a variant based on frequency and coverage criteria, and only sites that are called as variants in both sequencing replicates are kept for further analysis.


**OPEN ISSUES**

None currently recorded.