========================
ReplicateVariability
========================

This directory compares the two sequencing replicates for each sample, identifies low-replicability samples, and combines this information with information about sequencing depth to generate a list of low-quality samples that are excluded from downstream analyses.

**DATA**

I used BAM summary files that were concatenated across all sequenced samples for each patient. The summary files exclude bases that did not meet a certain quality threshold.

**ANALYSIS**

Using CallVariants.R, I calculated the frequency of variants at each position in the genome for each sequenced sample (this includes the two replicates of each sample). I called positions as variant in each replicate if a non-consensus base (i.e. a base that did not match the original consensus sequence) reached a frequency of at least 0.01, with a coverage at that site of at least 200, in BOTH library replicates. I exported the combined set of sites across all of the patients as "Variants.data".

Using AnalyzeReplicability.R, I calculated a metric of variability based on the difference in replicate frequencies at each timepoint. The metric I calculated was the distance between the point (freq.x, freq.y) from the y=x line, which is, as it turns out, a scaled metric of the difference in frequency between the two points. I set an arbitrary threshold of 0.05 and excluded timepoints for which variability exceeded this threshold. I export this metric for all sequenced samples, and I also export a list of high- and low-replicability timepoints.

Using SummarizeLowQuality.R, I combine the lists of samples that are called as high-quality based on both sequencing coverage (see SequencingDepth analysis) and replicability, and I output a list of high-quality samples that are used for downstream analyses.

Using AnalyzeVariantCorrespondance.R, I identify the correspondance between longitudinal variant calls between the two replicates, using the 5% frequency and 200x coverage thresholds used later in the LongitudinalFrequencies analysis.

**OPEN ISSUES**

None at this time.