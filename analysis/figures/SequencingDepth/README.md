========================
SequencingDepth
========================

This directory analyzes sequencing depth in all samples and identifies those with low coverage.

**DATA**

I used BAM summary files that were concatenated across all sequenced samples for each patient. The summary files exclude bases that did not meet a certain quality threshold.

**ANALYSIS**

Using SummarizeCoverage.R, I calculated coverage in 50bp bins for each gene and outputted the data for all patients, samples, and genes in the file SampleCoverage.data. I use the script IdentifyLowCoverage.R to identify samples that have more than 16 bins with an average coverage below 200x (to simulate the two ends of the eight gene segments, which are expected to have low coverage because we performed tagmentation using Nextera). I also use the script SummarizeGeneLengths.R to take in one of the BAM summary files and determine the number of nucleotides and codons in each chromosome and gene in the reference genome.

**OPEN ISSUES**

None at this time.