========================
GlobalFrequencies
========================

This directory calculates the global frequencies of each amino acid at each site in the flu genome each year for all H3N2 sequences. Note that under the GISAID terms and conditions, these sequences cannot be made publicly available, so this directory is not self-contained.

**DATA**

*within-host data* I used the set of within-host variants called in LongitudinalFrequencies/LongitudinalVariants.data.

*between-host data* I downloaded all human H3N2 sequences for each flu gene segment using GISAID. This was a large number of sequences, particularly for genes like HA.

**ANALYSIS**

The script CallLongitudinalVariants.R (LongitudinalFrequencies directory) takes in the BAM summary file. It identifies all variants that reach a frequency of at least 0.05 in both sequencing replicates of at least one timepoint. It excludes data from low-quality samples, then calculates metrics like coverage and variant frequency at each site in the genome for all samples. It also determines the initial consensus base at each position in the genome at the first sequenced timepoint, and all variants are called relative to this initial consensus. Each site is called separately in each sample and replicate as a variant based on frequency and coverage criteria, and only sites that are called as variants in both sequencing replicates are kept for further analysis.

I used the AlignFilterSummarize.sh pipeline to pairwise align the sequences that I downloaded from GISAID to the H3N2-Brisbane-2007 reference using needle (EMBOSS 6.6.0). I use the bin/CalculateSequenceDistances.py in the bin directory to calculate the amino-acid distance of each sequence to the reference, and then I use the script FilterOutlierSequences.R to exclude all outlier sequences whose distance from the reference is significantly (more than five interquartile ranges, or more than five if the interquartile range is 0, from the median distance of all sequences in that year). I then use the script bin/AnalyzeAASiteFrequencies-2.0.py to summarize the amino-acid frequencies in each year, excluding sequences that contain indels or have been annotated as outliers.

I used the script ExtractGlobalVariableSites.R to extract amino-acid sites that are variable, i.e. some variant is present at a frequency of at least 0.05 in at least two of the years from 2000 to 2016. Because the number of available sequences per year is low until about 2000, I restrict most of my analyses to sites that show variation between 2000 and 2015. These sites are exported in gene-site-base form in "H3N2-GISAID-sites.data", and the allele frequencies each year for those sites are exported in "H3N2-GISAID-frequencies.data".

I then used the script ExtractParallelSites.R to compare the within-host and between-host datasets and export the sites (ParallelSites.data) as well as the specific amino acids (ParallelMuts.data) that are called in both.

**OPEN ISSUES**

Note that under the GISAID terms and conditions, these sequences cannot be made publicly available, so this directory is not self-contained.