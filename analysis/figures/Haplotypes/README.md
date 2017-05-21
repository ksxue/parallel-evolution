========================
Haplotypes
========================

This directory contains scripts that infer haplotypes at specified sites for HA in patients A and C and summarize haplotype frequency at each sequenced timepoint.

**DATA**

This analysis begins with raw, unsorted BAM files for each sequenced timepoint and replicate. There are also input files "A/A-4-HA-sites.data" and "C/C-4-HA-sites.data" that list the sites of interest (see below for specifications), their amino-acid positions, and the ancestral and derived bases at that position. These sites were chosen manually based on inspection of the allele frequency trajectories for multiple patients.

**ANALYSIS**

*haplotype calling* The script "bin/CountHaplotypes-2.0" takes in an unsorted BAM file, a chromosome name, and an ordered list of one-indexed sites of interest on that chromosome. Note that these sites are by base position, not amino acid position. The script identifies paired-end reads that span the sites of interest and records the bases in each read at the sites of interest. If the read does not cover a site or the coverage is too low, then the script records 'N.' It outputs a .haplotype file with one haplotype per line, with tabs separating the bases recorded at each site. I then use basic bash tools to concatenate these base records into multi-base haplotypes (i.e. "AGTA") and to count how many were observed in each sequenced sample. This information is record in a .hapsummary file. The script Run.sh submits jobs to call haplotypes in all sequenced samples for patients and genes of interest and concatenates the haplotype summaries calculated from each sample. It requires a file like one specified above listing the sites of interest.

*haplotype frequencies and plotting* The R script CalculateFrequencies.R takes in a concatenated -summary.data file listing the counts of each haplotype at each timepoint. It excludes low-quality timepoints, removes incomplete haplotypes, and converts nucleotide haplotypes like "AGTA" to character haplotypes like 0120 using the information about ancestral and derived alleles above. It excludes all haplotypes that include a third allele, none of which are represented at high frequency in the overall population. It outputs a -frequency.data file summarizing the frequency of each haplotype at each timepoint. Crucially for plotting, it also "squares" the haplotype matrix; that is, it adds the equivalent of a pseudocount for haplotypes that are originally absent at any given timepoint. This prevents ggplot2 from plotting gaps in the frequency plot.


**OPEN ISSUES**

None currently recorded.