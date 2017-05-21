========================
Permutations
========================

**DATA**

*within-host data* I used the set of within-host variants called in LongitudinalFrequencies/LongitudinalVariants.data.

*global frequency data* I used the global frequency data for each site in each flu gene inferred from human H3N2 sequences in GISAID, which is saved as "H3N2-GISAID.data".

**ANALYSIS**

The purpose of this directory is to perform randomization tests to determine how likely it is to see the observed parallelism 1) between patients in our study, and 2) between our patients and the global flu population, purely by chance. This directory also implements a test 3) to determine whether the within-host variants are significantly more likely to be present in antigenic sites.

*1) testing the amount of parallelism between patients in the study*

I used the script ParallelWithinPermutationTest.R to test how likely it would be to see the observed patterns of parallel mutation between patients, purely due to chance. For instance, in HA, there are currently 26 variable sites called among the four patients, and those variants occur at 19 unique sites in the protein (three mutations occur in two patients each, and two mutations occur in three patients each). The script uses previously calculated information from the SequencingDepth analysis about the number of sites in each gene and also varies the proportion of sites in HA that are considered mutable, anywhere from 5 to 100%, since we do not expect that all sites are equally likely to give rise to mutations. The script then performs permutation tests, in which sets of mutations are drawn without replacement from this set of sites for each patient to match the number of mutations observed (i.e. 4 in patient A, 10 in patient E), and it calculates the number of unique sites represented by those mutations. For instance, if there are 26 mutations observed, then the maximum number of unique sites is 26, but any parallel mutation will decrease that number. The script runs a given number of simulations, currently 100000, and from the results, it calculates the proportion of simulations in which mutations fall in fewer than the empirically observed number of sites. I run these simulations for a range of possible proportions of variable sites for each gene and determine the p-value for each proportion of variable sites. (In genes other than HA, NA, and the NS segment there is no parallelism observed, and the script is written so that data for these genes does not show up in the later figures.) 

I also use the script GlobalVariableSites.R and the global flu sequences aligned and analyzed in the GlobalFrequencies directory to calculate the proportion of sites in each protein that at which at least two sequences differ from the consensus identity for all (filtered) human H3N2 sequences between 2000 and 2015. I run a set of permutations using these values and calculate a p-value for this constrained set of mutable sites.

*2) testing the amount of parallelism between patients and the global flu population*

I used the script ParallelScalesPermutationTest.R to test how likely it would be to see the observed overlap in sites of mutation within patients in our study and at the global scale (i.e. overlap within and between hosts). For instance, in HA, there are 19 unique sites of mutation within patients and 79 sites of variation identified between patients, in the global population. Of those 19 unique sites of within-host mutation, 9 of them are also variable in the global population (i.e. there is an overlap of 9 between scales). Like the script above, this script uses information about the number of sites in each gene and varies the proportion of sites that are considered mutable, since we do not expect that all sites are equally likely to give rise to mutations. The script then performs permutation tests: if 79 sites are drawn at random from the length of the protein, then how much overlap do we expect wtih the 19 sites of mutation that are identified within the patients? The script runs a given number of simulations, currently 100000, and from the results, it calculates the proportions of simulations in which there is an equivalent or greater amount of overlap between these two sets of sites than was empirically observed for that particular gene. I run these simulations for a range of possible proportions of variable sites for each gene and determine the p-value for each proportion of variable sites. Note that for analysis, I bin the variants into three categories: HA, NA, and other genes, since there are fewer variants in the other flu genes.

Like above, I also calcualte the proportion of sites in each protein that at which at least two sequences differ from the consensus identity for all (filtered) human H3N2 sequences between 2000 and 2015. I run a set of permutations using these values and calculate a p-value for this constrained set of mutable sites.

**OPEN ISSUES**

None at this time.