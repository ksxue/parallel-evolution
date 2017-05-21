========================
ProteinStructures
========================

**DATA**

I used the set of within-host variant calls from LongitudinalFrequencies/LongitudinalVariants.data. I used the structure 4hmg, of the Aichi/68 HA, and 2bat, of an N2 NA, for the subsequent analysis. Those structures are contained within the "reference" folder from the top-level Github repository.

**ANALYSIS**

Using the script SummarizeVariableSites.R, I take in the within-host variant calls and format them for plotting on the protein structures. I determine the unique sites that are mutated in HA and NA, and I determine the mutations that occur in the patients are synonymous or nonsynonymous, as well as the number of patients that had mutations at that site. There were no sites in the current call set that had synonymous mutations in one patient and nonsynonymous mutations in another. I output the lists of synonymous, nonsynonymous, and parallel sites of mutation for HA and NA.

Using the scripts LabelHA.py and LabelNA.py, I label synonymous, nonsynonymous, and parallel nonsynonymous sites on the protein structure and export a .png image to the main figures directory. Note that these scripts are not run from the "Run" script for the directory and must be run manually from the PyMOL command line, with the top level of the Github repository as the working directory.

Using the scripts CalculateDistancesHA.py and CalculateDistancesNA.py, I use BioPDB to calculate the distance of each residue in the protein from the receptor-binding site (that is, the minimum distance of each residue in the protein from any atom in the sialic acid molecule). I output these distances in 4HMG-SialicAcidDistances.data and 2BAT-SialicAcidDistances.data. Note that this script makes use of the PDB structures stored in the "reference" folder.

Using the script PermuteDistances.R and the list of within-host variant calls, I identify the sites of mutation within patients and determine their distances from the receptor-binding site, based on the list of distances previously calculated. I randomly draw sets of sites from the protein to match the number of variants observed within patients to create a distribution of expected median distances from the receptor-binding site for a random set of sites.

**OPEN ISSUES**

Note that the LabelHA.py and LabelNA.py scripts must be run from the PyMOL command line.