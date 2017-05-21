# parallel-evolution

This repository contains code and small intermediate data files associated with the manuscript "Parallel evolution of influenza across multiple spatiotemporal scales." All scripts are meant to be run from the top level of the repository.

The directory is organized as follows:

**reference** - This folder contains the A/Brisbane/10/2007 (H3N2) reference sequence in FASTA format, the associated BED annotation file, and crystal structures for HA and NA.

**scripts** - This folder contains custom C++ scripts written to process mapped reads. These scripts tally base frequencies, annotate variants within a gene, and identify and tally reads spanning a haplotype. This code should be compiled before running the sequence analysis pipelines.

**pipelines** - This folder contains shell scripts that process raw sequencing data by filtering out reads that map to the human genome, trimming adapters, mapping reads to a reference sequence, identifying common variants, and annotating them. Raw sequencing data is available in the SRA as BioProject PRJNA364676. See the "Read mapping" section of Materials and Methods for more information.

**analysis** - This folder contains custom R, python, and shell scripts to perform analyses in the manuscript. Each analysis (for instance, an analysis of sequencing depth or HA phylogeny) is given its own folder, which contains the code along with input and output files where these files were small enough to be tracked by Github. Each folder contains a Run.sh script that performs all analyses along with a README file that explains conceptually how the analyses were done. Global influenza sequences were downloaded from the GISAID database and are not included here, per use restrictions. Other small intermediate and output data files are included where possible. See the "Quality filtering," "Variant calling and annotation," "Phylogenetic analysis," "Haplotype inference," "Analysis of global variation," and "Statistical tests of parallelism" sections for more information.