========================
HAPhylogeny
========================

This directory generates a phylogeny to compare the patient sequences at the initial consensus timepoints to the other flu sequences from the USA around the same time period. Note that due to the GISAID terms of use, these other sequences cannot be made publicly available, so this directory is not self-contained.

**DATA**

I downloaded sequences from GISAID on November 30, 2016, with the following search parameters:

type: A
subtype: H3N2
host: human
location: North America - United States
collection date: 2004-01-01 to 2007-12-31
required segments: HA, only complete

I downloaded the 503 resulting isolates as FASTA DNA sequences with the following FASTA header:
Isolate name | DNA Accession no. | Type | Lineage | Collection date | Segment | Passage details/history 

The date format was "Year fraction". I saved the file as "H3N2-GISAID-USA-2004-2007-HA.fasta". Because sequences from GISAID are downloaded with Windows line endings, I also changed these to Unix line endings. I also downloaded the GISAID acknowledgement table as "GISAID-Acknowledgement-Table."

As reference points, I also included HA sequences for the A/Wisconsin/67/2005 and A/Brisbane/10/2007 reference strains. The Brisbane/2007 reference is the reference sequence to which all other sequences in this study were aligned. The Wisconsin/2005 sequence was Genbank accession KM821341.1, one of the many HA sequences listed under that strain name.

Using the custom script SummarizeBAM, I also generated FASTA consensus sequence files for the four patients based on the first sequenced timepoint. The HA sequences are contained in "H3N2-SCCA-HA.fasta"

**ANALYSIS**

Using ConstructPhylogeny.py, I pairwise align the sequences to the H3N2-Brisbane-2007 HA reference sequences using needle. I parse the sequence metadata to retain only GISAID sequences with "unpassaged", "p0", or "original" as their passage histories, and then I construct a phylogeny using RAxML. I manually use TempEstv1.5 to estimate the root of the best tree (RAxML_bestTree...) to form the phylogeny that I plot (RAxML_bestTree...rooted).

**OPEN ISSUES**

The sequences can not be made available here under the GISAID terms of use, so this directory is not self-contained.