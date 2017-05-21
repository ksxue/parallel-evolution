# Script to identify sites that change in the global flu population.
# This script is meant to be run from the top level of the Github repository.

# Current directory.
dir="analysis/figures/GlobalFrequencies"

# Other relevant paths.
refdir="reference"
seqsdir="data/GISAID"
aligneddir="data/GISAID/alignments"
outdir=${dir}/intermediates

####################################################################################
# Align sequences from GISAID to the H3N2-Brisbane-2007 reference
# for each flu gene. Do this for each exon of M2 and NS2 separately.
# The alignment script uses needle to pairwise align each sequence
# to the reference sequence, then excludes sequences that contain indels
# or outliers that have a distance that is too small or large from the reference
# sequence based on other sequences from that year.
####################################################################################
Genes=( "1-PB2" "2-PB1" "3-PA" "4-HA" "5-NP" "6-NA" "7-M1" "7-M2-exon2" "8-NS1" "8-NEP-exon2" )
Chromosomes=( "1-PB2" "2-PB1" "3-PA" "4-HA" "5-NP" "6-NA" "7-M" "7-M" "8-NS" "8-NS" )
NumGenes=${#Genes[@]}

for i in $(seq 0 $((${NumGenes}-1)) )
do
  gene=${Genes[$i]}
  chr=${Chromosomes[$i]}
  
  ${dir}/AlignFilterSummarize.sh \
    ${refdir}/H3N2-Brisbane-2007-${gene}.fasta \
	${seqsdir}/H3N2-${chr}-GISAID.fasta \
	${aligneddir}/H3N2-${gene}-GISAID.aligned \
	${outdir} \
	GISAID \
	H3N2-${gene}-GISAID
done


####################################################################################
# Parse the sequence alignment summaries
# and extract the sites that vary globally.
####################################################################################

# Other relevant paths.
withindata="analysis/figures/LongitudinalFrequencies/LongitudinalVariants.data"

# Script to extract sites that vary in the global flu population.
# Script takes arguments script <summary file> <minfreq> <minyear, inclusive>
ExtractGlobalVariableSites="analysis/figures/GlobalFrequencies/ExtractGlobalVariableSites.R"

# Script to extract sites that evolve in parallel within and between hosts.
# Script takes arguments script <summary file> <withinsites> <betweensites>
ExtractParallelSites="analysis/figures/GlobalFrequencies/ExtractParallelSites.R"

# Concatenate multiple-alignment summary files for each gene.
# Note that I am currently using summary files generated in a previous analysis,
# as opposed to generating them as part of this script.
for gene in "${Genes[@]}"
do
  if [[ ${gene} == ${Genes[0]} ]]
  then
    head -1 ${outdir}/H3N2-${gene}-GISAID-frequencies.data | \
      sed "s/$/\tGene/" \
      > ${dir}/H3N2-GISAID.data
  fi
  tail -n +2 -q ${outdir}/H3N2-${gene}-GISAID-frequencies.data | \
    sed "s/$/\t${gene}/" >> ${dir}/H3N2-GISAID.data
done

# Extract global variable sites.
Rscript ${ExtractGlobalVariableSites} ${dir}/H3N2-GISAID.data 0.1 2000

# Extract the set of parallel sites,
# i.e. those that are shared between the within-host and between-host data.
Rscript ${ExtractParallelSites} ${withindata} ${dir}/H3N2-GISAID-sites.data