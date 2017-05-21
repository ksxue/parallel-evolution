# Script to analyze haplotypes at several sites of interest in HA in a strand exchange control.
# This script is meant to be run from the top level of the Github repository.

# Current directory.
dir="analysis/figures/StrandExchange"

# Software path for extracting haplotypes.
CountHaplotypes="bin/CountHaplotypes-2.0"

Replicates=( "1" "2" )

for replicate in "${Replicates[@]}"
do

	# Take in arguments.
	f="nobackup/SCCA/MIX1A-NW-"${replicate}".bam" # BAM file of interest
	chr="4-HA" # chromosome of interest
	sites=${dir}/"MIX1A-4-HA-sites.data" # sites of interest
	outdir=${dir} # output directory

	# Run the haplotype inference script.
	${CountHaplotypes} -s <(samtools view ${f}) \
	  -c ${chr} \
	  -i <(cut -f1 ${sites}) \
	  -o ${outdir}/MIX1A-NW-${replicate}-${chr}.haplotypes
	  
	# Summarize the number of reads corresponding to each haplotype.
	# Eliminate paired-end reads that give no information at any of the sites of interest.
	sed $'s/\t//g' ${outdir}/MIX1A-NW-${replicate}-${chr}.haplotypes | \
	 grep -v "N" | sort | uniq -c | \
	 sed -e 's/ *//' -e 's/ /\t/' | sort -n | \
	 sed "s/$/\t${replicate}/" \
	 > ${outdir}/MIX1A-NW-${replicate}-${chr}.hapsummary

done

# Concatenate haplotype summaries for the two replicate libraries.
cat ${outdir}/MIX1A-NW-1-4-HA.hapsummary ${outdir}/MIX1A-NW-2-4-HA.hapsummary \
  > ${outdir}/MIX1A.hapsummary
  
# Remove intermediate files.
rm ${outdir}/MIX1A-NW-1-4-HA.hapsummary
rm ${outdir}/MIX1A-NW-2-4-HA.hapsummary
rm ${outdir}/MIX1A-NW-1-4-HA.haplotypes
rm ${outdir}/MIX1A-NW-2-4-HA.haplotypes

# Calculate the cumulative recombination frequency along the length of the selected sites.
Rscript ${dir}/CalculateFrequencies.R