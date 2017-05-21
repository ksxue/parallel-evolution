# Script to analyze haplotypes at several sites of interest in HA in patients A and C.
# This script is meant to be run from the top level of the Github repository.

# Current directory.
dir="analysis/figures/Haplotypes"

# stdio directory.
intdir="nobackup/SCCA/sge"

# Location of waiting script.
waitscript="analysis/figures/Wait.sh"

# Script to extract and summarize haplotypes at sites of interest from a given BAM file.
# Script is run as script <BAMfile> <gene> <sitefile> <outdir>
SummarizeHaplotypes="analysis/figures/Haplotypes/SummarizeHaplotypes.sh"

# For each patient of interest, run the script to summarize haplotypes at HA sites of interest.
Patients=( "A" "C" )
Genes=( "4-HA" )
for patient in "${Patients[@]}"
do
  for gene in "${Genes[@]}"
  do
    for f in nobackup/SCCA/${patient}*.bam
	do
	  sample=${f##*/}
      sample=${sample%%.*}
	  qsub -cwd -N "haplotypes" \
	    -o ${intdir}/${sample}.o -e ${intdir}/${sample} \
	    ${SummarizeHaplotypes} \
        ${f} ${gene} \
        ${dir}/${patient}/${patient}-${gene}-sites.data \
		${dir}/${patient}/
	done
  done
done

# Wait until the previous jobs are complete to continue.
sleep 5
${waitscript} haplotypes 30


# Concatenate the haplotype summaries for each patient and gene.
# Also use the R script CalculateFrequencies.R to calculate the haplotype frequencies.
# CalculateFrequencies.R takes the patient and gene as arguments.
for patient in "${Patients[@]}"
do
  for gene in "${Genes[@]}"
  do
    # Concatenate files.
    cat ${dir}/${patient}/*-${gene}.hapsummary > \
	  ${dir}/${patient}/${patient}-${gene}-summaries.data
	# Calculate haplotype frequencies.
	Rscript ${dir}/CalculateFrequencies.R ${patient} ${gene}
	# Remove intermediate files.
	rm -f ${dir}/${patient}/*-${gene}.hapsummary
	rm -f ${dir}/${patient}/*-${gene}.haplotypes
  done
done

echo "Haplotypes done."
