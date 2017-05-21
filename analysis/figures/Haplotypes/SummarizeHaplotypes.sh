# Given a BAM file of interest, a list of sites in a chromosome,
# the name of the chromosome, and output directory,
# this script infers the haplotypes at the sites of interest using paired-end reads
# and outputs the complete list of haplotypes, one per line, representing each paired-end read,
# as well as a summary of the haplotypes and their counts.
# This script is meant to be run from the top-level directory of the Github repository.

# Software path for extracting haplotypes.
CountHaplotypes="bin/CountHaplotypes-2.0"

# Take in arguments.
f="$1" # BAM file of interest
chr="$2" # chromosome of interest
sites="$3" # sites of interest
outdir="$4" # output directory

# Parse the sample name.
sample=${f##*/}
sample=${sample%%.*}
patient=${sample:0:1}
timepoint=${sample:1:2}
site=${sample:5:2}
aliquot=${sample:3:1}
replicate=${sample: -1}

# Run the haplotype inference script.
${CountHaplotypes} -s <(samtools view ${f}) \
  -c ${chr} \
  -i <(cut -f1 ${sites}) \
  -o ${outdir}/${sample}-${chr}.haplotypes
  
# Summarize the number of reads corresponding to each haplotype.
# Eliminate paired-end reads that give no information at any of the sites of interest.
sed $'s/\t//g' ${outdir}/${sample}-${chr}.haplotypes |
 grep -vw "N" | sort | uniq -c | \
 sed -e 's/ *//' -e 's/ /\t/' | sort -n | \
 sed "s/$/\t${sample}\t${patient}\t${timepoint}\t${site}\t${aliquot}\t${replicate}/" \
 > ${outdir}/${sample}-${chr}.hapsummary