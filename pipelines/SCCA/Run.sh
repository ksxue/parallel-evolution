# Driver script to run SCCA longitudinal analysis pipeline.
# Pipeline trims raw reads, aligns to the influenza reference genome,
# converts to a BAM file, summarizes the base calls at each position,
# and annotates the bases in terms of the amino-acid changes they create.
# Script is designed to be run from the top-level directory of the Github repository.

# Location of pipeline script.
pipeline="pipelines/SCCA/AlignSummarizeAnnotate.sh"
samplesheet="pipelines/SCCA/SCCA-H3N2.samples"

# Run script for all samples.
while read sample reference
do
  qsub -cwd -pe serial 4 \
  -N ${sample} -o nobackup/SCCA/sge/${sample}.o -e nobackup/SCCA/sge/${sample}.e \
  ${pipeline} nobackup/SCCA/${sample}-filtered.1.fastq.gz nobackup/SCCA/${sample}-filtered.2.fastq.gz ${reference}
done < ${samplesheet}