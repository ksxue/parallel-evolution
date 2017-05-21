# Driver script to run SCCA longitudinal analysis pipeline.
# Script is designed to be run from the top-level directory of the Github repository.

# Location of pipeline script.
pipeline="pipelines/SCCA/FilterOutHumanReads.sh"
samplesheet="pipelines/SCCA/SCCA-H3N2.samples"
reference="/net/shendure/vol10/nobackup/shared/alignments/bowtie2-2.0.2/human_g1k_hs37d5/hs37d5"

# Run script for all samples.
while read sample ref
do
  qsub -cwd -pe serial 4 \
  -N ${sample} -o nobackup/SCCA/sge/${sample}.o -e nobackup/SCCA/sge/${sample}.e \
  ${pipeline} raw/SCCA/${sample}_R1.fastq.gz raw/SCCA/${sample}_R2.fastq.gz ${reference}
done < ${samplesheet}