# Script for analyzing deep sequencing data from influenza clinical samples.
# Deep sequencing of longitudinal flu samples, in duplicate.
# Script is designed to be run from the top level of the Github repository.

# Load modules.
module load modules modules-init modules-gs
module load python/2.7.3
module load cutadapt/1.8.3
module load bowtie2/2.2.3
module load samtools/1.3
module load picard/1.43
module load bcftools/1.3.1

# Folder in which to save intermediate files.
dir="nobackup"
projectdir="SCCA"

# Other software paths.
SummarizeBAM="bin/SummarizeBAM-1.21"
CallVariants="pipelines/SCCA/CallVariants.r"
AnnotateVariants="bin/AnnotateVariants-1.1"

# Folder in which to save small output files.
outdir="data"

# Raw sequence reads and reference sequence.
fastq1="$1"
fastq2="$2"
reference="$3"

# Parse the sample name.
sample=${fastq1##*/}
sample=${sample%-*}
patient=${sample:0:1}
timepoint=${sample:1:2}
site=${sample:5:2}
aliquot=${sample:3:1}
replicate=${sample: -1}
echo ${sample}

# Set separate names for the control samples.
if [[ ${sample} =~ "PLASMID" ]]
then
  patient="PLASMID"
fi

if [[ ${sample} =~ "WSN" ]]
then
  patient="WSN"
fi

if [[ ${sample} =~ "PLASMID" ]] || [[ ${sample} =~ "WSN" ]]
then
  timepoint="00"
  site="CL"
  aliquot="A"
fi


# Trim adapter sequences and bases below a quality threshold of 25.
# Also remove all reads that are shorter than 20 bases after trimming.
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
    -q 25 -m 20 \
    -o ${dir}/${projectdir}/${sample}_trimmed-R1.fastq.gz -p ${dir}/${projectdir}/${sample}_trimmed-R2.fastq.gz \
    ${fastq1} ${fastq2} \
    > ${dir}/${projectdir}/${sample}.cutadapt.log


# Align trimmed reads to the appropriate references using Bowtie2.
# Map reads as paired-end reads.
# Use very sensitive settings for end-to-end alignment.
echo "Align reads."
bowtie2 --very-sensitive-local --un-conc-gz ${dir}/${projectdir}/${sample}-unmapped \
	-p 4 \
    -X 2300 \
    -x ${reference%%.*} \
	-k 2 \
    -1 ${dir}/${projectdir}/${sample}_trimmed-R1.fastq.gz \
    -2 ${dir}/${projectdir}/${sample}_trimmed-R2.fastq.gz \
    -S ${dir}/${projectdir}/${sample}.sam \
    2> ${dir}/${projectdir}/${sample}.bt2.log
  
# Use SAMtools to convert SAM files to BAM files.
# Sort the .bam files and also create an index.
echo "Convert SAM files to BAM files."
samtools view -b ${dir}/${projectdir}/${sample}.sam -o ${dir}/${projectdir}/${sample}.bam

# Summarize base frequencies in the sorted BAM file.
echo "Summarize base frequencies."
${SummarizeBAM} -i <(samtools view ${dir}/${projectdir}/${sample}.bam) \
  -f ${reference} -o ${dir}/${projectdir}/${sample}.summary

# Annotate variants as synonymous, nonsynonymous, etc.
echo "Annotate variants."
${AnnotateVariants} -i ${dir}/${projectdir}/${sample}.summary -f ${reference} \
  -b ${reference%%.*}.bed -o ${dir}/${projectdir}/${sample}-annotated.summary
 
 # Annotate the annotation files with the sample name.
sed -i "s/$/\t${sample}\t${patient}\t${timepoint}\t${site}\t${aliquot}\t${replicate}/" \
	${dir}/${projectdir}/${sample}-annotated.summary

# Remove intermediate files after the pipeline finishes.
echo "Remove intermediate files."
rm ${dir}/${projectdir}/${sample}.sam
