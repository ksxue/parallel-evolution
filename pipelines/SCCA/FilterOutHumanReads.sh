# Script for analyzing deep sequencing data from influenza clinical samples.
# Designed to analyze data from sequencing run on August 21, 2016.
# Deep sequencing of longitudinal flu samples from two patients, in duplicate.
# Script is designed to be run from the top level of the Github repository.
# Script is designed to be run from the top-level directory of the Github repository.

# Load modules.
module load modules modules-init modules-gs
module load bowtie2/2.2.3
# Folder in which to save intermediate files.
dir="nobackup"
projectdir="SCCA"

# Raw sequence reads and reference sequence.
fastq1="$1"
fastq2="$2"
reference="$3"

# Parse the sample name.
sample=${fastq1##*/}
sample=${sample%%_*}
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

# Filter out all reads that align to the human genome.
# Retain unmapped reads for further analysis.
# Use a bowtie2 index built for human_g1k_hs37d5, a version of the hg19 reference genome.
:<<END
bowtie2 --very-fast-local --un-conc-gz ${dir}/${projectdir}/${sample}-filtered \
  -p 4 \
  -X 1000 \
  -x ${reference} \
  -1 ${fastq1} \
  -2 ${fastq2} \
  -S ${dir}/${projectdir}/${sample}-filtering.sam \
  2> ${dir}/${projectdir}/${sample}-filtering.bt2.log
END
  
# Rename output untrimmed read files to have a .fastq.gz extension.
mv ${dir}/${projectdir}/${sample}-filtered.1 ${dir}/${projectdir}/${sample}-filtered.1.fastq.gz
mv ${dir}/${projectdir}/${sample}-filtered.2 ${dir}/${projectdir}/${sample}-filtered.2.fastq.gz