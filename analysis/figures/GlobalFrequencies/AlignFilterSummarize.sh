# Script designed to align a set of homologous sequences to a reference,
# clean the alignment to remove outlier sequences,
# and summarize the alignment into a set of frequencies for each amino-acid
# in each year.

# Takes in arguments for
# 1 - reference sequence in FASTA format
# 2 - homologous sequences in multi-FASTA format
# 3 - directory path for alignment output files
# 4 - directory path for final output files
# 5 - sequence source: GISAID or Genbank
ref=$1
seqs=$2
aligned=$3
outdir=$4
source=$5
filestem=$6

echo ${filestem}

# Align sequences to the given reference using needle.
# Check first if the alignment file already exists.
if [ ! -f ${aligned} ]
then
  needle -asequence ${ref} -bsequence ${seqs} \
    -gapopen 10.0 -gapextend 0.5 \
    -outfile ${aligned} -aformat fasta
fi

# Parse the alignment file and calculate the distance
# of each sequence from the reference.
python bin/CalculateSequenceDistances-1.0.py ${ref} ${aligned} \
  ${outdir}/${filestem}-distances.data ${source}

# Analyze the distances and exclude sequences that are more or less than
# five interquartile ranges from the median amino-acid distance
# from the median distance to the reference sequence,
# binned by year.
Rscript analysis/figures/GlobalFrequencies/FilterOutlierSequences.R \
  ${outdir}/${filestem}-distances.data
  
# Analyze the sequence alignment, excluding sequences that contain indels
# or are on the list of sequence exclusions set previously.
python bin/AnalyzeAASiteFrequencies-2.0.py ${ref} ${aligned} \
  ${outdir}/${filestem}-distances-exclusions.data ${outdir}/${filestem}-frequencies.data ${source}