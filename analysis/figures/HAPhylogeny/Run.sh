# Script to generate HA phylogeny.
# This script is meant to be run from the top level of the Github repository.

dir="analysis/figures/HAPhylogeny"

# Concatenate sequences downloaded from GISAID with the patient consensus sequences.
cat ${dir}/H3N2-SCCA-HA.fasta \
  ${dir}/H3N2-GISAID-USA-2004-2007-HA.fasta > \
  ${dir}/H3N2-GISAID-USA-2004-2007-patients-HA.fasta

# Construct the HA phylogeny. 
# Use needle to pairwise align the sequences to the Brisbane-2007 HA reference.
# Then, use RAxML to construcct the phylogeny.
python ${dir}/ConstructPhylogeny.py \
  reference/H3N2-Brisbane-2007-4-HA.fasta \
  ${dir}/H3N2-GISAID-USA-2004-2007-patients-HA.fasta \
  ${dir}

# Remove the intermediate FASTA file.
rm ${dir}/H3N2-GISAID-USA-2004-2007-patients-HA.fasta

# Move the RAxML output to the main directory.
mv RAxML* ${dir}

# Extract the names of the sequences in the phylogeny.
for f in ${dir}/*.reduced
do
  sample=${f##*/}
  sample=${sample%%.*}
  tail -n +2 ${f} | cut -f1 -d' ' > ${dir}/${sample}.data
done