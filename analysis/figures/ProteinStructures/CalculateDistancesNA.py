# This script uses BioPDB to analyze the 4HMG HA structure
# and calculate the distance of each residue
# from the receptor-binding pocket,
# as represented to the sialic-acid residue.

from Bio.PDB import *

AminoAcids = ["ALA", "ARG", "ASN", "ASP", "CYS",
              "GLU", "GLN", "GLY", "HIS", "ILE",
              "LEU", "LYS", "MET", "PHE", "PRO",
              "SER", "THR", "TRP", "TYR", "VAL"]
ThreeToOne = {"ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C',
              "GLU":'E', "GLN":'Q', "GLY":'G', "HIS":'H', "ILE":'I',
              "LEU":'L', "LYS":'K', "MET":'M', "PHE":'F', "PRO":'P',
              "SER":'S', "THR":'T', "TRP":'W', "TYR":'Y', "VAL":'V'}

# Import protein structure.
parser = PDBParser()
structure = parser.get_structure('HA', 'reference/N2-2bat.pdb')


# Calculate the distance of all residues from the sialic-acid molecule.
# First, find the sialic acid residue in chain A.
for residue in structure[0]['A']:
    if residue.get_resname() == "SIA":
        sialic=residue

# Iterate through all of the residues in chain A (one monomer)
# Verify that they are amino acids and not heteromolecules.
# If so, then calculate their distance from every atom in the sialic acid
# and print the minimum to another file.
with open('analysis/figures/ProteinStructures/2BAT-SialicAcidDistances.data','w') as f:
    f.write("AANumber\tChainAANumber\tResidue\tDistance\n")
    chains = ["A"]
    aanumber = 1
    for chain in chains:
        chainresnumber = 1
        for residue in structure[0][chain]:
            if residue.get_resname() in AminoAcids:
                mindistance=10000000000
                for atom in sialic:
                    mindistance=min(mindistance, atom-residue['CA'])
                f.write("%d\t%d\t%s\t%f\n" %
                        (aanumber, chainresnumber,
                         ThreeToOne[residue.get_resname()], mindistance))
                aanumber += 1
                chainresnumber += 1
