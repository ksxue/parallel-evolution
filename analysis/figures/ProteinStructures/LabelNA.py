"""
This script takes in a set of sites in HA colors them on the protein structure.
It is intended to be run from the top level of the Github repository.
"""

ParallelColor="orange"
NSColor="orange"
SColor="lightteal"

# Import protein structure and perform basic formatting operations.
def LabelNA(SiteFile, Color, OutFile):

    # Display one NA subunit.
    cmd.fetch("2bat")
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background",0)
    cmd.show("cartoon","chain A+B")
    cmd.color("gray", "chain A+B")
    cmd.rotate("y", -75)
    cmd.rotate("x", -35)

    # Import data on sites to be labeled.
    with open(SiteFile, "r") as f:
        Sites=f.read().splitlines()

    # Split lines of file.
    for i in range(0,len(Sites)):
        Sites[i]=Sites[i].split()

    # Iterate through the list of sites and color
    # based on whether the mutation is synonymous or nonsynonymous
    # and whether it occurs in multiple patients.
    for i in range(0,len(Sites)):

        # Determine the chain identity and residue number.
        residue=Sites[i][0]
        chain="A"


        # Send the command to color the residues of interest
        # and display them as spheres.
        selection = "resi " + residue + " and chain " + chain
        cmd.color(Color, selection)
        cmd.show("spheres", selection)

    # Export a png of the sites of interest.
    cmd.png(OutFile, width=1000, dpi=200, ray=1)

    # Reinitialize pymol.
    cmd.reinitialize()

LabelNA("analysis/figures/ProteinStructures/6-NA-NS.data",NSColor,
        "analysis/figures/ProteinStructures/6-NA-NS.png")
LabelNA("analysis/figures/ProteinStructures/6-NA-S.data",SColor,
        "analysis/figures/ProteinStructures/6-NA-S.png")
LabelNA("analysis/figures/ProteinStructures/6-NA-parallel.data",ParallelColor,
        "analysis/figures/ProteinStructures/6-NA-parallel.png")
