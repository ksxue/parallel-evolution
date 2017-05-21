library(dplyr)

args<-commandArgs(TRUE)
NumSimulations<-as.integer(args[1])

set.seed(0)

# Read in list of longitudinal sites for each patient.
Data <- read.table(
  "analysis/figures/LongitudinalFrequencies/LongitudinalVariants.data",
  header=TRUE, stringsAsFactors = FALSE)

# For each variant, summarize whether it is synonymous or nonsynonymous.
# Note that HA numbering has already been shifted to match the H3 numbering,
# so no further modification is needed.
Data <- Data %>% mutate(AAChange=ifelse(InitAA==DerAA,"S","NS")) %>%
  dplyr::select(Patient, Gene, Codon, AAChange) %>% distinct() %>%
  filter(Gene %in% c("4-HA","6-NA"))

# Read in list of distances for HA and NA.
# These are the distances of each residue from sialic acid.
# Note that the NA residues are numbered beginning at 82.
Distances <- rbind(read.table(
  "analysis/figures/ProteinStructures/4HMG-SialicAcidDistances.data",
  header=TRUE, stringsAsFactors = FALSE) %>% mutate(Gene="4-HA"),
  read.table(
    "analysis/figures/ProteinStructures/2BAT-SialicAcidDistances.data",
    header=TRUE, stringsAsFactors = FALSE) %>% 
    mutate(Gene="6-NA", AANumber=AANumber+81))

# Write a function that, given a gene name and codon number,
# returns the distance of that residue from sialic acid,
# using the Distances dataframe above.
ReturnDistance <- function(gene, codon){
  return((Distances %>% filter(Gene==gene, AANumber==codon))$Distance[1])
}
Data$Distance <- mapply(ReturnDistance, Data$Gene, Data$Codon)

# Write a function to permute the distances in a given protein
# by drawing from the distribution in line with the number
# of variants per patient.
# Compare the median distance to the median for the patient variants.
# Do this for both synonymous and nonsynonymous variants.
RunPermutation <- function(gene, aachange, numpermutations){
  
  # Select the patient sites of the specified gene and aachange type.
  # Ensure that they are sites that are also in the crystal structure.
  GeneDistances <- Distances %>% filter(Gene==gene)
  Within <- Data %>% 
    filter(Gene==gene, AAChange==aachange, 
           Codon %in% GeneDistances$AANumber)
  
  # Determine the number of sites present in each patient.
  NumSitesPerPatient <- (Within %>% group_by(Patient) %>%
                           summarize(NumSites=n()))$NumSites
  
  # Draw a set of random sites for each patient
  SampleSites <- function(){
    return(median(unlist(
      sapply(NumSitesPerPatient,
             function(x) sample(GeneDistances$Distance, x,
                                replace=FALSE))
    )))
  }
  
  # Run the specified number of permutations.
  Permutations <- sort(replicate(numpermutations, SampleSites()))
  
  # Determine the p-value - that is, the proportion of permutations
  # that have a median distance equal to or lower than what is
  # empirically observed.
  pvalue <- length(Permutations[Permutations<=median(unlist(Within$Distance))]) /
    numpermutations
  
  return(Within %>% mutate(NumPermutations=numpermutations,
                           pvalue=pvalue) %>%
           dplyr::select(Gene, AAChange, NumPermutations, pvalue) %>%
           distinct())
}

# Run permutation tests for HA and NA for NS and S sites.
Genes <- c("4-HA", "6-NA")
AAChange <- c("NS","S")
PermuteAll <- do.call(rbind,
  mapply(RunPermutation,
         rep(Genes, length(AAChange)),
         rep(AAChange, each=length(Genes)),
         rep(NumSimulations, length(Genes)*length(AAChange)),
         SIMPLIFY=FALSE)
)

# Also run permutation tests for parallel sites in HA.
Within <- Data %>% filter(Gene=="4-HA") %>%
  group_by(Gene, Codon) %>% filter(n()>1) %>%
  ungroup() %>% dplyr::select(-Patient) %>% distinct()
GeneDistances <- Distances %>% filter(Gene=="4-HA")
SampleSites <- function(){
  return(median(sample(GeneDistances$Distance, nrow(Within), replace=FALSE)))
}
  
Permutations <- unlist(replicate(NumSimulations, SampleSites()))

# Determine the p-value - that is, the proportion of permutations
# that have a median distance equal to or lower than what is
# empirically observed.
pvalue <- length(Permutations[Permutations<=
                                median(unlist(Within$Distance))]) /
  NumSimulations

PermuteAll <- rbind(PermuteAll,
                    Within %>% mutate(NumPermutations=NumSimulations,
                                      pvalue=pvalue, AAChange="parallel") %>%
                      dplyr::select(Gene, AAChange, NumPermutations, pvalue) %>%
                      distinct())

# Output table of p-values.
write.table(PermuteAll, 
            "analysis/figures/ProteinStructures/PermuteDistances.data",
            row.names=FALSE, quote=FALSE)
