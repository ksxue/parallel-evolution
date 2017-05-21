library(dplyr)

args<-commandArgs(TRUE)
NumSimulations<-as.integer(args[1])

set.seed(0)

# Read in list of longitudinal sites for each patient.
Data <- read.table(
  "analysis/figures/LongitudinalFrequencies/LongitudinalVariants.data",
  header=TRUE, stringsAsFactors = FALSE)

# Read in list of genes and sites.
GenomeSummary <- read.table(
  "analysis/figures/SequencingDepth/GenomeGeneSummary.data",
  header=TRUE, stringsAsFactors = FALSE
)


# Write a function that, given a gene, its number of codons,
# and a specified proportion of mutable sites,
# assess the patient data and draw the appropriate number of nonsynonymous sites
# per patient for the specified gene from the set of mutable sites, 
# given the number that appear in the actual data.
RunPermutation <- function(gene, genelength, numpermutations, siteproportion){
  
  # Determine the number of nonsynonymous sites in the specified gene per patient,
  # and the total number of unique nonsynonymous sites in the specified gene
  # across all patients.
  NumSitesPerPatient <- (Data %>% group_by(Patient) %>%
                           filter(Gene==gene, InitAA!=DerAA) %>%
                           summarize(NumSites=length(unique(Codon))))$NumSites
  NumUniqueSitesTotal <- length(unique((Data %>% 
                                          filter(Gene==gene, InitAA!=DerAA))$Codon))
  
  # If the number of unique sites equals the number of total sites,
  # then exit the function without performing simulations,
  # since the p-value will be 1.
  if(sum(NumSitesPerPatient)==NumUniqueSitesTotal){
    return()
  }
  
  # Draw a set of random sites for each patient.
  # For each patient independently, draw the same number of sites, without replacement,
  # as were annotated in the real data.
  # Calculate the total number of unique sites drawn as a result.
  SampleSites <- function(x){
    return(length(unique(unlist(
      sapply(NumSitesPerPatient, 
             function(x) sample(seq(1,round(genelength*siteproportion)),
                                x, replace=FALSE))))))
  }
  Permutations <- as.data.frame(table(replicate(numpermutations,SampleSites())))
  colnames(Permutations) <- c("NumUniqueSimulated","Frequency")
  Permutations <- Permutations %>% 
    mutate(Gene=gene, NumUniqueData=NumUniqueSitesTotal, 
           NumSitesTotal=sum(NumSitesPerPatient)) %>%
    dplyr::select(Gene, NumSitesTotal,
                  NumUniqueData, NumUniqueSimulated, Frequency)
  Permutations$NumUniqueData <- as.numeric(as.character(Permutations$NumUniqueData))
  Permutations$NumUniqueSimulated <- 
    as.numeric(as.character(Permutations$NumUniqueSimulated))
  
  # Calculate a p-value, which is the proportion of simulations
  # that give the same number or fewer unique sites.
  pvalue <- (Permutations %>%
      filter(NumUniqueSimulated > NumUniqueData) %>%
      group_by(Gene) %>% 
      mutate(pvalue=1-(sum(Frequency))/numpermutations))$pvalue[1]
  Permutations <- Permutations %>%
    mutate(Gene=gene,
           GeneLength=genelength,
           SiteProportion=siteproportion,
           NumSimulations=numpermutations,
           pvalue=pvalue)
    
  return(Permutations)
}

# Set a vector of site proportions to be tested.
# The proportion restricts the number of sites considered to be mutable.
# For instance, a site proportion of 0.2 means that only 20%
# of the sites in the protein are considered mutable and
# are included in the permutation test.
SiteProportions <- seq(0.05,1,0.05)

# Run permutation tests for each gene, across the specified
# proportions of mutable sites.
# The way that the permutation function is written,
# a gene will not show up if the p-value is 1.
PermuteAll <- do.call(rbind, mapply(RunPermutation, 
                      rep(GenomeSummary$Gene, each=length(SiteProportions)), 
                      rep(GenomeSummary$NumCodons, each=length(SiteProportions)), 
                      rep(NumSimulations, nrow(GenomeSummary)*length(SiteProportions)),
                      rep(SiteProportions, nrow(GenomeSummary)),
                      SIMPLIFY=FALSE)
)

# Export the summaries of all of the simulations.
write.table(PermuteAll, 
            "analysis/figures/Permutations/ParallelWithinPermutation.data",
            row.names=FALSE, quote=FALSE)

# Read in proportion of globally mutable sites for each gene
# and calculate the appropriate p-value.
# This is the proportion of sites at which there is not a 100% consensus.
VariableSites <- read.table("analysis/figures/Permutations/GlobalVariableSites.data",
                            header=TRUE, stringsAsFactors = FALSE)

PermuteGlobal <- do.call(rbind, mapply(RunPermutation, 
                                       GenomeSummary$Gene,
                                       GenomeSummary$NumCodons,
                                       rep(NumSimulations, nrow(GenomeSummary)),
                                       VariableSites$VariableProportion,
                                       SIMPLIFY=FALSE))

write.table(PermuteGlobal, 
            "analysis/figures/Permutations/ParallelWithinPermutation-GlobalValues.data",
            row.names=FALSE, quote=FALSE)