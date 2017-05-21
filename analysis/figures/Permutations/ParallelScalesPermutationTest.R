library(dplyr)

args<-commandArgs(TRUE)
NumSimulations<-as.integer(args[1])

set.seed(0)

# Read in list of variable sites within patients.
Within <- read.table(
  "analysis/figures/LongitudinalFrequencies/LongitudinalVariants.data",
  header=TRUE, stringsAsFactors = FALSE)
Between <- read.table("analysis/figures/GlobalFrequencies/H3N2-GISAID-sites.data",
                      header=TRUE, stringsAsFactors = FALSE)
ParallelScales <- read.table(
  "analysis/figures/GlobalFrequencies/ParallelSites.data",
  header=TRUE, stringsAsFactors = FALSE)

# Read in list of genes and sites.
GenomeSummary <- read.table(
  "analysis/figures/SequencingDepth/GenomeGeneSummary.data",
  header=TRUE, stringsAsFactors = FALSE
)


# Read in proportion of globally mutable sites for each gene.
# This is the proportion of sites at which there is not a 100% consensus.
VariableSites <- read.table("analysis/figures/Permutations/GlobalVariableSites.data",
                            header=TRUE, stringsAsFactors = FALSE)
VariableSites <- merge(GenomeSummary, VariableSites, by=c("Gene"))

# Bin the variants and gene lengths into HA, NA, and other.
WithinVariants <- Within %>% group_by(Gene) %>%
  filter(InitAA!=DerAA) %>%
  summarize(NumVariants=length(unique(Codon))) %>%
  ungroup() %>% 
  mutate(Gene=ifelse(Gene %in% c("4-HA","6-NA"), Gene, "other")) %>%
  group_by(Gene) %>% summarize(NumVariants=sum(NumVariants))
BetweenVariants <- Between %>% group_by(Gene) %>%
  summarize(NumVariants=length(unique(Position))) %>%
  ungroup() %>% 
  mutate(Gene=ifelse(Gene %in% c("4-HA","6-NA"), Gene, "other")) %>%
  group_by(Gene) %>% summarize(NumVariants=sum(NumVariants))
ParallelVariants <- ParallelScales %>%
  mutate(Gene=ifelse(Gene %in% c("4-HA","6-NA"), Gene, "other")) %>%
  group_by(Gene) %>% summarize(NumVariants=n())
GenomeSummary <- GenomeSummary %>% 
  mutate(Gene=ifelse(Gene %in% c("4-HA","6-NA"), Gene, "other")) %>%
  group_by(Gene) %>% summarize(NumCodons=sum(NumCodons))
VariableSites <- VariableSites %>%
  mutate(Gene=ifelse(Gene %in% c("4-HA","6-NA"), Gene, "other")) %>%
  group_by(Gene) %>% 
  summarize(VariableProportion=sum(VariableSites)/sum(NumCodons))

# Write a function that, given a gene, its number of codons,
# and a specified proportion of mutable sites,
# assesses the within- and between-patient data.
# Draw sites randomly from the set of mutable sites in that gene,
# calculate the overlap with the number of within-host sites,
# and compare that overlap to the empirical overlap.
RunPermutation <- function(gene, genelength, numpermutations, siteproportion,
                           withinsites, betweensites, overlapsites){
  
  # If the proportion of mutable sites is lower than the 
  # number that need to be drawn, then quit.
  if(genelength*siteproportion<betweensites){
    return()
  }
  
  # Draw a set of random sites from the mutable sites without replacement
  # with set size indicated by betweensites.
  SampleSites <- function(){
    simbetween <- sample(seq(1,round(genelength*siteproportion)), 
                         betweensites, replace=FALSE)
    return(length(simbetween[simbetween<=withinsites]))
  }
  
  # Run many permutations and summarize the results.
  Permutations <- as.data.frame(table(replicate(numpermutations, SampleSites())))
  colnames(Permutations) <- c("OverlapSimulated", "Frequency")
  Permutations$OverlapSimulated <- as.numeric(
    as.character(Permutations$OverlapSimulated))
  
  # Calculate the p-value, which is the proportion of simulations
  # that give the same number or greater overlap than the
  # empirical data.
  pvalue <- (Permutations %>% mutate(Gene=gene) %>%
               group_by(Gene) %>%
               filter(OverlapSimulated < overlapsites) %>%
               mutate(pvalue=1-(sum(Frequency)/numpermutations)))$pvalue[1]
  if(is.na(pvalue)){
    pvalue<-1
  }
  return(Permutations %>% 
           mutate(Gene=gene,
                  GeneLength=genelength,
                  SiteProportion=siteproportion,
                  NumSimulations=numpermutations,
                  OverlapEmpirical=(overlapsites),
                  pvalue=pvalue))
}
  
# Set a vector of site proportions to be tested.
# The proportion restricts the number of sites considered to be mutable.
# For instance, a site proportion of 0.2 means that only 20%
# of the sites in the protein are considered mutable and
# are included in the permutation test.
SiteProportions <- seq(0.05,1,0.05)

# Run permutation tests for each gene, across the specified
# proportions of mutable sites.
PermuteAll <- do.call(rbind,
                      mapply(RunPermutation,
                             rep(GenomeSummary$Gene, each=length(SiteProportions)),
                             rep(GenomeSummary$NumCodons, each=length(SiteProportions)),
                             rep(NumSimulations, nrow(GenomeSummary)*length(SiteProportions)),
                             rep(SiteProportions, nrow(GenomeSummary)),
                             rep(WithinVariants$NumVariants, each=length(SiteProportions)),
                             rep(BetweenVariants$NumVariants, each=length(SiteProportions)),
                             rep(ParallelVariants$NumVariants, each=length(SiteProportions)),
                             SIMPLIFY=FALSE))

# Export the summaries of all of the simulations.
write.table(PermuteAll, 
            "analysis/figures/Permutations/ParallelScalesPermutation.data",
            row.names=FALSE, quote=FALSE)

PermuteGlobal <- do.call(rbind, mapply(RunPermutation, 
                                       GenomeSummary$Gene,
                                       GenomeSummary$NumCodons,
                                       rep(NumSimulations, nrow(GenomeSummary)),
                                       VariableSites$VariableProportion,
                                       WithinVariants$NumVariants,
                                       BetweenVariants$NumVariants,
                                       ParallelVariants$NumVariants,
                                       SIMPLIFY=FALSE))

write.table(PermuteGlobal, 
            "analysis/figures/Permutations/ParallelScalesPermutation-GlobalValues.data",
            row.names=FALSE, quote=FALSE)
