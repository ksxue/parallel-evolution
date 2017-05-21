library(dplyr)

args<-commandArgs(TRUE)
NumSimulations<-as.integer(args[1])

set.seed(0)

# Read in list of longitudinal sites for each patient.
Data <- read.table(
  "analysis/figures/LongitudinalFrequencies/LongitudinalVariants.data",
  header=TRUE, stringsAsFactors = FALSE)

# Analyze only variants that are in HA, and classify variants
# as synonymous or nonsynonymous.
Data <- Data %>% mutate(AAChange=ifelse(InitAA==DerAA,"S","NS")) %>%
  dplyr::select(Patient, Gene, Codon, AAChange) %>% distinct() %>%
  filter(Gene %in% c("4-HA")) %>% ungroup()

# Read in list of antigenic site annotations.
# Note that these annotations are already in HA numbering.
AntigenicSites <- (read.table(
  "analysis/figures/Permutations/ChenLeeAntigenicSites.txt",
  header=TRUE, stringsAsFactors = FALSE
) %>% filter(!is.na(Antigenic_site)))$AA

# Calculate the number of nonsynonymous within-host mutations 
# that fall in antigenic sites.
# Count each variant separately based on its number of occurrences,
# i.e. do not collapse multiple variants that occur at the same site.
Data <- Data %>% mutate(AntigenicSite=ifelse(Codon %in% AntigenicSites, 1, 0))
WithinAntigenic <- (Data %>% group_by(AAChange) %>% 
                      summarize(NumAntigenic=sum(AntigenicSite==1),
                                NumTotal=n()))

# Read in the total length of the HA gene, in codons.
HALength <- (read.table("analysis/figures/SequencingDepth/GenomeGeneSummary.data",
                        header=TRUE, stringsAsFactors = FALSE) %>%
               filter(Gene=="4-HA"))$NumCodons

# For each type of site, synonymous and nonsynonymous,
# draw sites without replacement for each patient from the set of HA sites,
# matching the distribution of variants per patient.
# Determine what proportion fall within antigenic sites under this classification.
RunPermutation <- function(sitetype){
  
  # Filter for the desired site type.
  DataSite <- Data %>% filter(AAChange==sitetype) 
  
  # Determine the distribution of sites per patient.
  NumSitesPerPatient <- (DataSite %>% group_by(Patient) %>%
                           summarize(NumSites=length(unique(Codon))))$NumSites
  
  # Sample sites and return number that are found in antigenic sites.
  SampleSites <- function(){
    Sites <- unlist(sapply(NumSitesPerPatient, 
             function(x) sample(seq(1,HALength),x, replace=FALSE)))
    return(length(Sites[Sites %in% AntigenicSites]))
  }
  
  # Run the specified number of permutations.
  Permutations <- as.data.frame(table(replicate(NumSimulations, SampleSites())))
  colnames(Permutations) <- c("NumAntigenic","Frequency")
  Permutations$NumAntigenic <- as.integer(as.character(Permutations$NumAntigenic))
  Permutations$Frequency <- as.integer(as.character(Permutations$Frequency))
  
  # Calculate a p-value, which is the number of simulations with equal or
  # greater overlap with antigenic sites than observed in the empirical data.
  NumEmpirical <- (WithinAntigenic  %>% filter(AAChange==sitetype))$NumAntigenic
  pvalue <- (Permutations %>% 
               filter(NumAntigenic < NumEmpirical) %>% 
               mutate(pvalue=1-sum(Frequency)/NumSimulations) %>%
               mutate(Gene="4-HA", AAChange=sitetype, 
                      NumSimulations=NumSimulations) %>%
               dplyr::select(Gene, AAChange, NumSimulations, pvalue)) %>% distinct()
  return(pvalue)
}

PermuteAll <- do.call(rbind, lapply(c("NS","S"), RunPermutation))

# Export the data.
write.table(PermuteAll, "analysis/figures/Permutations/AntigenicSitesPermutation.data",
            row.names=FALSE, quote=FALSE)

