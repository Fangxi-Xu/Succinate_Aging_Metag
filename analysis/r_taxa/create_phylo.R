# Author: Fangxi Xu
# Date: 10/05/2023

# Set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")

# Load required packages
library(readxl)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(ape)

############################################################################
# STEP 1 - Prepare dataset
############################################################################

# Read and process the count data
profile <- read_excel("/Volumes/T7/Analysis/aging_manuscript_oral/from_hpc/Metaphlan4/sgb_profile_merged_r_uc.xlsx") %>%
  as.data.frame() %>%
  subset(grepl('t__', clade_name)) %>%
  separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB"), sep = ";") %>%
  mutate(clade_taxid = gsub("t__|SGB|_group", "", SGB),
         clade_taxid = gsub(pattern = "t__", replacement = "", clade_taxid)) %>%
  select(clade_taxid, everything())

# Create and clean taxonomy data table
taxonomy <- profile[, 1:9] %>%
  column_to_rownames(var = "clade_taxid") %>%
  mutate(Class = ifelse(grepl("CFG", Class), paste0(Class, "_", Phylum), Class),
         Order = ifelse(grepl("OFG", Order), paste0(Order, "_", Class), Order),
         Family = ifelse(grepl("FGB", Family), paste0(Family, "_", Class), Family),
         Genus = ifelse(grepl("GGB", Genus), paste0(Genus, "_", Family), Genus),
         Species = ifelse(grepl("SGB", Species), paste0(Species, "_", Genus), Species)) %>%
  as.matrix() %>%
  t()

# Clean and reshape count table
biom_table <- profile[, -c(1:8)] %>%
  column_to_rownames(var = "clade_taxid") %>%
  rename_with(~ gsub("metaphlan_r_uc_|_S*..", "", .)) %>%
  sapply(as.numeric) %>%
  as.data.frame(row.names = rownames(.))

# Read in phylogenetic tree and metadata
tree <- read.tree("/Volumes/T7/Analysis/phylo_tree/vOct22/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")
meta <- read_excel("/Volumes/T7/Analysis/aging_manuscript_oral/from_hpc/mapping/mapping_250k.xlsx") %>%
  column_to_rownames(var = "SampleID")

############################################################################
# STEP 2 - Import all as phyloseq objects
############################################################################

# Construct phyloseq objects
OTU <- otu_table(biom, taxa_are_rows = TRUE)
TAX <- tax_table(taxon)
META <- sample_data(meta)

# Merge and create phyloseq object
phylo <- phyloseq(OTU, TAX, META)
phylo <- merge_phyloseq(phylo, phy_tree(tree))

############################################################################
# STEP 3 - Subset phylo based on SampleType and cleanup
############################################################################

# Subset to bacteria only and reorder factor levels
phylo.bacteria <- subset_taxa(phylo, Kingdom == "k__Bacteria") %>%
  transform_sample_counts(function(x) (x / sum(x)) * 100) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  filter_taxa(., function(x) sum(x > 0) > 1, prune = TRUE)

# Reorder groups for analysis and plotting
sample_data(phylo.bacteria)$Group <- factor(sample_data(phylo.bacteria)$Group, 
                                            levels = c("WT-Young", "WT-Old", "KO-Young", "KO-Old"))

# Save the cleaned-up phyloseq object
saveRDS(phylo.bacteria, "./Robjects/myphyloseq.rds")

############################################################################
# STEP 4 - Subset phylo for pairwise comparisons
############################################################################

subset_and_save <- function(phylo, condition, filename) {
  subset <- subset_samples(phylo, condition) %>%
    prune_taxa(taxa_sums(.) > 0, .)
  saveRDS(subset, file.path("./Robjects", filename))
  return(subset)
}

phylo.oral.WT <- subset_and_save(phylo.bacteria, Genotype == "WT", "myphyloseq1.rds")
phylo.oral.KO <- subset_and_save(phylo.bacteria, Genotype == "KO", "myphyloseq2.rds")
phylo.oral.Young <- subset_and_save(phylo.bacteria, Aging == "3-month", "myphyloseq3.rds")
phylo.oral.Old <- subset_and_save(phylo.bacteria, Aging == "Old", "myphyloseq4.rds")
phylo.oral.KOo.WTy <- subset_and_save(phylo.bacteria, Group %in% c("WT-Young", "KO-Old"), "myphyloseq5.rds")
