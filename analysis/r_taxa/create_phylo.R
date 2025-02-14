#Author: Fangxi Xu
#Date: 10/05/2023
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")
#load required packages
library(readxl)
library(dplyr)
library(tidyverse)
library(phyloseq)
############################################################################
#STEP1 - prepare dataset
############################################################################
#read the count data(convert txt to xlsx, delete first row, and replace | with ; before importing data)
profile <- read_excel("/Volumes/T7/Analysis/aging_manuscript_oral/from_hpc/Metaphlan4/sgb_profile_merged_r_uc.xlsx")
#convert to dataframe
profile = as.data.frame(profile)
#extract sgb level data only
profile_sgb <- subset(profile, grepl('t__', clade_name))
#separate the clade names 
profile_sgb = as.data.frame(separate(profile_sgb, clade_name, 
                                     into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","SGB"), sep=";"))
#add a column representing SGB ID
profile_sgb$clade_taxid <- profile_sgb$SGB
#move the "clade_taxid" column up to the front
profile_sgb <- profile_sgb %>% dplyr::select("clade_taxid", everything())
profile_sgb$clade_taxid <-  gsub(pattern = "t__", replacement = "", x = profile_sgb$clade_taxid) 
profile_sgb$clade_taxid <- gsub("SGB|_group","",profile_sgb$clade_taxid)
####
#create taxonomy data table
taxonomy <- profile_sgb[,c(1:9)]
taxonomy <- taxonomy %>% remove_rownames %>% column_to_rownames(var="clade_taxid")
#clean the species column - remove GGB#
taxonomy$Species = gsub(pattern = "GGB.*_", replacement = "", x = taxonomy$Species)
#clean the taxonomy data table
taxonomy_clean <- taxonomy %>% 
  mutate(Class = ifelse(grepl("CFG",taxonomy$Class),paste0(Class, "_",Phylum),Class)) %>%
  mutate(Order = ifelse(grepl("OFG",taxonomy$Order),paste0(Order, "_",Class),Order)) %>%
  mutate(Family = ifelse(grepl("FGB",taxonomy$Family),paste0(Family, "_",Class),Family)) %>%
  mutate(Genus = ifelse(grepl("GGB",taxonomy$Genus),paste0(Genus, "_",Family),Genus)) %>%
  mutate(Species = ifelse(grepl("SGB",taxonomy$Species),paste0(Species, "_",Genus),Species))
taxon <- as.matrix(taxonomy_clean)
####
#clean count table so the SampleID match with the mapping file
#clean count table so the SampleID match with the mapping file
biom_table <- profile_sgb
biom_table <- biom_table %>% remove_rownames %>% column_to_rownames(var="clade_taxid")
biom_table <- biom_table[,-c(1:8)]
names(biom_table) = gsub(pattern = "metaphlan_r_uc_", replacement = "", x = names(biom_table))
names(biom_table) = gsub(pattern = "_S*..", replacement = "", x = names(biom_table))
biom = as.data.frame(sapply(biom_table, as.numeric),row.names = rownames(biom_table))
####
#Metaphlan sgb phylogenetic tree
library("ape")
tree=ape::read.tree("/Volumes/T7/Analysis/phylo_tree/vOct22/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")
#tree$tip.label
####
# Read in metadata
meta <- read_excel("/Volumes/T7/Analysis/aging_manuscript_oral/from_hpc/mapping/mapping_250k.xlsx")
#make SampleID as row names
meta <- meta %>% remove_rownames %>% column_to_rownames(var="SampleID")

############################################################################
#STEP2 - Import all as phyloseq objects
############################################################################
OTU <- otu_table(biom, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxon)
META <- sample_data(meta)
# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)
# Same sample names
sample_names(OTU)
sample_names(META)
#merge and create phyloseq object
phylo <- phyloseq(OTU, TAX, META)
phylo
####
phy_tree(phylo)=tree
phylo
############################################################################
#STEP3 - Subset phylo based on SampleType
#clean the phyloseq object to remove controls, ineligible samples, and unwanted taxa
############################################################################
##STEP3.1 select bacteria only
phylo
phylo.bacteria <- subset_taxa(phylo, Kingdom == "k__Bacteria")
phylo.bacteria
###
#Change the order of levels for analysis and plotting
sample_data(phylo.bacteria)$Group <- factor(sample_data(phylo.bacteria)$Group)
is.factor(sample_data(phylo.bacteria)$Group)
levels(sample_data(phylo.bacteria)$Group)
# Change the order of levels
sample_data(phylo.bacteria)$Group <- factor(sample_data(phylo.bacteria)$Group, 
                                          levels = c("WT-Young","WT-Old","KO-Young","KO-Old"),
                                          labels = c("WT-Young","WT-Old","KO-Young","KO-Old"))
sample_data(phylo.bacteria)$Group
############################################################################
###
##STEP3
#remove unwanted sample
phylo.oral <- subset_samples(phylo.bacteria, SeqID != "MYWT-6_S91" & SeqID != "MKO10_S29")
phylo.oral
#prune OTUs that are not present in at least one sample
phylo.oral <- prune_taxa(taxa_sums(phylo.oral) > 0, phylo.oral)
phylo.oral
#filter out taxa present in less than 2 samples 
#phylo.oral <- filter_taxa(phylo.oral, function (x) {sum(x > 0) > 1}, prune=TRUE)
#phylo.oral
#re normalize the relative abundance so it sums to 100%
phylo.oral.rel = transform_sample_counts(phylo.oral, function(x) (x / sum(x))*100 )
phylo.oral.rel
saveRDS(phylo.oral.rel, "./Robjects/myphyloseq.rds")
############################################################################
#STEP4 - subset phylo for pairwise comparison - Oral
############################################################################
####
#Oral
#subset WT
phylo.oral.WT <- subset_samples(phylo.oral.rel, Genotype == "WT")
phylo.oral.WT
phylo.oral.WT <- prune_taxa(taxa_sums(phylo.oral.WT) > 0, phylo.oral.WT)
phylo.oral.WT
saveRDS(phylo.oral.WT, "./Robjects/myphyloseq1.rds")
#subset KO
phylo.oral.KO <- subset_samples(phylo.oral.rel, Genotype == "KO")
phylo.oral.KO
phylo.oral.KO <- prune_taxa(taxa_sums(phylo.oral.KO) > 0, phylo.oral.KO)
phylo.oral.KO
saveRDS(phylo.oral.KO, "./Robjects/myphyloseq2.rds")
#subset Young
phylo.oral.Young <- subset_samples(phylo.oral.rel, Aging == "3-month")
phylo.oral.Young
phylo.oral.Young <- prune_taxa(taxa_sums(phylo.oral.Young) > 0, phylo.oral.Young)
phylo.oral.Young
saveRDS(phylo.oral.Young, "./Robjects/myphyloseq3.rds")
#subset Old
phylo.oral.Old <- subset_samples(phylo.oral.rel, Aging == "Old")
phylo.oral.Old
phylo.oral.Old <- prune_taxa(taxa_sums(phylo.oral.Old) > 0, phylo.oral.Old)
phylo.oral.Old
saveRDS(phylo.oral.Old, "./Robjects/myphyloseq4.rds")
#subset KO-Old and WT-Young
phylo.oral.KOo.WTy <- subset_samples(phylo.oral.rel, Group == "WT-Young"|Group == "KO-Old")
phylo.oral.KOo.WTy
phylo.oral.KOo.WTy <- prune_taxa(taxa_sums(phylo.oral.KOo.WTy) > 0, phylo.oral.KOo.WTy)
phylo.oral.KOo.WTy
saveRDS(phylo.oral.KOo.WTy, "./Robjects/myphyloseq5.rds")
############################################################################
