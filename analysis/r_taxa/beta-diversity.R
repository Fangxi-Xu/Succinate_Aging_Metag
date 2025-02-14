# Author: Fangxi Xu
# Date: 10/05/2023

# Set working directory and load required libraries
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")
library("ggpubr")

# Load phyloseq object
phylo <- readRDS("./Robjects/myphyloseq.rds")

# Source functions for beta diversity analysis
source("./functions/functions-beta-diversity.R")

# Define plot path
path <- "./plots"

# Assign colors and shapes for plotting
my.colors <- c("WT-Young"="#f39a8f", "WT-Old"="#d95475", "KO-Young"="#f2d58c", "KO-Old"="#ef811e")
group.shapes <- c("WT-Young"=17, "WT-Old"=19, "KO-Young"=17, "KO-Old"=19)

##############################################
# Beta Diversity Analysis
##############################################

# Bray-Curtis Beta diversity plot
p <- plot.bray.curtis(phylo, group = "Group")
pdf(file.path(path, "beta_diversity_bray_curtis_all_groups.pdf"), width=4.5, height=3)
print(p)
dev.off()

# Load additional phyloseq objects
phyloseq_obj1 <- readRDS("./Robjects/myphyloseq1.rds")
phyloseq_obj2 <- readRDS("./Robjects/myphyloseq2.rds")
phyloseq_obj3 <- readRDS("./Robjects/myphyloseq3.rds")
phyloseq_obj4 <- readRDS("./Robjects/myphyloseq4.rds")
phyloseq_obj5 <- readRDS("./Robjects/myphyloseq5.rds")

# Pairwise Bray-Curtis plots
p1 <- plot.bray.curtis.pairwise(phyloseq_obj1, group = "Group")
p2 <- plot.bray.curtis.pairwise(phyloseq_obj2, group = "Group")
p3 <- plot.bray.curtis.pairwise(phyloseq_obj3, group = "Group")
p4 <- plot.bray.curtis.pairwise(phyloseq_obj4, group = "Group")
p5 <- plot.bray.curtis.pairwise(phyloseq_obj5, group = "Group")

# Save pairwise plots to PDF
pdf(file.path(path, "beta_diversity_bray_curtis_pairwise.pdf"), width=14, height=3)
ggarrange(p1, p2, p3, p4, p5, ncol=5, widths=c(0.2, 0.2, 0.2, 0.2, 0.2))
dev.off()

###############################################
# Permanova Tests
###############################################

# Run Permanova on all phyloseq objects and display significance levels
sig <- test_permanova(phylo, group = "Group")
print(sig) # Print significance level for main object

# Repeat for other objects
sig1 <- test_permanova(phyloseq_obj1, group = "Group")
sig2 <- test_permanova(phyloseq_obj2, group = "Group")
sig3 <- test_permanova(phyloseq_obj3, group = "Group")
sig4 <- test_permanova(phyloseq_obj4, group = "Group")
sig5 <- test_permanova(phyloseq_obj5, group = "Group")

###############################################
# Weighted-Unifrac Beta Diversity Analysis
###############################################

# Plot Weighted-Unifrac for all groups and save
p <- plot.wuni(phylo, group = "Group")
pdf(file.path(path, "beta_diversity_weighted_unifrac_all_groups.pdf"), width=4.5, height=3)
print(p)
dev.off()

# Plot pairwise Weighted-Unifrac
p1 <- plot.wuni.pairwise(phyloseq_obj1, group = "Group")
p2 <- plot.wuni.pairwise(phyloseq_obj2, group = "Group")
p3 <- plot.wuni.pairwise(phyloseq_obj3, group = "Group")
p4 <- plot.wuni.pairwise(phyloseq_obj4, group = "Group")
p5 <- plot.wuni.pairwise(phyloseq_obj5, group = "Group")

# Save pairwise Weighted-Unifrac plots to PDF
pdf(file.path(path, "beta_diversity_weighted_unifrac_pairwise.pdf"), width=14, height=3)
ggarrange(p1, p2, p3, p4, p5, ncol=5, widths=c(0.2, 0.2, 0.2, 0.2, 0.2))
dev.off()
