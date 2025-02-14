#set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")
# create plot folder
path <- "./plots"
#load phyloseq object
phylo <- readRDS("./Robjects/myphyloseq.rds")
#load required library
library("ggpubr")
##############################################
source("./functions/functions-beta-diversity.R")
##############################################
#Beta Diversity
#assign colors
my.colors <- c("WT-Young"="#f39a8f", "WT-Old"="#d95475", "KO-Young"="#f2d58c", "KO-Old"="#ef811e")
#assign shapes
group.shapes <- c("WT-Young"=17, "WT-Old"=19, "KO-Young"=17, "KO-Old"=19)
################################################
#plot for beta diversity - Oral
#Bray-Curtis Beta diversity
p <- plot.bray.curtis(phylo, group = "Group")
p
#save plot to pdf
pdf(file.path(path, "beta_diversity_bray_curtis_all_groups.pdf"),width=4.5,height=3)
p
dev.off()
###############################################
phyloseq_obj1 <- readRDS("./Robjects/myphyloseq1.rds")
phyloseq_obj2 <- readRDS("./Robjects/myphyloseq2.rds")
phyloseq_obj3 <- readRDS("./Robjects/myphyloseq3.rds")
phyloseq_obj4 <- readRDS("./Robjects/myphyloseq4.rds")
phyloseq_obj5 <- readRDS("./Robjects/myphyloseq5.rds")
p1 <- plot.bray.curtis.pairwise(phyloseq_obj1, group = "Group")
p1
p2 <- plot.bray.curtis.pairwise(phyloseq_obj2, group = "Group")
p2
p3 <- plot.bray.curtis.pairwise(phyloseq_obj3, group = "Group")
p3
p4 <- plot.bray.curtis.pairwise(phyloseq_obj4, group = "Group")
p4
p5 <- plot.bray.curtis.pairwise(phyloseq_obj5, group = "Group")
p5
pdf(file.path(path, "beta_diversity_bray_curtis_pairwise.pdf"),width=14,height=3)
ggarrange(p1, p2, p3, p4, p5, ncol = 5, widths = c(0.2,0.2,0.2,0.2,0.2))
dev.off()

################################################
#Permanova
BC=TRUE
sig <- test_permanova(phylo, group = "Group")
sig#0.001

sig1<- test_permanova(phyloseq_obj1, group = "Group")
sig1#0.001

sig2<- test_permanova(phyloseq_obj2, group = "Group")
sig2#0.01

sig3<- test_permanova(phyloseq_obj3, group = "Group")
sig3#0.6

sig4<- test_permanova(phyloseq_obj4, group = "Group")
sig4#0.002

sig5<- test_permanova(phyloseq_obj5, group = "Group")
sig5#0.001

################################################
################################################
#plot for beta diversity - Oral
#Weighted-Unifrac Beta diversity
p <- plot.wuni(phylo, group = "Group")
p
#save plot to pdf
#png(file.path(path1, "beta_diversity_oral.pdf"), width=22,height=20,units="cm",res=1200)
pdf(file.path(path, "beta_diversity_weighted_unifrac_all_groups.pdf"),width=4.5,height=3)
p
dev.off()

p1 <- plot.wuni.pairwise(phyloseq_obj1, group = "Group")
p1
p2 <- plot.wuni.pairwise(phyloseq_obj2, group = "Group")
p2
p3 <- plot.wuni.pairwise(phyloseq_obj3, group = "Group")
p3
p4 <- plot.wuni.pairwise(phyloseq_obj4, group = "Group")
p4
p5 <- plot.wuni.pairwise(phyloseq_obj5, group = "Group")
p5
pdf(file.path(path, "beta_diversity_weighted_unifrac_pairwise.pdf"),width=14,height=3)
ggarrange(p1, p2, p3, p4, p5, ncol = 5, widths = c(0.2,0.2,0.2,0.2,0.2))
dev.off()

BC=FALSE
sig <- test_permanova(phylo, group = "Group")
sig#0.001

sig1<- test_permanova(phyloseq_obj1, group = "Group")
sig1#0.003

sig2<- test_permanova(phyloseq_obj2, group = "Group")
sig2#0.03

sig3<- test_permanova(phyloseq_obj3, group = "Group")
sig3#0.6

sig4<- test_permanova(phyloseq_obj4, group = "Group")
sig4#0.002

sig5<- test_permanova(phyloseq_obj5, group = "Group")
sig5#0.006
