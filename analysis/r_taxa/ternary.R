#set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R-genefamily")
# create plot folder
#dir.create("./plots")
path <- "./plots"
#load required library
library(phyloseq)
library(microbiomeutilities)
library(microbiome)
library(knitr)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggtern)
library(viridis)
library(ggalt)
#load phyloseq object
phylo <- readRDS("./Robjects/myphyloseq.rds")
phylo.new <- subset_samples(phylo, Group  != "KO-Young")
#####
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(phylo.new),
               MARGIN = ifelse(taxa_are_rows(phylo.new), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phylo.new),
                    tax_table(phylo.new))
# Define prevalence threshold as 75% of total samples
prevalenceThreshold = 0.75 * nsamples(phylo.new)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keep = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
phylo.new = prune_taxa(keep, phylo.new)
phylo.new
#################################################################################

#################################################################################
tern <- prep_ternary(phylo.new, abund.thres= 0.0001, group="Group", level= "lowest")
names(tern)[names(tern) == 'OTUID'] <- 'feature'
names(tern)[names(tern) == 'KO-Old'] <- 'KO_Old'
names(tern)[names(tern) == 'WT-Old'] <- 'WT_Old'
names(tern)[names(tern) == 'WT-Young'] <- 'WT_Young'

library(readr)
res_WT_final_tern <- read_csv("res_WT_final_tern.csv")
View(res_WT_final_tern)
  
tern.1 <- merge(tern, res_WT_final_tern, by ="feature",all=T)
library(tidyr)
tern.1 %>% drop_na(KO_Old) -> tern.2
tern.2$tern <- tern.2$tern %>% replace_na('Others')
tern.2$tern <- factor(tern.2$tern, levels = c("WT_Old enriched gene families","WT_Old depleted gene families","Others"))
levels(tern.2$tern)
#head(tern_v1)
my.colors <- c("WT_Old depleted gene families"="#29b895", "WT_Old enriched gene families"="#b8294c", "Others"="#ffbf00")
my.shapes<- c("WT_Old depleted gene families"=17, "WT_Old enriched gene families"=18, "Others"=20)
my.sizes <- c("WT_Old depleted gene families"=8, "WT_Old enriched gene families"=8, "Others"=4)


p <- ggtern(data=tern.2, aes(KO_Old, WT_Old, WT_Young)) + 
  geom_point(aes(fill  = tern, color = tern, shape=tern), 
             size=8,
             alpha=0.6, 
             show.legend=T) + scale_color_manual(values = my.colors) + 
  scale_shape_manual(values = my.shapes)+ 
  theme(panel.border = element_rect(color = "black"),
        panel.grid =  element_line(color = "black",size=1)) +
  theme_showarrows() +  
  theme(plot.title = element_text(size = 12, color = "black"),
        axis.text.x=element_text(size=12, color = "black"),
        axis.title.x=element_text(size=12, color = "black"),
        axis.title.y=element_text(size=12, color = "black"),
        axis.text.y=element_text(size=12,  color = "black"),
        legend.title=element_blank(), 
        legend.text=element_text(size=12, color = "black"),
        legend.position="top",)

p
pdf(file.path(path, "tern_genefamily_WT-Old_new.pdf"), width = 8, height = 6)
p
dev.off()
###############################################################################
res_KO_final_tern <- read_csv("res_KO_final_tern.csv")
View(res_KO_final_tern)

tern.1 <- merge(tern, res_KO_final_tern, by ="feature",all=T)
library(tidyr)
tern.1 %>% drop_na(KO_Old) -> tern.2
tern.2$tern <- tern.2$tern %>% replace_na('Others')
tern.2$tern <- factor(tern.2$tern, levels = c("KO-Old enriched gene families","KO-Old depleted gene families","Others"))
levels(tern.2$tern)
my.colors2 <- c("KO-Old depleted gene families"="#0e72cc", "KO-Old enriched gene families"="#610081", "Others"="#ffbf00")
my.shapes2 <- c("KO-Old depleted gene families"=17, "KO-Old enriched gene families"=18, "Others"=20)
my.sizes2 <- c("KO-Old depleted gene families"=8, "KO-Old enriched gene families"=8, "Others"=4)


p <- ggtern(data=tern.2, aes(KO_Old, WT_Old, WT_Young)) + 
  geom_point(aes(fill  = tern, color = tern, shape=tern), 
             alpha=0.6, size=8,
             show.legend=T) + scale_color_manual(values = my.colors2) + 
  scale_shape_manual(values = my.shapes2 )+ 
  scale_size_manual(values = my.sizes2) +
  theme(panel.border = element_rect(color = "black"),
        panel.grid =  element_line(color = "black",size=1)) +
  theme_showarrows() +  
  theme(plot.title = element_text(size = 12, color = "black"),
        axis.text.x=element_text(size=12, color = "black"),
        axis.title.x=element_text(size=12, color = "black"),
        axis.title.y=element_text(size=12, color = "black"),
        axis.text.y=element_text(size=12,  color = "black"),
        legend.title=element_blank(), 
        legend.text=element_text(size=12, color = "black"),
        legend.position="top",)

p
pdf(file.path(path, "tern_genefamily_KO-Old_new.pdf"), width = 8, height = 6)
p
dev.off()