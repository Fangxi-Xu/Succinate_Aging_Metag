#set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")
# create plot folder
#dir.create("./plots")
#create sub folders for different analysis
path <- "./plots"
#load phyloseq object
phylo <- readRDS("./Robjects/myphyloseq.rds")
##############################################
#load required library
library(phyloseq)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
##############################################
#species level heatmap (group average based)
taxa = c("Species")
group = c("Group")
phylo_melt <- tax_glom(phylo, taxa) %>%#agglomerate at the domain taxonomic rank 
  psmelt()# melt phylo object into large data frame

df = subset(phylo_melt, select = c(group,taxa, "Abundance"))#select columns to use; put value variable at the end

#calculate group mean
df_summary <- df %>%
  group_by(.data[[taxa]], .data[[group]]) %>%
  dplyr::summarise(mean_Abun = mean(Abundance))
#merge duplicated rows if any
#df_summary %>%
#  group_by(.data[[group]],.data[[taxa]]) %>%
#  dplyr::summarise_all(sum) %>%
#  data.frame() -> df_summary

phylo_reshape <- reshape2::dcast(df_summary, Group ~ Species, value.var='mean_Abun') %>% #transform data 
  remove_rownames %>% column_to_rownames(var="Group") %>% #column "Sample" to row name
  t()%>% as.data.frame()

#prepare row annotation - species
row_ann <- phylo_melt[,c(30:36)]
row_ann <- row_ann[!duplicated(row_ann), ] %>%
  remove_rownames %>% column_to_rownames(var="Species")
row_ann <- row_ann[,2, drop=FALSE]

#prepare column annotation - group
col_ann <-  as.data.frame(t(phylo_reshape))
library(tibble)
col_ann_1 <- tibble::rownames_to_column(col_ann, "Group")
row.names(col_ann_1) <- row.names(col_ann)
col_ann_1 <- col_ann_1[,1, drop=FALSE]


df <- as.data.frame(t(phylo_reshape))
df_ordered <- df[rownames(col_ann), ]#reorder matrix based on annotation

df_zscore <- log2(df_ordered +1)#log2 transformation
summary(df_zscore)

library(RColorBrewer)
myColor <- rev(brewer.pal(11, "Spectral"))
#myColor <- colorRamp2(c(0, 0.01,0.1), hcl_palette = "Blue-Red 2")
ann_colors = list(Group = c("WT-Young"="#f39a8f", "WT-Old"="#d95475", "KO-Young"="#f2d58c", "KO-Old"="#ef811e"),
                  Phylum = c("p__Firmicutes"="#62bf82", "p__Actinobacteria"="#629fbf", "p__Proteobacteria"="#bf629f","p__Bacteroidetes"="#bf8262",
                             "p__Candidatus_Saccharibacteria" = "#9fbf62","p__Verrucomicrobia"="#bfb162"))#specify colors for annotation
p<-pheatmap(t(df_zscore), color = myColor,
            cellheight = 12,
            cellwidth = 25,
            fontsize_row=8, 
            fontsize_col=10, 
            show_colnames = F,show_rownames = T, annotation_col = col_ann_1, annotation_row = row_ann, annotation_colors = ann_colors,
            cluster_cols = F, border_color = NA)

p
pdf(file.path(path,"species_heatmap_group_avg_new.pdf"), width = 8, height =12)
p
dev.off()











