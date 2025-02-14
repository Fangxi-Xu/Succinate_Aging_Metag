#set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")
#load required library
library("phyloseq")
library("igraph")
library("ggplot2")
library("tidyverse")
# create plot folder
path <- "./plots"
#load phyloseq object
phylo <- readRDS("./Robjects/myphyloseq.rds")
physeq_filtered <- subset_samples(phylo, Group %in% c("WT-Young", "WT-Old"))

##############################################
#assign colors
my.colors <- c("WT-Young"="#f39a8f", "WT-Old"="#d95475")
################################################
# Define the colors for the groups
group_colors <- c("WT-Young" = "#f39a8f", "WT-Old" = "#d95475", "unknown" = "grey")

# Extract sample data and convert to a data frame
sample_data_df <- as.data.frame(sample_data(physeq_filtered))

# Create a data frame of edges
edge_list <- as.data.frame(as.table(as.matrix(otu_table)))
edge_list <- edge_list[edge_list$Freq > 0, ] # Keep only non-zero edges
colnames(edge_list) <- c("species", "sample", "abundance")

# Merge edge list with sample data to get the group for each edge
edge_list <- merge(edge_list, sample_data_df, by.x = "sample", by.y = "row.names")

# Assign colors based on group
edge_list$color <- group_colors[edge_list$Group]

# Create a graph object
network <- graph_from_data_frame(d = edge_list[, c("species", "sample")], directed = FALSE)

# Add attributes to the graph object
E(network)$abundance <- edge_list$abundance
E(network)$color <- edge_list$color
V(network)$type <- ifelse(V(network)$name %in% rownames(sample_data_df), "sample", "species")
V(network)$group <- ifelse(V(network)$type == "sample", sample_data_df$Group[match(V(network)$name, rownames(sample_data_df))], NA)

# Ensure no NA in type attribute
V(network)$type[is.na(V(network)$type)] <- "species"

# Convert type to a factor
V(network)$type <- as.factor(V(network)$type)

# Ensure no NA in group attribute for plotting
V(network)$group[is.na(V(network)$group)] <- "unknown"

# Plot the bipartite network using ggraph
# Plot the network using ggraph without specifying the bipartite layout
edge_width_scale <- 0.01
edge_alpha <- 0.5
ggraph(network,layout = "fr") + 
  geom_edge_link(aes(width = abundance * edge_width_scale), color = color, alpha = edge_alpha) +
  geom_node_point(aes(color = factor(group), shape = type), size = 3) +
  scale_color_manual(values = c("#f39a8f", "#d95475", "grey"), na.value = "grey") +
  theme_void() +
  ggtitle("Sample-Species Network with Group A (red), Group B (blue), and Species (green)")

library(ggraph)

# Adjust edge width scaling factor
edge_width_scale <- 0.1
edge_alpha <- 0.5

# Plot with adjusted edge aesthetics and consistent edge colors
p <- ggraph(network, layout = "fr") + 
  geom_edge_link(aes(color = I(E(network)$color), width = abundance * edge_width_scale), alpha = edge_alpha) +
  geom_node_point(aes(color = factor(group), shape = type), size = 3,stroke = 0.5) +
  scale_color_manual(values = c("#CD8E87", "#B6556B","darkgrey"), na.value = "grey") +
  theme_void() +
  ggtitle("Sample-Species Network")

pdf(file.path(path, "Sample-Species Network WT.pdf"),width=8,height=6)
p
dev.off()

##############################################
#load required library
library(phyloseq)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
##############################################
# Subset the data for WT-Young and WT-Old
#taxa = c("Species")
group = c("Group")
phylo_melt <- tax_glom(physeq_filtered, taxa) %>%#agglomerate at the domain taxonomic rank 
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

#prepare column annotation - group
col_ann <-  as.data.frame(t(phylo_reshape))
library(tibble)
col_ann_1 <- tibble::rownames_to_column(col_ann, "Group")
row.names(col_ann_1) <- row.names(col_ann)
col_ann_1 <- col_ann_1[,1, drop=FALSE]

df <- as.data.frame(t(phylo_reshape))
df_ordered <- df[rownames(col_ann), ]#reorder matrix based on annotation

df_zscore <- log2(df_ordered +1)#Z score transformation
summary(df_zscore)

library(RColorBrewer)
myColor <- rev(brewer.pal(11, "Spectral"))
#myColor <- colorRamp2(c(0, 0.01,0.1), hcl_palette = "Blue-Red 2")
ann_colors = list(Group = c("WT-Young"="#f39a8f", "WT-Old"="#d95475"))#specify colors for annotation
p<-pheatmap(t(df_zscore), color = myColor,
            cellheight = 12,
            cellwidth = 25,
            fontsize_row=8, 
            fontsize_col=10, 
            show_colnames = F,show_rownames = T, annotation_col = col_ann_1, annotation_colors = ann_colors,
            cluster_cols = F, border_color = NA)

p
pdf(file.path(path,"WT_species_heatmap_group_avg_new.pdf"), width = 8, height =12)
p
dev.off()


# Save the heatmap to a PDF
pdf(file.path(path, "species_heatmap_individual_samples.pdf"), width = 16, height = 8)
print(p)
dev.off()



