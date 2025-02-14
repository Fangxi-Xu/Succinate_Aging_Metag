# Set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")

# Load required libraries
library(phyloseq)
library(igraph)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggraph)

# Create plot folder
path <- "./plots"

# Load phyloseq object
phylo <- readRDS("./Robjects/myphyloseq.rds")
physeq_filtered <- subset_samples(phylo, Group %in% c("WT-Young", "WT-Old"))

# Assign colors
group_colors <- c("WT-Young" = "#f39a8f", "WT-Old" = "#d95475", "unknown" = "grey")

# Extract sample data and convert to a data frame
sample_data_df <- as.data.frame(sample_data(physeq_filtered))

# Create and process edge list for network graph
edge_list <- as.data.frame(as.table(as.matrix(otu_table(physeq_filtered))))
edge_list <- subset(edge_list, Freq > 0)
colnames(edge_list) <- c("species", "sample", "abundance")
edge_list <- merge(edge_list, sample_data_df, by.x = "sample", by.y = "row.names")
edge_list$color <- group_colors[edge_list$Group]

# Create and plot network graph
network <- graph_from_data_frame(d = edge_list[, c("species", "sample")], directed = FALSE)
E(network)$abundance <- edge_list$abundance
E(network)$color <- edge_list$color
V(network)$type <- ifelse(V(network)$name %in% rownames(sample_data_df), "sample", "species")
V(network)$group <- ifelse(V(network)$type == "sample", sample_data_df$Group[match(V(network)$name, rownames(sample_data_df))], "unknown")
V(network)$type <- as.factor(V(network)$type)

g <- ggraph(network, layout = "fr") + 
  geom_edge_link(aes(color = I(E(network)$color), width = abundance * 0.1), alpha = 0.5) +
  geom_node_point(aes(color = factor(group), shape = type), size = 3, stroke = 0.5) +
  scale_color_manual(values = group_colors, na.value = "grey") +
  theme_void() +
  ggtitle("Sample-Species Network")

pdf(file.path(path, "Sample-Species Network WT.pdf"), width = 8, height = 6)
print(g)
dev.off()
