#functions for beta diversity PCoA
#author: Fangxi Xu
#load required packages
library(phyloseq)
library(ggplot2)
library(vegan)

plot.bray.curtis <- function(phyloseq_obj, group) {
  phylo_object <- phyloseq_obj #input data has to be a phyloseq object
  bray<- ordinate(phylo_object, "PCoA", "bray")#
  p<- plot_ordination(phylo_object, 
                       bray, color=group, shape = group) + 
    stat_ellipse(geom = "polygon", type="norm", linewidth=0.5, alpha=0.3, aes(color=.data[[group]], fill=.data[[group]]))+
    scale_fill_manual(values = my.colors, guide = "none")+
    scale_color_manual(values = my.colors)+
    scale_shape_manual(values = group.shapes)+
    ggtitle("Bray-Curtis") + geom_point(size = 3.5)+
    theme_classic() + 
    theme(plot.title = element_text(size = 12, color = "black"),
          axis.text.x=element_text(size=10, color = "black"),
          axis.title.x=element_text(size=12, color = "black"),
          axis.title.y=element_text(size=12, color = "black"),
          axis.text.y=element_text(size=10,  color = "black"),
          legend.title=element_blank(), 
          legend.text=element_text(size=8, color = "black"),
          legend.position="left",
          plot.background = element_rect(fill = "transparent",colour = NA),)}

plot.bray.curtis.pairwise <- function(file, group) {
  phylo_object <- file #input data has to be a phyloseq object
  bray<- ordinate(phylo_object, "PCoA", "bray")#
  p <- plot_ordination(phylo_object, 
                       bray, color=group, shape = group) + 
    stat_ellipse(geom = "polygon", type="norm", linewidth=0.5, alpha=0.3, aes(color=.data[[group]], fill=.data[[group]]))+
    scale_fill_manual(values = my.colors, guide = "none")+
    scale_color_manual(values = my.colors)+
    scale_shape_manual(values = group.shapes)+
    ggtitle("Bray-Curtis") + geom_point(size = 3.5)+
    theme_classic() + 
    theme(plot.title = element_text(size = 12, color = "black"),
          axis.text.x=element_text(size=10, color = "black"),
          axis.title.x=element_text(size=12, color = "black"),
          axis.title.y=element_text(size=12, color = "black"),
          axis.text.y=element_text(size=10,  color = "black"),
          legend.title=element_blank(), 
          legend.text=element_text(size=8, color = "black"),
          legend.position="bottom",
          plot.background = element_rect(fill = "transparent",colour = NA),)}




plot.wuni <- function(file, group) {
  phylo_object <- file #input data has to be a phyloseq object
  wuni <- ordinate(phylo_object, "PCoA", "unifrac", weighted=T)
  p <- plot_ordination(phylo_object, 
                       wuni, color=group,shape = group) + 
    stat_ellipse(geom = "polygon", type="norm", linewidth=0.5, alpha=0.3, aes(color=.data[[group]], fill=.data[[group]]))+
    scale_fill_manual(values = my.colors, guide = "none")+
    scale_color_manual(values = my.colors)+
    scale_shape_manual(values = group.shapes)+
    ggtitle("Weighted UniFrac") + geom_point(size = 3.5)+
    theme_classic() + 
    theme(plot.title = element_text(size = 12, color = "black"),
          axis.text.x=element_text(size=10, color = "black"),
          axis.title.x=element_text(size=12, color = "black"),
          axis.title.y=element_text(size=12, color = "black"),
          axis.text.y=element_text(size=10,  color = "black"),
          legend.title=element_blank(), 
          legend.text=element_text(size=8, color = "black"),
          legend.position="left",
          plot.background = element_rect(fill = "transparent",colour = NA),)}

plot.wuni.pairwise <- function(file, group) {
  phylo_object <- file #input data has to be a phyloseq object
  wuni <- ordinate(phylo_object, "PCoA", "unifrac", weighted=T)
  p <- plot_ordination(phylo_object, 
                       wuni, color=group,shape = group) + 
    stat_ellipse(geom = "polygon", type="norm", linewidth=0.5, alpha=0.3, aes(color=.data[[group]], fill=.data[[group]]))+
    scale_fill_manual(values = my.colors, guide = "none")+
    scale_color_manual(values = my.colors)+
    scale_shape_manual(values = group.shapes)+
    ggtitle("Weighted UniFrac") + geom_point(size = 3.5)+
    theme_classic() + 
    theme(plot.title = element_text(size = 12, color = "black"),
          axis.text.x=element_text(size=10, color = "black"),
          axis.title.x=element_text(size=12, color = "black"),
          axis.title.y=element_text(size=12, color = "black"),
          axis.text.y=element_text(size=10,  color = "black"),
          legend.title=element_blank(), 
          legend.text=element_text(size=8, color = "black"),
          legend.position="bottom",
          plot.background = element_rect(fill = "transparent",colour = NA),)}



test_permanova <- function(file, group) {
  library(vegan)
  if (BC == TRUE){
    dist <- phyloseq::distance(file, method = "bray")
    f <- as.formula(paste0("dist ~ ", group))
    perma <- adonis2(f, data = as(sample_data(file), "data.frame"), permutations = 999)
    perma
  }
  else{
    
    dist <- UniFrac(file, 
                    weighted = T, 
                    normalized = TRUE,  
                    parallel = FALSE, 
                    fast = TRUE)
    f <- as.formula(paste0("dist ~ ", group))
    perma <- adonis2(f, data = as(sample_data(file), "data.frame"), permutations = 999)
    perma
  }
}