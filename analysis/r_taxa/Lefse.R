#set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R")
# create plot folder
path <- "./plots"
#load phyloseq object
phyloseq_obj1 <- readRDS("./Robjects/myphyloseq1.rds")
phyloseq_obj2 <- readRDS("./Robjects/myphyloseq2.rds")
phyloseq_obj3 <- readRDS("./Robjects/myphyloseq3.rds")
phyloseq_obj4 <- readRDS("./Robjects/myphyloseq4.rds")
phyloseq_obj5 <- readRDS("./Robjects/myphyloseq5.rds")
#load required library
library(microbiomeMarker)
####################################
my.themes <- theme(axis.text.x=element_text(size=12, color="black"),
                                    axis.title.x=element_text(size=12, color="black"),
                                    axis.text.y=element_text(size=12, color="black"),
                                    axis.title.y=element_text(size=10, color="black"),
                                    legend.position="top",
                                    legend.text=element_text(size=12, color="black"),
                                    legend.title = element_blank(),) +
                   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"))
tax_level = "Species"
group = "Group"
########################################################################
#WT young vs WT old
phylo_object <- phyloseq_obj1
sample_data(phylo_object)$Group <- factor(sample_data(phylo_object)$Group, 
                                            levels = c("WT-Old","WT-Young"))
tax_glom <- tax_glom(phylo_object, taxrank = tax_level)
tax_glom@tax_table <- subset(tax_glom@tax_table, select = -c(Kingdom, SGB))
#######

res_lefse1 <- run_lefse(
  tax_glom, group = group, wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 4)
(scipen=999)
res_lefse1 <- as.data.frame(res_lefse1@marker_table)
res_lefse1 <- res_lefse1
res_lefse1$compare <- c("WT")

p1 <- plot_ef_bar(res_lefse1) + scale_fill_manual(values = my.colors) + my.themes
p1
pdf(file.path(path, "lefse_WT_groups.pdf"), width = 5, height =6)
p1
dev.off()
####################################
#oral KO
phylo_object <- phyloseq_obj2
sample_data(phylo_object)$Group <- factor(sample_data(phylo_object)$Group, 
                                          levels = c("KO-Old","KO-Young"))
tax_glom <- tax_glom(phylo_object, taxrank = tax_level)
tax_glom@tax_table <- subset(tax_glom@tax_table, select = -c(Kingdom, SGB))
#######
res_lefse2 <- run_lefse(
  tax_glom, group = group, wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 4)
res_lefse2 <- as.data.frame(res_lefse2@marker_table)
res_lefse2<- res_lefse2
res_lefse2$compare <- c("KO")
p2 <- plot_ef_bar(res_lefse2) + scale_fill_manual(values = my.colors) + my.themes
p2
pdf(file.path(path, "lefse_KO_groups.pdf"), width = 5, height =4)
p2
dev.off()
####################################
#oral Young
phylo_object <- phyloseq_obj3
sample_data(phylo_object)$Group <- factor(sample_data(phylo_object)$Group, 
                                          levels = c("KO-Young","WT-Young"))
tax_glom <- tax_glom(phylo_object, taxrank = tax_level)
tax_glom@tax_table <- subset(tax_glom@tax_table, select = -c(Kingdom, SGB))
#######
res_lefse3 <- run_lefse(
  tax_glom, group = group, wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 4)
#Warning message:
#  No marker was identified 
####################################
#oral Old
phylo_object <- phyloseq_obj4
sample_data(phylo_object)$Group <- factor(sample_data(phylo_object)$Group, 
                                          levels = c("KO-Old","WT-Old"))
tax_glom <- tax_glom(phylo_object, taxrank = tax_level)
tax_glom@tax_table <- subset(tax_glom@tax_table, select = -c(Kingdom, SGB))
#######
res_lefse4 <- run_lefse(
  tax_glom, group = group, wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 4)
res_lefse4 <- as.data.frame(res_lefse4@marker_table)
res_lefse4 <- res_lefse4
res_lefse4$compare <- c("OLD")
p4 <- plot_ef_bar(res_lefse4) + scale_fill_manual(values = my.colors) + my.themes
p4
pdf(file.path(path, "lefse_Old_groups.pdf"), width = 5, height =2)
p4
dev.off()
####################################
res_lefse_all <- rbind(res_lefse1,res_lefse2,res_lefse4)
#plot
################
res_lefse_all$feature <- sub(".*\\|", "\\1", res_lefse_all$feature)
res_lefse_all$compare <- factor(res_lefse_all$compare, levels = c("WT","KO","OLD"))
my.colors <- c("WT-Young"="#f39a8f", "WT-Old"="#d95475", "KO-Young"="#f2d58c", "KO-Old"="#ef811e")
#df_summary$unique_group_code <- as.character(df_summary$unique_group_code)
################
p <- ggplot(res_lefse_all, aes(x=ef_lda, y=reorder(feature,-ef_lda), fill=enrich_group,size=-log(pvalue))) + 
 # geom_dotplot(position = position_dodge2(width = 0.75),outlier.shape=NA) +  
  geom_point(position = position_dodge(0.75), shape = 21) +
  scale_fill_manual(values = my.colors) +  facet_wrap(~compare,  ncol=1, strip.position = "right",scales = "free")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=10, color="black"),
        axis.title.x=element_text(size=10, color="black"),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        plot.title = element_text(hjust = 0.5, size=12, face = "bold"),
        legend.position="top",
        legend.text=element_text(size=10, color="black"),
        legend.title = element_text(size=12, color="black"),
        plot.background = element_rect(fill = "transparent",colour = NA),
  )

p

pdf(file.path(path, "lefse_all_groups_.pdf"), width = 7, height =8)
p
dev.off()