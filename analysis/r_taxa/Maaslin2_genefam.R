#set working directory
setwd("/Volumes/T7/Analysis/aging_manuscript_oral/R-genefamily")
# create plot folder
#dir.create("./plots")
path <- "./plots"
# create results folder
#dir.create("./Maaslin2_res")
#Load MaAsLin2 package into the R environment
library("Maaslin2")
library(phyloseq)
library(dplyr)
library(tidyverse)
#?Maaslin2
#load phyloseq object
phylo <- readRDS("./Robjects/myphyloseq.rds")
################################################################################
#ternary plot needs to be colored based on masslin2 results
#so perform masslin2 first
taxa = c("Gene")

# create dataframe from phyloseq object 
phylo_melt <- tax_glom(phylo, taxa) %>%#agglomerate at the domain taxonomic rank 
  psmelt()# melt phylo object into large data frame

df = subset(phylo_melt, select = c("Sample",taxa,"Abundance"))#select columns to use; put value variable at the end
#merge duplicated rows
df %>%
  group_by(Sample,.data[[taxa]]) %>%
  summarise_all(sum) %>%
  data.frame() -> df_agg
#transform data
df_reshape <- reshape2::dcast(df_agg, Sample ~ Gene, value.var='Abundance') 
df_reshape <- df_reshape %>% remove_rownames %>% column_to_rownames(var="Sample")
df_reshape <- as.data.frame(t(df_reshape))
#prepare metadata for comparison
meta <- data.frame(sample_data(phylo))

fit_data = Maaslin2(
  input_data = df_reshape, 
  input_metadata = meta, 
  output = "./Maaslin2_res", 
  fixed_effects = c("Group"),
  reference = c("Group,WT-Young"),
  max_significance = 0.1)

saveRDS(fit_data, "./Robjects/fit_data_maaslin2.rds")

################################################################################
#extract significant results
sig_res <- as.data.frame(filter(fit_data$results, qval < 0.1))
#extract WT Old enriched and depleted gene families
sig_WT <- subset(sig_res, value == "WT-Old")

sig_WT_plus <- subset(sig_WT, sig_WT$coef > 0)
sig_WT_plus$tern1 <- "WT_Old enriched gene families"

sig_WT_minus <- subset(sig_WT, sig_WT$coef < 0)
sig_WT_minus$tern1 <- "WT_Old depleted gene families"

res_WT_final_tern <- rbind(sig_WT_minus,sig_WT_plus)

write.csv(res_WT_final_tern, file = "res_WT_final_tern.csv")

#extract KO Old enriched and depleted gene families
sig_KO <- subset(sig_res, value == "KO-Old")

sig_KO_plus <- subset(sig_KO, sig_KO$coef > 0)
sig_KO_plus$tern1 <- "KO-Old enriched gene families"

sig_KO_minus <- subset(sig_KO, sig_KO$coef < 0)
sig_KO_minus$tern1 <- "KO-Old depleted gene families"

res_KO_final_tern <- rbind(sig_KO_minus,sig_KO_plus)

write.csv(res_KO_final_tern, file = "res_KO_final_tern.csv")
################################################################################
