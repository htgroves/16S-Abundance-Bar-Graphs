library(Biostrings) 
library(ggplot2)
library(plyr)
library(stringr)
library(gridExtra)
library(reshape2)
library(knitr)
library(vegan)
library(DESeq2)
library(phyloseq)
library(scales)
library(ggthemes)
library(xlsx)
library(rJava)

setwd("M:/2016/Workspaces/M_final_files_for_r_HG0115")

load("HG0115_bargraphs") # this is unrarefied data

total = median(sample_sums(RSV_a))

standf = function(x, t=total) round(t * (x / sum(x)))      

RSV_a = transform_sample_counts(RSV_a, standf)

mRSV_a <- merge_samples(RSV_a, "DatePostInfection")
sample_data(mRSV_a)$DatePostInfection <- factor(sample_names(mRSV_a))

mRSV_a_fam = tax_glom(mRSV_a, "Family")

mRSV_a_fam_rel_abun  = transform_sample_counts(mRSV_a_fam, function(x) x / sum(x) ) 

family.sum = tapply(taxa_sums(mRSV_a_fam_rel_abun), tax_table(mRSV_a_fam_rel_abun)[, "Family"], sum, na.rm=TRUE)
top9fam = names(sort(family.sum, TRUE))[1:8]
mRSV_a_fam_rel_abun_top9 = prune_taxa((tax_table(mRSV_a_fam_rel_abun)[, "Family"] %in% top9fam), mRSV_a_fam_rel_abun)

mRSV_a_fam_rel_abun_top9_per = transform_sample_counts(mRSV_a_fam_rel_abun_top9, function(x) 100 * x/sum(x))

mRSV_aphy = tax_glom(mRSV_a, "Phylum") # merging OTUs at the phylum level

mRSV_aphy_rel_abun  = transform_sample_counts(mRSV_aphy, function(x) x / sum(x) )  # transforming to relative abundance (before removing the less abundant OTUs - makes very little diff at phyla level

phylum.sum = tapply(taxa_sums(mRSV_aphy_rel_abun), tax_table(mRSV_aphy_rel_abun)[, "Phylum"], sum, na.rm=TRUE)
top2phylamRSV_a = names(sort(phylum.sum, TRUE))[1:2]
mRSV_aphy_rel_abun_top2 = prune_taxa((tax_table(mRSV_aphy_rel_abun)[, "Phylum"] %in% top2phylamRSV_a), mRSV_aphy_rel_abun)

mRSV_a_fam_rel_abun_top9_per_phy_top2 = prune_taxa((tax_table(mRSV_a_fam_rel_abun_top9_per)[, "Phylum"] %in% top2phylamRSV_a), mRSV_a_fam_rel_abun_top9_per)

fam_top2phyla_RSV = psmelt(mRSV_a_fam_rel_abun_top9_per_phy_top2)

family_factor = fam_top2phyla_RSV %>%
  mutate(Family=factor(Family, levels=c("vadinBB60", "Ruminococcaceae", "Lactobacillaceae", "Lachnospiraceae", "Bacteroidaceae", "Rikenellaceae", "Porphyromonadaceae", "S24_7"), ordered=TRUE))

ggplot(family_factor, aes(x=DatePostInfection, y=Abundance)) +
  geom_bar(aes(fill=Family), stat="identity", colour = "black") +
  scale_x_discrete(limits=c("D0", "D4", "D7")) +
  xlab("\nDay After Infection") +
  ylab("% relative abundance\n") +
  ylim(0, 100) +
  ggtitle("RSV") +
  facet_grid(. ~ Phylum) +
  scale_fill_manual(values=c("#CC79A7", "#E69F00", "#999999", "#009E73", "#F0E442", "aquamarine", "#D55E00", "cadetblue3")) +
  theme_base() +
  theme(strip.text.x = element_text(size=18)) +
  theme(plot.title = element_text(face = "bold", size = 25)) + 
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  theme(plot.title = element_text(margin=margin(0,0,15,0))) +
  theme(axis.title.x = element_text(margin=margin(0,5,0,0))) +
  theme(axis.title.y = element_text(margin=margin(5,0,0,0))) +
  theme(plot.title = element_text(hjust=0.5))  +
  theme(plot.margin = unit(c(0,0,0,1), "cm")) +
  theme(plot.background=element_blank())
