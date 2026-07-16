# set to OWN working directory
rm(list=ls())
setwd("/scratch/project_2007483/FINAL_material/microbio")

.libPaths(c("/projappl/pwatts/rpackages_440", .libPaths()))
libpath <- .libPaths()[1]
.libPaths() 

# setup environment
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)
library(bitops)
library(ggpubr)
library(hrbrthemes)
library(tidyverse)
library(devtools)
library(knitr)
library(rstatix)
library(ggsci)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(microbiome)
library(decontam)
library(ANCOMBC) #v2.9.1
library(metacoder)
library(ggpubr)
library(rlang)
library(tibble)
library(matrixStats)
library(lmerTest)
library(foreach)
library(modeest)
library(pheatmap)
library(tidyr)
library(stringr)
library(rmarkdown)
library(pander)
library(tinytex)
library(scales)
library(ape)
library(ggrepel)
library(parallel)
library(aplot)
library(yaml)
library(biomformat)
library(Biostrings)
library(phangorn) 
library(gcookbook)
library(kableExtra)
library(breakaway)
library(gcookbook)
library(ggh4x)
library(qiime2R)

library(gghalves)
library(microbiomeutilities)

#####DATA ANALYSES
#####These data have been exported from QIIME2
#Read in OTU tableotu_table <- read.csv("otu_tableFINAL.txt", sep = "\t", row.names = 1)
otu_table <- read.csv("otu_tableFINAL.txt", sep = "\t", row.names = 1)
otu_table <- as.matrix(otu_table)
otu_table

# read in taxonomy
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", header = TRUE, row.names = 1)
taxonomy <- as.matrix(taxonomy)

# read in metadatametadata 
metadata <- read.table("FINALFINAL_microbiometa.txt", sep = "\t", header = TRUE, row.names = 1)


# import as phyloseq objects
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
MET <- sample_data(metadata)

# sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)

sample_names(OTU)
sample_names(MET)

# merge the objects to create a phyloseq object
ps <- phyloseq(OTU, TAX, MET)

## sanity check - do you have the correct amount of data?
ps

#SAVE and skip ps production
save.image(file = "ps.RData")
load(file = "ps.RData")

#ps:n pilkkominen lajeihin
psPF <- subset_samples(ps, species == "PF")
psGTa <- subset_samples(ps, lifestage == "adult")
psGT <- subset_samples(ps, species == "GT")
psGTN <- subset_samples(psGT, lifestage == "nestling")

###########################################################
## ANCOM-BC2
############################################# PHILL

psGTN.genus <- tax_glom(psGTN, taxrank = 'genus', NArm = FALSE)
psGTa.genus <- tax_glom(psGTa, taxrank = 'genus', NArm = FALSE)
psPF.genus <- tax_glom(psPF, taxrank = 'genus', NArm = FALSE)

#factors
sample_data(psGTN.genus)$zone <- factor(sample_data(psGTN.genus)$zone,
                                             levels = c("C", "H"))
sample_data(psGTa.genus)$zone <- factor(sample_data(psGTa.genus)$zone,
                                        levels = c("C", "H"))
sample_data(psPF.genus)$zone <- factor(sample_data(psPF.genus)$zone,
                                        levels = c("C", "H"))

#clean the taxonomy
tax.2 <- data.frame(tax_table(psGTN.genus))
tax.2 <- data.frame(tax_table(psGTa.genus))
tax.2 <- data.frame(tax_table(psPF.genus))
tax.clean.2 <- data.frame(row.names = row.names(tax.2),
                          kingdom = str_replace(tax.2[,1], "unculture.*",""),
                          phylum = str_replace(tax.2[,2], "unculture.*",""),
                          class = str_replace(tax.2[,3], "unculture.*",""),
                          order = str_replace(tax.2[,4], "unculture.*",""),
                          family = str_replace(tax.2[,5], "unculture.*",""),
                          genus = str_replace(tax.2[,6], "unculture.*",""),
                          species = str_replace(tax.2[,7], "unculture.*",""),
                          stringsAsFactors = FALSE)

tax.clean.2[is.na(tax.clean.2)] <- ""
# all uncultured text removed and NAs -> empty cells
for (i in 1:7){ tax.clean.2[,i] <- as.character(tax.clean.2[,i])}

#FILL HOLES
tax.clean.2[is.na(tax.clean.2)] <- ""
for (i in 1:nrow(tax.clean.2)){
  #  fill in missing taxonomy with last assigned taxonomic level
  if (tax.clean.2[i,2] == ""){
    kingdom <- paste("k_", tax.clean.2[i,1], sep = "")
    tax.clean.2[i, 2:7] <- kingdom
  } else if (tax.clean.2[i,3] == ""){
    phylum <- paste("p_", tax.clean.2[i,2], sep = "")
    tax.clean.2[i, 3:7] <- phylum
  } else if (tax.clean.2[i,4] == ""){
    class <- paste("c_", tax.clean.2[i,3], sep = "")
    tax.clean.2[i, 4:7] <- class
  } else if (tax.clean.2[i,5] == ""){
    order <- paste("o_", tax.clean.2[i,4], sep = "")
    tax.clean.2[i, 5:7] <- order
  } else if (tax.clean.2[i,6] == ""){
    family <- paste("f_", tax.clean.2[i,5], sep = "")
    tax.clean.2[i, 6:7] <- family
  } else if (tax.clean.2[i,7] == ""){
    tax.clean.2$Species[i] <- paste("g_",tax.clean.2$Genus[i], sep = "_")
  }
}

# insert the cleaned taxonomy table back to phyloseq object
tax_table(psGTN.genus) <- as.matrix(tax.clean.2)
tax_table(psGTa.genus) <- as.matrix(tax.clean.2)
tax_table(psPF.genus) <- as.matrix(tax.clean.2)

#get reads per sample
# get reads per sample
sample_sums(psGTN.genus)
sort(sample_sums(psGTN.genus))
hist(sample_sums(psGTN.genus), main = NULL, xlab="total reads (N)", 
     border="black", col="#C1000EFF", las=1, breaks=20)

#sample_sums(psGTa.genus)
#sort(sample_sums(psGTa.genus))
#hist(sample_sums(psGTa.genus), main = NULL, xlab="total reads (adults)", 
#     border="black", col="#C1000EFF", las=1, breaks=20)

#SAVE
save.image(file = "ps.cleanedtaxancombc.RData")
load(file = "ps.cleanedtaxancombc.RData")

output.N = ancombc2(data = psGTN.genus,
                  assay_name = "counts",
                  tax_level = "family",
                  fix_formula = "zone",
                  rand_formula = NULL,
                  p_adj_method = "fdr",
                  pseudo = FALSE, #do not adds pseudo count
                  pseudo_sens = TRUE,
                  prv_cut = 0.1, #increase if too man NAs
                  lib_cut = 5900, # default value 1000 to include libraries lose 6 samples<<<<<<<<<
                  group = "zone", #discrete group - set to NULL if only two categories
                  s0_perc = 0.05,
                  struc_zero = TRUE,
                  neg_lb = FALSE,
                  alpha = 0.05,
                  n_cl = 4, #parallel clusters
                  verbose = TRUE,
                  global = FALSE,
                  pairwise = FALSE,
                  dunnet = FALSE,
                  trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20,
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL,
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = NULL)

res_prim.N = output.N$res
res_prim.N
library(WriteXLS)
WriteXLS("res_prim.N","ANCON_GTN_FAMILY.xlsx")

## waterfall for continuous
## "#71A7C4", "#590007"
df_N = res_prim.N %>%
  dplyr::select(taxon, ends_with("zoneH"))

df_fig_N = df_N %>%
  dplyr::filter(diff_zoneH == 1) %>% 
  dplyr::arrange(desc(lfc_zoneH)) %>%
  dplyr::mutate(direct = ifelse(lfc_zoneH > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(diff_zoneH, "black", "grey"))
df_fig_N$taxon = factor(df_fig_N$taxon, levels = df_fig_N$taxon)
df_fig_N$direct = factor(df_fig_N$direct, 
                           levels = c("Positive LFC", "Negative LFC"))
#colours
#"#001260" "#71A7C4" "#D29773" "#590007"
#scale_color_manual(values = c(GTnC="#365C8DFF",GTnH="#FAC127FF")) +
# Plot
fig_N = df_fig_N %>%
  ggplot(aes(x = taxon, y = lfc_zoneH, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(
    ymin = lfc_zoneH - se_zoneH,
    ymax = lfc_zoneH + se_zoneH),
    width = 0.2, color = "black") +
  scale_fill_manual(values = c("Positive LFC" = "#FAC127FF",   
                               "Negative LFC" = "#365C8DFF")) + 
  labs(x = NULL, y = "log fold change\n") +
#  coord_cartesian(ylim = c(-0.05, 0.03)) +  # y-axis limits here
  theme_classic(base_size = 14) +  # Removes background gridlines
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size=14,
                               color = df_fig_N$color),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.text.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

fig_N

#### ORDER
output.N.order = ancombc2(data = psGTa.genus,
                    assay_name = "counts",
                    tax_level = "order",
                    fix_formula = "zone",
                    rand_formula = NULL,
                    p_adj_method = "fdr",
                    pseudo = FALSE, #do not adds pseudo count
                    pseudo_sens = TRUE,
                    prv_cut = 0.1, #increase if too man NAs
                    lib_cut = 5900, # default value 1000 to include libraries lose 6 samples<<<<<<<<<
                    group = "zone", #discrete group - set to NULL if only two categories
                    s0_perc = 0.05,
                    struc_zero = TRUE,
                    neg_lb = FALSE,
                    alpha = 0.05,
                    n_cl = 4, #parallel clusters
                    verbose = TRUE,
                    global = FALSE,
                    pairwise = FALSE,
                    dunnet = FALSE,
                    trend = FALSE,
                    iter_control = list(tol = 1e-5, max_iter = 20,
                                        verbose = FALSE),
                    em_control = list(tol = 1e-5, max_iter = 100),
                    lme_control = NULL,
                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                    trend_control = NULL)

res_prim.N.order = output.N.order$res
res_prim.N

WriteXLS("res_prim.N","ANCON_GTN_ORDER.xlsx")
#SAVE
save.image(file = "res_prim.N.GTNorder")
load(file = "res_prim.N.GTNorder")

## waterfall for continuous
## "#71A7C4", "#590007"
df_N.order = res_prim.N.order %>%
  dplyr::select(taxon, ends_with("zoneH"))

df_fig_N.order = df_N.order %>%
  dplyr::filter(diff_zoneH == 1) %>% 
  dplyr::arrange(desc(lfc_zoneH)) %>%
  dplyr::mutate(direct = ifelse(lfc_zoneH > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(diff_zoneH, "black", "grey"))
df_fig_N.order$taxon = factor(df_fig_N.order$taxon, levels = df_fig_N.order$taxon)
df_fig_N.order$direct = factor(df_fig_N.order$direct, 
                         levels = c("Positive LFC", "Negative LFC"))
#colours
#"#001260" "#71A7C4" "#D29773" "#590007"

# Plot
fig_N.order = df_fig_N.order %>%
  ggplot(aes(x = taxon, y = lfc_zoneH, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(
    ymin = lfc_zoneH - se_zoneH,
    ymax = lfc_zoneH + se_zoneH),
    width = 0.2, color = "black") +
  scale_fill_manual(values = c("Positive LFC" = "#FAC127FF",   
                               "Negative LFC" = "#365C8DFF")) +
  labs(x = NULL, y = "log fold change\n") +
  #  coord_cartesian(ylim = c(-0.05, 0.03)) +  # y-axis limits here
  theme_classic(base_size = 14) +  # Removes background gridlines
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size=18,
                               color = df_fig_N.order$color),
    axis.title.y = element_text(size = 18, face = "plain"),
    axis.text.y = element_text(size = 18),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

fig_N.order

#SAVE
save.image(file = "ps.ancombc2.1.RData")
#load(file = "ps.ancombc2.1.RData")




################################>>>>>>>>>>>>>>>>PF
#### order PF
output.N.order = ancombc2(data = psPF,
                           assay_name = "counts",
                           tax_level = "order",
                           fix_formula = "zone",
                           rand_formula = NULL,
                           p_adj_method = "fdr",
                           pseudo = FALSE, #do not adds pseudo count
                           pseudo_sens = TRUE,
                           prv_cut = 0.1, #increase if too man NAs
                           lib_cut = 5900, # default value 1000 to include libraries lose 6 samples<<<<<<<<<
                           group = "zone", #discrete group - set to NULL if only two categories
                           s0_perc = 0.05,
                           struc_zero = TRUE,
                           neg_lb = FALSE,
                           alpha = 0.05,
                           n_cl = 4, #parallel clusters
                           verbose = TRUE,
                           global = FALSE,
                           pairwise = FALSE,
                           dunnet = FALSE,
                           trend = FALSE,
                           iter_control = list(tol = 1e-5, max_iter = 20,
                                               verbose = FALSE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = NULL,
                           mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                           trend_control = NULL)

res_prim.PFN.order = output.N.order$res
res_prim.PFN.order

## waterfall for continuous
## "#71A7C4", "#590007"
df_N.phylum = res_prim.N.phylum %>%
  dplyr::select(taxon, ends_with("zoneH"))

df_fig_N.phylum = df_N.phylum %>%
  dplyr::filter(diff_zoneH == 1) %>% 
  dplyr::arrange(desc(lfc_zoneH)) %>%
  dplyr::mutate(direct = ifelse(lfc_zoneH > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(diff_zoneH, "black", "grey"))
df_fig_N.phylum$taxon = factor(df_fig_N.phylum$taxon, levels = df_fig_N.phylum$taxon)
df_fig_N.phylum$direct = factor(df_fig_N.phylum$direct, 
                                levels = c("Positive LFC", "Negative LFC"))
#colours
#"#001260" "#71A7C4" "#D29773" "#590007"

# Plot
fig_N.phylum = df_fig_N.phylum %>%
  ggplot(aes(x = taxon, y = lfc_zoneH, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(
    ymin = lfc_zoneH - se_zoneH,
    ymax = lfc_zoneH + se_zoneH),
    width = 0.2, color = "black") +
  scale_fill_manual(values = c("Positive LFC" = "#FAC127FF",   
                               "Negative LFC" = "#365C8DFF")) + 
  labs(x = NULL, y = "log fold change\n") +
  #  coord_cartesian(ylim = c(-0.05, 0.03)) +  # y-axis limits here
  theme_classic(base_size = 14) +  # Removes background gridlines
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size=18,
                               color = df_fig_N.phylum$color),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.text.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

fig_N.phylum

#SAVE
save.image(file = "res_prim.PFN.order")
#load(file = "ps.ancombc2.1.RData")


