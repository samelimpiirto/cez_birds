# CHERNOBYL
rm(list=ls())
# own packages - in project libpath
.libPaths(c("/projappl/project_2007483/rpackages_440", .libPaths())) 
#.libPaths(c("/projappl/pwatts/rpackages_421", .libPaths()))
libpath <- .libPaths()[1]

# set working directory
setwd("/scratch/project_2007483/FINAL_material/microbio")

# setup environment
library(phyloseq) 
library(ggplot2)
library(ggpubr)
library(vegan)
#library(GenomeInfoDb)
#library(GenomeInfoDbData)
library(bitops)
library(dplyr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse) 
library(devtools)
library(knitr)
library(rstatix)
library(ggsci)
library(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
#install.packages("devtools")
#devtools::install_github("adw96/breakaway")
#library(breakaway)
library(microbiome)
library(ANCOMBC)
library(eulerr)
#BiocManager::install("decontam", lib = libpath)
library(decontam)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ANCOMBC", force = TRUE)
#remotes::install_github("R-CoderDotCom/ggcats@main")
#library(ANCOMBC2)
#install.packages("paletteer")
library(paletteer)
library(viridis)           
library(scales)
########################################################################
#show_col(viridis_pal()(20))
## DECONTAM 
## https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
###############################################################################

#inspect library sizes
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=categoryHC)) + geom_point()

#identify contaminants - based on prevalence (frequency needs accurate? quant)
#search in column labelled *sample_control* with control...
sample_data(ps)$is.neg <- sample_data(ps)$species == "neg_cont"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#remember that the number of SVs must equal that in the phyloseq object

head(which(contamdf.prev$contaminant))

#change the prevalence threshold to 0.5 - more *aggressive* detection
#I used this one to decontam
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$species == "neg_cont", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_type == "sample", ps.pa)


# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

## https://rpubs.com/microbiotic/833244
##ordinate samples to see if contaminants are separate
a.ord <- ordinate(ps, "PCoA", "bray")
plot_ordination(ps, a.ord, type="samples", color="sample_type")

#grab the row indices that correspond with identified contaminants to locate taxonomic information in the corresponding OTU file
row_indices <- which(contamdf.prev05$contaminant) 

# prune contaminant taxa using phyloseqs prune_taxa function.
ps.final_biom <- prune_taxa(!contamdf.prev05$contaminant, ps)
ps.final_biom

## get table of SVs - remove these in QIIME 
## examine if contaminants are rare or common SVs

sd = data.frame(sample_data(ps.pa))
sd.neg = data.frame(sample_data(ps.pa.neg))
sd = data.frame(sample_data(ps.pa))

droplist <- contamdf.prev05[which(contamdf.prev05$"contaminant" == "TRUE"),]
drop2 <- subset(droplist, select = c(contaminant))
library(WriteXLS)
WriteXLS("drop2","contASVs.xlsx")

########################################################################
## DATA ANALYSES
## These data have been exported from QIIME2
otu_table <- read.csv("otu_tableFINAL.txt", sep = "\t", row.names = 1)
otu_table <- as.matrix(otu_table)
otu_table

# read in taxonomy
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", header = TRUE, row.names = 1)
taxonomy <- as.matrix(taxonomy)

# read in metadata 
metadata <- read.table("METADATA.txt", sep = "\t", header = TRUE, row.names = 1)

# Read in tree
#tree <- read_tree("tree.nwk")

# import as phyloseq objects
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
MET <- sample_data(metadata)

# sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)
#taxa_names(tree)

sample_names(OTU)
sample_names(MET)

# merge the objects to create a phyloseq object
#pstree <- phyloseq(OTU, TAX, MET, tree)
ps <- phyloseq(OTU, TAX, MET)

# sanity check - do you have the correct amount of data?
ps

# check taxonomy table
phyloseq::tax_table(ps)[1:10, 1:7] # check table and show 1-10 features, and 1-7 ranks

# CLEAN TAXONOMY TABLE
# turn the ps tax_table as a data frame
tax <- data.frame(tax_table(ps))

# remove any leading taxonomy text
tax.clean.1 <- data.frame(row.names = row.names(tax),
                        domain = str_replace(tax[,1], "unculture.*",""),
                        phylum = str_replace(tax[,2], "unculture.*",""),
                        class = str_replace(tax[,3], "unculture.*",""),
                        order = str_replace(tax[,4], "unculture.*",""),
                        family = str_replace(tax[,5], "unculture.*",""),
                        genus = str_replace(tax[,6], "unculture.*",""),
                        species = str_replace(tax[,7], "unculture.*",""),
                        stringsAsFactors = FALSE)


tax.clean <- data.frame(row.names = row.names(tax.clean.1),
                        domain = str_replace(tax.clean.1[,1], "metagenom.*",""),
                        phylum = str_replace(tax.clean.1[,2], "metagenom.*",""),
                        class = str_replace(tax.clean.1[,3], "metagenom.*",""),
                        order = str_replace(tax.clean.1[,4], "metagenom.*",""),
                        family = str_replace(tax.clean.1[,5], "metagenom.*",""),
                        genus = str_replace(tax.clean.1[,6], "metagenom.*",""),
                        species = str_replace(tax.clean.1[,7], "metagenom.*",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
# all NAs are now empty cells
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

# fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  # fill in missing taxonomy with last assigned taxonomic level
  if (tax.clean[i,2] == ""){
    kingdom <- paste("domain", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$species[i] <- paste("genus",tax.clean$genus[i], sep = "_")
  }
}

# insert the cleaned taxonomy table back to phyloseq object
tax_table(ps) <- as.matrix(tax.clean)
ps

# check taxonomy table
phyloseq::tax_table(ps)[1:10, 1:7] # check table and show 1-10 features, and 1-7 ranks

set.seed(111) # keep result reproducable data
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=5000, replace=F)
ps.rarefied
#ps.comp <- transform(ps.rarefied, "compositional")

# SAVE and skip ps production
#save.image(file = "ps.RData")
load(file = "ps.RData")

######################################## 
## work on pooled samples
ps.pool <- merge_samples(ps.rarefied, "categoryHC")
physeq3.p = transform_sample_counts(ps.pool, function(x) x / sum(x) )
physeq3.p

glom.p <- tax_glom(physeq3.p, taxrank = 'phylum')
glom.p # should list # taxa as # phyla
data_glom.p <- psmelt(glom.p) # create dataframe from phyloseq object
data_glom.p$phylum <- as.character(data_glom.p$phylum) #convert to character

#simple way to rename phyla with < 5% abundance
data_glom.p$phylum[data_glom.p$Abundance < 0.05] <- "<5% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom.p$phylum))
Count

unique(data_glom.p$phylum)
# 5% abundance
data_glom.p$phylum <- factor(data_glom.p$phylum, 
                             levels = c("Actinobacteriota",
                                        "Desulfobacterota",
                                        "Firmicutes",
                                        "Proteobacteria",
                                        "<5% abund."))

# "final" version, use this, needs: labels for species, decide what to use for 'treatment'
CHorder<-c("PFnC", "PFnH", "GTnC", "GTnH", "GTaC", "GTaH")
#species.labs <- c("pf nestling", "gt nestling", "gt adult")
#names(species.labs) <- c("PF", "GT", "GT")

##############################################################################################


FIG1a <-  ggplot(data=data_glom.p, aes(x = Sample, y = Abundance, fill = phylum)) +
    geom_bar(stat = "identity", size = 1) +   # border thickness
    labs(x = "\n\n", y = "\nRelative abundance\n") +
  #  scale_fill_viridis_d("phylum", option = "inferno" ) +
    scale_fill_viridis_d("Phylum", option = "viridis" ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 20),
          axis.text.y = element_text(angle = 360, hjust = 0.5, size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) + #kills facet headers (= strip)
    theme(panel.background = element_blank()) +
    scale_x_discrete(limits=CHorder, labels=c("PF nest  - Un",
                                              "PF nest  - Ct",
                                              "GT nest  - Un",
                                              "GT nest  - Ct",
                                              "GT adult - Un",
                                              "GT adult - Ct")) +
    theme(axis.text=element_text(size=20), 
          legend.text = element_text(size=20), 
          legend.title = element_text(size=20)) 

FIG1a

##############################################################################################
## FACET WRAP
data_glom.p$Group <- dplyr::case_when(
  data_glom.p$Sample %in% c("PFnC","PFnH") ~ "PF nest",
  data_glom.p$Sample %in% c("GTnC","GTnH") ~ "GT nest",
  data_glom.p$Sample %in% c("GTaC","GTaH") ~ "GT adult"
)

# safer?
#data_glom.p <- data_glom.p %>%
#  mutate(Group = dplyr::case_when(
#    Sample %in% c("PFnC", "PFnH") ~ "PF nest",
#    Sample %in% c("GTnC", "GTnH") ~ "GT nest",
#    Sample %in% c("GTaC", "GTaH") ~ "GT adult",
#    TRUE ~ "Other"
#  ))

# define custom x-axis labels
x_labels <- c(
  "PFnC" = "Unc",
  "PFnH" = "Con",
  "GTnC" = "Unc",
  "GTnH" = "Con",
  "GTaC" = "Unc",
  "GTaH" = "Con"
)
                                

# custom facet labels
facet_labels <- c(
  "PF nest" = "F. hypoleuca nestlings",
  "GT nest" = "P. major nestlings",
  "GT adult" = "P. major adults"
)

#plot with faceting
FIG1a_facet <- ggplot(data_glom.p, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x",
             labeller = labeller(Group = facet_labels)) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "\n\n", y = "\nRelative abundance\n") +
  theme_bw() +  # start from a clean theme
  theme(
    panel.spacing.x = unit(0.3, "cm"),          # gap between facets
    panel.grid = element_blank(),               # remove grid lines
    panel.border = element_blank(),              # no border
    strip.background = element_rect(
      fill = "white",                       # facet background color
      color = "white",
      size = 0.5
    ),
    strip.text.x = element_text(
      size = 16,
      face = "plain",
      color = "black"
    ),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    panel.background = element_blank()
  )

FIG1a_facet


# ORDER
glom.o <- tax_glom(physeq3.p, taxrank = 'order')
glom.o # should list # taxa as # order
data_glom.o <- psmelt(glom.o) # create dataframe from phyloseq object
data_glom.o$order <- as.character(data_glom.o$order) #convert to character

#simple way to rename phyla with < 5% abundance
data_glom.o$order[data_glom.o$Abundance < 0.05] <- "<5% abund."

#Count # phyla to set color palette
Count.o = length(unique(data_glom.o$order))
Count.o #returns 12

unique(data_glom.o$order)

# 5% abundance
data_glom.o$order <- factor(data_glom.o$order, 
                             levels = c( "Clostridiales" ,     
                                         "Corynebacteriales"  ,
                                         "Desulfovibrionales",
                                         "Diplorickettsiales" ,
                                         "Enterobacterales"  , 
                                         "Frankiales" ,
                                         "Lactobacillales" ,
                                         "Micrococcales"   ,  
                                         "Mycoplasmatales"  , 
                                         "Rhizobiales"  ,      
                                         "Solirubrobacterales",
                                         "<5% abund."))

FIG1b <-  ggplot(data=data_glom.o, aes(x = Sample, y = Abundance, fill = order)) +
  geom_bar(stat = "identity", size = 1) +   # border thickness
  labs(x = "\n\n", y = "\nRelative abundance\n") +
  #  scale_fill_viridis_d("Order", option = "inferno" ) +
  scale_fill_viridis_d("Order", option = "viridis" ) + #makes me dizzy<<<<----gotta fix
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 20),
        axis.text.y = element_text(angle = 360, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + #kills facet headers (= strip)
  theme(panel.background = element_blank()) +
  scale_x_discrete(limits=CHorder, labels=c("PF nest  - Un",
                                            "PF nest  - Ct",
                                            "GT nest  - Un",
                                            "GT nest  - Ct",
                                            "GT adult - Un",
                                            "GT adult - Ct")) +
  theme(axis.text=element_text(size=20), 
        legend.text = element_text(size=20), 
        legend.title = element_text(size=20)) 

FIG1b

# SAVE and skip ps production
save.image(file = "ps.bar.RData")
load(file = "ps.bar.RData")

#################################################################################
# ALPHA DIVERSITY
# https://microbiome.github.io/tutorials/Alphadiversity.html
# get an overview of diversity in samples for all available indices of alpha diversity
alpha_table <- microbiome::alpha(ps.rarefied, index = "all")
kable(head(alpha_table))

# can also search for (1) richness, (2) dominance, (3) rarity, (4) coverage, (5) core, (6) gini (inequality), and (7) evenness
# as an example
rarity_table <- rarity(ps.rarefied, index = "all")
kable(head(rarity_table))

alpha_df <- alpha_table %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

meta_df <- as(sample_data(ps.rarefied), "data.frame") %>%
  rownames_to_column(var = "SampleID")

alpha_with_meta <- alpha_df %>%
  left_join(meta_df, by = "SampleID")


color_map <- c(
  "GTnC" = "#365C8DFF",  
  "GTnH" = "#277F8EFF",
  "PFnC" = "#4AC16DFF",
  "PFnH" = "#9FDA3AFF",
  "GTaC" = "#F57D15FF",
  "GTaH" = "#FAC127FF"
)

obs <- ggplot(alpha_with_meta, aes(
  x = factor(categoryHC, levels = names(color_map)),
  y = observed,
  fill = categoryHC,
  color = categoryHC)) +    # also map color for dots outline
  scale_x_discrete(limits=CHorder, labels=c("F. hypo nest  - Unc",
                                            "F. hypo nest- Con",
                                            "P. majo nest  - Unc",
                                            "P. majo nest  - Con",
                                            "P. majo adult - Unc",
                                            "P. majo adult - Con")) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 1.0, size = 2.2, shape = 21) +  # shape 21 has fill + border
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +  # keep dot border same as fill
  labs(x = "\nsample\n", y = "\nObserved features\n") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),   # removes grey grid
    axis.line = element_line(color = "black", size = 0.8),  # black axis lines with thickness 0.8
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

obs

sha <- ggplot(alpha_with_meta, aes(
  x = factor(categoryHC, levels = names(color_map)),
  y = diversity_shannon,
  fill = categoryHC,
  color = categoryHC)) +    # also map color for dots outline
  scale_x_discrete(limits=CHorder, labels=c("F. hypo nest  - Unc",
                                            "F. hypo nest- Con",
                                            "P. majo nest  - Unc",
                                            "P. majo nest  - Con",
                                            "P. majo adult - Unc",
                                            "P. majo adult - Con")) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 1.0, size = 2.2, shape = 21) +  # shape 21 has fill + border
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +  # keep dot border same as fill
  labs(x = "\nsample\n", y = "\nShannon index\n") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),   # removes grey grid
    axis.line = element_line(color = "black", size = 0.8),  # black axis lines with thickness 0.8
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

sha

# combining 2 plots above
FIG2 <- ggarrange(
  obs, sha, 
  labels = c(" ", " "),
  common.legend = FALSE, 
  legend = FALSE
)

FIG2

# SAVE and skip ps production
#save.image(file = "ps.alpha.RData")
#load(file = "ps.alpha.RData")

#######################################################################################################
## BETA DIVERSITY
## distances = jaccard, bray, weighted unifrac (unifrac, weight = T), unweighted unifrac 
## ord_meths = "DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA"

## basic ordination plots
## https://github.com/joey711/phyloseq/issues/895
## https://github.com/vegandevs/vegan/issues/153

jaccard_dist = phyloseq::distance(ps.rarefied, method = "jaccard", binary = T)
bray_dist = phyloseq::distance(ps.rarefied, method = "bray", weighted = T)

## ordination options - type="samples", color="group_site", shape = "group"
ord.jacc = ordinate(ps.rarefied, 
                    method = "PCoA", 
                    distance = jaccard_dist)

# base ordination
jaccard <- plot_ordination(ps.rarefied, ord.jacc, color = "categoryHC") +
  theme(aspect.ratio = 1)

CHorder<-c("PFnC", "PFnH", "GTnC", "GTnH", "GTaC", "GTaH")

color_map <- c(
  "GTnC" = "#365C8DFF",  
  "GTnH" = "#FAC127FF",
  "PFnC" = "#365C8DFF",
  "PFnH" = "#FAC127FF",
  "GTaC" = "#365C8DFF",
  "GTaH" = "#FAC127FF"
)


jaccard$data$categoryHC <- factor(as.character(jaccard$data$categoryHC),
                                     levels = CHorder)

jaccard_final <- jaccard +
  geom_point(size = 2.5, alpha = 0.5, shape = 21,
             colour = "black", stroke = 0.4,
             aes(fill = categoryHC)) +
  stat_ellipse(
    data = jaccard$data,            
    aes(colour = categoryHC),       
    geom = "path",
    type = "norm",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  scale_x_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),        # tick positions
    labels = scales::number_format(accuracy = 0.1)  # decimal places
  ) +
  scale_y_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),
    labels = scales::number_format(accuracy = 0.1)
  ) +
theme_minimal(base_size = 16) +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(size =20, vjust = -0.5),
    axis.title.y = element_text(size = 20, vjust = 1),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
#  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) + # better options above
  labs(x = "axis 1 (3.5% variation)\n",
       y = "axis 2 (2.2% variation)")

jaccard_final

# FACET
jaccard_faceted <- jaccard_final + 
  facet_wrap(
    vars(zone), 
    ncol = 2, 
    labeller = labeller(zone = c(
      "C" = "uncontaminated",
      "H" = "contaminated")
  )) +
  theme(
    strip.background = element_rect(
      fill = "white",  
      colour = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 16,
      face = "plain",
      colour = "black"
    )
  )

jaccard_faceted

# BRAY
# ordination options - type="samples", color="group_site", shape = "group"
ord.bray = ordinate(ps.rarefied, 
                    method = "PCoA", 
                    distance = bray_dist)


bray <- plot_ordination(ps.rarefied, ord.bray, color = "categoryHC") +
  theme(aspect.ratio = 1)

bray$data$categoryHC <- factor(as.character(bray$data$categoryHC),
                                  levels = CHorder)

bray_final <- bray +
  geom_point(size = 2.5, alpha = 0.5, shape = 21,
             colour = "black", stroke = 0.4,
             aes(fill = categoryHC)) +
  stat_ellipse(
    data = bray$data,            
    aes(colour = categoryHC),       
    geom = "path",
    type = "norm",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  scale_x_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),        # tick positions
    labels = scales::number_format(accuracy = 0.1)  # decimal places
  ) +
  scale_y_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(size = 16, vjust = -0.5),
    axis.title.y = element_text(size = 16, vjust = 1),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) +
  labs(x = "axis 1 (6.2% variation)\n",
       y = "axis 2 (4.9% variation)")

bray_final

bray_faceted <- bray_final + 
  facet_wrap(
    vars(zone), 
    ncol = 2, 
    labeller = labeller(zone = c(
      "C" = "uncontaminated",
      "H" = "contaminated")
    )) +
  theme(
    strip.background = element_rect(
      fill = "white",  
      colour = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 16,
      face = "plain",
      colour = c("black")
    )
  )

bray_faceted

# combine 2 plots.
FIG3 <- ggarrange(
  jaccard_faceted, bray_faceted, 
  labels = c(" ", " "),
  ncol = 1,        
  nrow = 2,       
  common.legend = FALSE, 
  legend = FALSE
)

FIG3

# alternate faceting
jaccard_final.1 <- jaccard +
  geom_point(size = 2.5, alpha = 0.5, shape = 21,
             colour = "black", stroke = 0.4,
             aes(fill = categoryHC)) +
  stat_ellipse(
    data = jaccard$data,            
    aes(colour = categoryHC),       
    geom = "path",
    type = "norm",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  scale_x_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),        # tick positions
    labels = scales::number_format(accuracy = 0.1)  # decimal places
  ) +
  scale_y_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(size = 16, vjust = -0.5),
    axis.title.y = element_text(size = 16, vjust = 1),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) +
  labs(x = "axis 1 (3.5% variation)\n",
       y = "axis 2 (2.2% variation)")

jaccard_final.1

jaccard_faceted.1 <- jaccard_final.1 + ######TÄMÄ
  facet_wrap(
    vars(category), 
    ncol = 3,  
    labeller = labeller(category = c(
      "PFn" = "F. hypoleuca nestlings",
      "GTn" = "P. major nestlings",
      "GTa" = "P. major adults")
    )) +
  theme(
    strip.background = element_rect(
      fill = "grey80",  
      colour = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 16,
      face = "plain",
      colour = c("black")
    )
  )

jaccard_faceted.1

bray_final.1 <- bray +
  geom_point(size = 2.5, alpha = 0.5, shape = 21,
             colour = "black", stroke = 0.4,
             aes(fill = categoryHC)) +
  stat_ellipse(
    data = bray$data,            
    aes(colour = categoryHC),       
    geom = "path",
    type = "norm",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_map) +
  scale_colour_manual(values = color_map) +
  scale_x_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),        # tick positions
    labels = scales::number_format(accuracy = 0.1)  # decimal places
  ) +
  scale_y_continuous(
    limits = c(-0.4, 0.4),
    breaks = seq(-0.4, 0.4, by = 0.4),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(size = 16, vjust = -0.5),
    axis.title.y = element_text(size = 16, vjust = 1),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4)) +
  labs(x = "axis 1 (6.2% variation)\n",
       y = "axis 2 (4.9% variation)")

bray_final.1

bray_faceted.micro <- bray_final.1 + 
  facet_wrap(
    vars(category), 
    ncol = 3,  
    labeller = labeller(category = c(
      "PFn" = "F. hypoleuca nestlings",
      "GTn" = "P. major nestlings",
      "GTa" = "P. major adults")
    )) +
  theme(
    strip.background = element_rect(
      fill = "grey80",  
      colour = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 16,
      face = "plain",
      colour = c("black")
    )
  )

bray_faceted.micro

# combine 2 plots.
FIG3.alt <- ggarrange(
  jaccard_faceted.1, bray_faceted.1, 
  labels = c(" ", " "),
  ncol = 1,        
  nrow = 2,       
  common.legend = FALSE, 
  legend = FALSE
)

FIG3.alt


# SAVE and skip ps production
#save.image(file = "ps.beta.RData")
#load(file = "ps.beta.RData")

ggarrange(
  bray_faceted.diet, bray_faceted.micro, labels = c("A", "B"), ncol = 1, nrow = 2, align = "v",
  common.legend = TRUE, legend = "top"
)

##############################################################################

# subset samples into three main groups
# slice ps to species/lifestage
psPF <- subset_samples(ps.rarefied, species == "PF")
psGT <- subset_samples(ps.rarefied, species == "GT")
psGTa <- subset_samples(psGT, lifestage == "adult")
psGTn <- subset_samples(psGT, lifestage == "nestling")

##############################################################################
# headers for sample location
# make facet label names for header variable
# https://www.datanovia.com/en/blog/how-to-change-ggplot-facet-labels/
#species.labs <- c("great tit", "pied flycatcher")
#names(species.labs) <- c("GT", "PF")
##############################################################################

# PERMANOVA testing (using adonis2) -> marginal estimates
help("distance", package = "phyloseq")

jaccard_dist.psGTa = phyloseq::distance(psGTa, method = "jaccard", binary = T)
bray_dist.psGTa = phyloseq::distance(psGTa, method = "bray", weighted = T)

# test all variates - hash out any that do not apply
adonis2(jaccard_dist.psGTa ~ sample_data(psGTa)$categoryHC,
                     #sample_data(psGTa)$sex +
                     #sample_data(psGTa)$forest_type +
                     #sample_data(psGTa)$body_condition +
                     #sample_data(psGTa)$brood_size,
                     na.action = na.exclude,
                     by = 'margin')

adonis2(bray_dist.psGTa ~ sample_data(psGTa)$categoryHC +
          #sample_data(psGTa)$sex +
          sample_data(psGTa)$forest_type +
          sample_data(psGTa)$body_condition +
          sample_data(psGTa)$brood_size,
          na.action = na.exclude,
          by = 'margin')

# GT nestlings
jaccard_dist.psGTn = phyloseq::distance(psGTn, method = "jaccard", binary = T)
bray_dist.psGTn = phyloseq::distance(psGTn, method = "bray", weighted = T)

adonis2(jaccard_dist.psGTn ~ sample_data(psGTn)$categoryHC +
          sample_data(psGTn)$age +
          sample_data(psGTn)$forest_type +
          #sample_data(psGTn)$body_condition +
          sample_data(psGTn)$hatch_value,
          #sample_data(psGTn)$brood_size,
          na.action = na.exclude,
          by = 'margin')

adonis2(bray_dist.psGTn ~ sample_data(psGTn)$categoryHC +
          #sample_data(psGTn)$age +
          sample_data(psGTn)$forest_type +
          sample_data(psGTn)$body_condition +
          #sample_data(psGTn)$hatch_value +
          sample_data(psGTn)$brood_size,
          na.action = na.exclude,
          by = 'margin')

# PF nestlings
jaccard_dist.psPF = phyloseq::distance(psPF, method = "jaccard", binary = T)
bray_dist.psPF = phyloseq::distance(psPF, method = "bray", weighted = T)

adonis2(jaccard_dist.psPF ~ sample_data(psPF)$categoryHC +
          #sample_data(psPF)$age +
          sample_data(psPF)$forest_type +
          sample_data(psPF)$body_condition,
          #sample_data(psPF)$hatch_value,
          #sample_data(psPF)$brood_size,
        na.action = na.exclude,
        by = 'margin')

adonis2(bray_dist.psPF ~ sample_data(psPF)$categoryHC +
          #sample_data(psPF)$age +
          sample_data(psPF)$forest_type +
          sample_data(psPF)$body_condition +
          #sample_data(psPF)$hatch_value +
          sample_data(psPF)$brood_size,
        na.action = na.exclude,
        by = 'margin')


# SAVE and skip ps production
#save.image(file = "ps.permanova.RData")
#load(file = "ps.RData")

######################################################
# Linear models of alpha diversity & ecology
# https://ourcodingclub.github.io/tutorials/mixed-models/
# lineaarinen malli

library(lme4)
library(lmerTest)
library(modelr)

#save residuals of linear models (used to determine body condition)
rare_res<-bclm$residuals
plot(fitted(bclm), rare_res)
abline(0,0)
pfcBC<-add_residuals(pfc, bclm, var="body_condition")

#cheching linear models with dharma
library(DHARMa, lib.loc = "/appl/soft/math/r-env/430/430-rpackages")
library(DHARMa)

resp = simulateResiduals(Malli1)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = Malli1, n = 1000)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

###ALPHA DIVERSITY SHANNON&observed
#Great tit nestling shannon
Malli1= lmer(Shannon ~ age+categoryHC+age:categoryHC+brood_size+hatch_value+forest_type+body_condition+(1|box)+(1|plot), data = GTNrare)
#Malli1= lmer(Shannon ~ age+radiation+(1|box)+(1|plot), data = GTNrare)

resp = simulateResiduals(Malli1)
plot(resp, rank = T)
summary(Malli1, ddf="Kenward-Roger")
anova(Malli1, ddf="Kenward-Roger") 

#observed
Malli1= lmer(Observedlog ~ age+categoryHC+age:categoryHC+hatch_value+forest_type+body_condition+(1|box)+(1|plot), data = GTNrare)
#Malli1= lmer(Observedlog ~ age+categoryHC+brood_size+(1|box)+(1|plot), data = GTNrare)

summary(Malli1, ddf="Kenward-Roger")
anova(Malli1, ddf="Kenward-Roger") 
resp = simulateResiduals(Malli1)
plot(resp, rank = T)

##########################################################
#Pied flycatcher
Malli1= lmer(Shannonlog ~ age+categoryHC+age:categoryHC+brood_size+hatch_value+forest_type+body_condition+(1|box)+(1|plot), data = PFrare)
#Malli1= lmer(Shannonlog ~ categoryHC+hatch_value+(1|box)+(1|plot), data = PFrare)

summary(Malli1, ddf="Kenward-Roger")
anova(Malli1, ddf="Kenward-Roger")

#obs
Malli1= lmer(Observedlog ~ age+categoryHC+brood_size+hatch_value+forest_type+body_condition+(1|box)+(1|plot), data = PFrare)
#Malli1= lmer(Observedlog ~ radiation+(1|box)+(1|plot), data = PFrare)

summary(Malli1, ddf="Kenward-Roger")
anova(Malli1, ddf="Kenward-Roger") 
resp = simulateResiduals(Malli1)
plot(resp, rank = T)

#########################################################
#Great Tit adults shannon
Malli1= lmer(Shannonlog ~ categoryHC+sex+forest_type+brood_size+(1|box)+(1|plot), data = GTArare)
#Malli1= lmer(Shannonlog ~ categoryHC+brood_size+(1|box)+(1|plot), data = GTArare)

#observed
Malli1= lmer(Shannonlog ~ categoryHC+sex+forest_type+brood_size+(1|box), data = GTArare)
#Malli1= lmer(Shannonlog ~ radiation+brood_size+(1|box), data = GTArare)

resp = simulateResiduals(Malli1)
plot(resp, rank = T)
summary(Malli1, ddf="Kenward-Roger")
anova(Malli1, ddf="Kenward-Roger") 

library(glmmTMB)



###############################################################################
