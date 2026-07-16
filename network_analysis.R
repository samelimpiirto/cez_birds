# clear environment
#rm(list=ls())

# own packages - in project libpath
.libPaths(c("/projappl/pwatts/rpackages_432", .libPaths()))
libpath <- .libPaths()[1]
# this command can be used to check that the folder is now visible:
.libPaths() 

## --- FORCE REINSTALL graphlayouts (required for ggraph) ---
#install.packages("graphlayouts", lib = libpath, dependencies = TRUE)

## --- RESTART SAFELY IN SESSION ---
#detach("package:graphlayouts", unload = TRUE, character.only = TRUE)

## LOAD IN CORRECT ORDER ---
library(igraph)
library(graphlayouts)
library(ggraph)

#packageVersion("graphlayouts")

# set working directory
setwd("/scratch/pwatts/bird")

# setup environment
library(phyloseq)
library(ggplot2)
library(dplyr)
library(BiocManager)
#BiocManager::install("UCSC.utils", lib = libpath, ask = FALSE, update = FALSE)
#library(UCSC.utils)
library(ggClusterNet)
#library(igraph)
library(tibble)
library(tidyverse) #if needed?
library(network)
#library(ggraph)
library(SpiecEasi)

#load(file = "cherno.bird.diet.2026.RData")
#h1 = 0.1
#ps.bird.gtn.01 = tip_glom(ps.bird.gtn, h = h1)

#h2 = 0.05
#ps.bird.gtn.005 = tip_glom(ps.bird.gtn, h = h2)

#h3 = 0.2
#ps.bird.gtn.02 = tip_glom(ps.bird.gtn, h = h3)

#ps.diet.gtn.01 = tip_glom(ps.diet.gtn, h = h1)
#ps.diet.gtn.005 = tip_glom(ps.diet.gtn, h = h2)
#ps.diet.gtn.02 = tip_glom(ps.diet.gtn, h = h3)

load(file = "cherno.bird.diet.agglomerate.RData")
## PS objects
#ps.diet.gtn
#ps.bird.gtn

## BACTERIA
bact_net_obj <- network.pip(
  ps = ps.bird.gtn.02,
  N = 0,
  ra = 0.05,
  big = TRUE,
  select_layout = TRUE,
  layout_net = "model_igraph",
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "pearson",
  label = FALSE,
  group = "zone",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 500,
  R = 50,
  ncpus = 4
)

#save.image(file = "cherno.bird.bact.net.2.RData")

## DIET
diet_net_obj <- network.pip(
  ps = ps.diet.gtn,
  N = 0,
  ra = 0.05,
  big = TRUE,
  select_layout = TRUE,
  layout_net = "model_igraph",
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "pearson",
  label = FALSE,
  group = "zone",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 500,
  R = 50,
  ncpus = 4
)

#save.image(file = "cherno.bird.diet.net.2.RData")
#load(file = "cherno.bird.bact.net.2.RData")


###############################################################################
## 1. BACTERIA NETWORK
###############################################################################

### EXTRACT NETWORK VISUALISATION AND CORRELATION MATRICES ---
plot <- bact_net_obj[[1]]
p0 <- plot[[1]]      
p0.1 <- plot[[2]]    
p0.2 <- plot[[3]]    

dat <- bact_net_obj[[2]]
cortab <- dat$net.cor.matrix$cortab

#install.packages("graphlayouts", lib = .libPaths()[1])
library(graphlayouts)
packageVersion("graphlayouts")

### NETWORK COMPARISON STATISTICS
#dat.compare <- module.compare.net.pip(
#  ps = NULL,
#  corg = cortab,
#  degree = TRUE,
#  zipi = TRUE,
#  r.threshold = 0.6,
#  p.threshold = 0.05,
#  method = "pearson",
#  padj = TRUE,
#  n = 4
#)
#res <- dat.compare[[1]]
#head(res)

### CUSTOMISED NETWORK VISUALISATION ---
node <- dat$net.cor.matrix$node
edge <- dat$net.cor.matrix$edge

p <- ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor),
               data = edge, size = 0.05, alpha = 0.8) +
  geom_point(aes(x = X1, y = X2, fill = Phylum, size = igraph.degree),
             pch = 21, data = node, color = "grey80") +
  facet_wrap(. ~ label, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("black", "blue")) +
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.grid = element_blank())

### CUSTOMISE NETWORKS ---
#library(tibble)
otu_mat <- phyloseq::otu_table(ps.bird.gtn.02)
if (phyloseq::taxa_are_rows(ps.bird.gtn.02)) {
  otu_mat <- t(otu_mat)
}
otu_mean <- colMeans(otu_mat)

otu_mean_df <- data.frame(elements = names(otu_mean),
                          mean = as.numeric(otu_mean),
                          stringsAsFactors = FALSE)

node <- dplyr::left_join(node, otu_mean_df, by = "elements")

top_phyla <- node %>%
  dplyr::filter(!is.na(Phylum)) %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(total = sum(mean, na.rm = TRUE)) %>%
  dplyr::arrange(desc(total)) %>%
  dplyr::slice_head(n = 3)

top_phyla_vec <- top_phyla$Phylum

custom_colors <- c(
  " Firmicutes" = "#440154FF",
  " Actinobacteriota" = "#238A8DFF",
  " Proteobacteria" = "#FDE725FF"
)

custom_colors <- custom_colors[names(custom_colors) %in% top_phyla_vec]
custom_colors <- c(custom_colors, Other = "grey80")

node$Phylum_colored <- ifelse(
  node$Phylum %in% names(custom_colors),
  node$Phylum,
  "Other"
)

colnames(node)

desired_order <- unique(node$Group)
node$label <- factor(node$Group, levels = desired_order)
edge$label <- factor(edge$Group, levels = desired_order)

p_custom <- ggplot() +
  geom_segment(aes(x = X1, y = Y1,
                   xend = X2, yend = Y2,
                   color = cor),
               data = edge, size = 0.5, alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,
                 fill = Phylum_colored,
                 size = igraph.degree),
             pch = 21, data = node, color = "grey90") +
  facet_wrap(. ~ label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = custom_colors) +
  scale_colour_manual(values = c("grey5", "black")) +
  scale_size(range = c(1, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray90", color = "gray90")
  )

p_custom <- ggplot() +
  geom_segment(aes(x = X1, y = Y1,
                   xend = X2, yend = Y2,
                   color = cor),
               data = edge, size = 0.4, alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,
                 fill = Phylum_colored,
                 size = igraph.degree),
             pch = 21, data = node, color = "grey30") +
  facet_wrap(. ~ label, scales = "fixed", nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_manual(values = custom_colors) +
  scale_colour_manual(values = c("grey5", "black")) +
  scale_size(range = c(2, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray90", color = "gray90")
  )

p_custom

p_custom <- p_custom +
  facet_wrap(
    ~ label,
    scales = "fixed",
    nrow = 1,
    labeller = labeller(
      label = c(
        "C" = "uncontaminated",
        "H" = "contaminated"
      )
    )
  )

### NETWORK PROPERTIES ---
node_df <- dat$net.cor.matrix$node
edge_df <- dat$net.cor.matrix$edge
groups <- unique(node_df$Group)

results <- list()
net_props_df <- NULL
graphs <- list()

for (group_id in groups) {
  
  nodes_group <- dplyr::filter(node_df, Group == group_id)
  edges_group <- dplyr::filter(edge_df, Group == group_id)
  
  g <- igraph::graph_from_data_frame(
    d = edges_group[, c("OTU_1", "OTU_2", "weight")],
    vertices = nodes_group[, c("elements", "Phylum")],
    directed = FALSE
  )
  
  props <- net_properties.4(g, n.hub = TRUE)
  props_df <- as.data.frame(t(props))
  props_df$Group <- group_id
  
  net_props_df <- dplyr::bind_rows(net_props_df, props_df)
  graphs[[group_id]] <- g
}

print(net_props_df)
print(graphs)

plot(graphs[[1]])

#save.image(file = "cherno.bird.diet.net.3.RData")
#load(file = "cherno.bird.diet.net.3.RData")

### SAMPLE-LEVEL NETWORK PROPERTIES ---
dat.f2 <- NULL
for (group_name in names(cortab)) {
  pst <- ps.bird.gtn.02 %>%
    subset_samples.wt("zone", group_name) %>%
    subset_taxa.wt("OTU", colnames(cortab[[group_name]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")
  
  dat.f <- netproperties.sample(pst = pst, cor = cortab[[group_name]])
  dat.f$Group <- group_name
  dat.f2 <- dplyr::bind_rows(dat.f2, dat.f)
}

dat.f2$Group <- as.factor(dat.f2$Group)

r <- ggplot(dat.f2, aes(x = Group, y = no.clusters, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank())

#####

bact_node <- bact_net_obj[[2]]$net.cor.matrix$node
bact_edge <- bact_net_obj[[2]]$net.cor.matrix$edge
bact_cortab <- bact_net_obj[[2]]$net.cor.matrix$cortab

bact_plot <- ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor),
               data = bact_edge, size = 0.3, alpha = 0.6) +
  geom_point(aes(x = X1, y = X2, fill = Phylum, size = igraph.degree),
             data = bact_node, shape = 21, color = "grey30") +
  facet_wrap(. ~ Group) +
  scale_colour_manual(values = c("black", "grey60")) +
  scale_size(range = c(1, 4)) +
  theme_void()

bact_compare <- module.compare.net.pip(
  ps = NULL,
  corg = bact_cortab,
  degree = TRUE,
  zipi = TRUE,
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "pearson",
  padj = FALSE,
  n = 4
)

bact_stats <- bact_compare[[1]]
bact_groups <- unique(bact_node$Group)
bact_props_list <- list()

for (g in bact_groups) {
  nodes_g <- bact_node %>% filter(Group == g)
  edges_g <- bact_edge %>% filter(Group == g)
  
  g_obj <- graph_from_data_frame(
    d = edges_g[, c("OTU_1", "OTU_2", "weight")],
    vertices = nodes_g[, c("elements", "Phylum")],
    directed = FALSE
  )
  
  props <- net_properties.4(g_obj, n.hub = TRUE)
  props_df <- as.data.frame(t(props))
  props_df$Group <- g
  
  bact_props_list[[g]] <- props_df
}

bact_props <- bind_rows(bact_props_list)

bact_plot
bact_stats
bact_props



###############################################################################
## 2. DIET NETWORK
###############################################################################

### EXTRACT NETWORK VISUALISATION AND CORRELATION MATRICES (DIET) ---
diet_plot <- diet_net_obj[[1]]
diet_p0 <- diet_plot[[1]]      
diet_p0.1 <- diet_plot[[2]]    
diet_p0.2 <- diet_plot[[3]]    

diet_dat <- diet_net_obj[[2]]
diet_cortab <- diet_dat$net.cor.matrix$cortab

### CUSTOMISED NETWORK VISUALISATION ---
diet_node <- diet_dat$net.cor.matrix$node
diet_edge <- diet_dat$net.cor.matrix$edge

diet_p <- ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor),
               data = diet_edge, size = 0.05, alpha = 0.8) +
  geom_point(aes(x = X1, y = X2, fill = Phylum, size = igraph.degree),
             pch = 21, data = diet_node, color = "grey80") +
  facet_wrap(. ~ Group, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("black", "blue")) +
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.grid = element_blank())

### CUSTOMISE NETWORKS ---
otu_mat_diet <- phyloseq::otu_table(ps.diet.gtn)
if (phyloseq::taxa_are_rows(ps.diet.gtn)) {
  otu_mat_diet <- t(otu_mat_diet)
}
otu_mean_diet <- colMeans(otu_mat_diet)

otu_mean_df_diet <- data.frame(
  elements = names(otu_mean_diet),
  mean = as.numeric(otu_mean_diet),
  stringsAsFactors = FALSE
)

diet_node <- dplyr::left_join(diet_node, otu_mean_df_diet, by = "elements")

diet_node$Order <- trimws(diet_node$Order)

top_phyla_diet <- diet_node %>%
  dplyr::filter(!is.na(Order)) %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(total = sum(mean, na.rm = TRUE)) %>%
  dplyr::arrange(desc(total)) %>%
  dplyr::slice_head(n = 3)

top_phyla_vec_diet <- top_phyla_diet$Order



custom_colors_diet <- c(
  "Lepidoptera" = "#440154FF",
  "Diptera" = "#238A8DFF",
  "unassigned" = "#FDE725FF"
)

custom_colors_diet <- custom_colors_diet[names(custom_colors_diet) %in% top_phyla_vec_diet]
custom_colors_diet <- c(custom_colors_diet, Other = "grey80")

diet_node$Phylum_colored <- ifelse(
  diet_node$Order %in% names(custom_colors_diet),
  diet_node$Order,
  "Other"
)

diet_node$label <- factor(diet_node$Group, levels = unique(diet_node$Group))
diet_edge$label <- factor(diet_edge$Group, levels = unique(diet_edge$Group))

diet_p_custom <- ggplot() +
  geom_segment(aes(x = X1, y = Y1,
                   xend = X2, yend = Y2,
                   color = cor),
               data = diet_edge, size = 0.4, alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,
                 fill = Phylum_colored,
                 size = igraph.degree),
             pch = 21, data = diet_node, color = "grey30") +
  facet_wrap(. ~ label, scales = "fixed", nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_manual(values = custom_colors_diet) +
  scale_colour_manual(values = c("grey5", "black")) +
  scale_size(range = c(2, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray90", color = "gray90")
  )

diet_p_custom

### NETWORK PROPERTIES ---
diet_node_df <- diet_dat$net.cor.matrix$node
diet_edge_df <- diet_dat$net.cor.matrix$edge
diet_groups <- unique(diet_node_df$Group)

diet_net_props_df <- NULL
diet_graphs <- list()

for (group_id in diet_groups) {
  
  nodes_group <- dplyr::filter(diet_node_df, Group == group_id)
  edges_group <- dplyr::filter(diet_edge_df, Group == group_id)
  
  g_diet <- igraph::graph_from_data_frame(
    d = edges_group[, c("OTU_1", "OTU_2", "weight")],
    vertices = nodes_group[, c("elements", "Phylum")],
    directed = FALSE
  )
  
  props <- net_properties.4(g_diet, n.hub = TRUE)
  props_df <- as.data.frame(t(props))
  props_df$Group <- group_id
  
  diet_net_props_df <- dplyr::bind_rows(diet_net_props_df, props_df)
  diet_graphs[[group_id]] <- g_diet
}

print(diet_net_props_df)
print(diet_graphs)

plot(diet_graphs[[1]])

### SAMPLE-LEVEL NETWORK PROPERTIES ---
diet_dat_f2 <- NULL

for (group_name in names(diet_cortab)) {
  
  pst_diet <- ps.diet.gtn %>%
    subset_samples.wt("zone", group_name) %>%
    subset_taxa.wt("OTU", colnames(diet_cortab[[group_name]])) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    scale_micro("rela")
  
  dat_f <- netproperties.sample(pst = pst_diet, cor = diet_cortab[[group_name]])
  dat_f$Group <- group_name
  
  diet_dat_f2 <- dplyr::bind_rows(diet_dat_f2, dat_f)
}

diet_dat_f2$Group <- as.factor(diet_dat_f2$Group)

diet_r <- ggplot(diet_dat_f2, aes(x = Group, y = no.clusters, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank())

diet_r

#save.image(file = "cherno.bird.diet.net.4.RData")
#load(file = "cherno.bird.diet.net.4.RData")

######


diet_node <- diet_net_obj[[2]]$net.cor.matrix$node
diet_edge <- diet_net_obj[[2]]$net.cor.matrix$edge
diet_cortab <- diet_net_obj[[2]]$net.cor.matrix$cortab

diet_plot <- ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor),
               data = diet_edge, size = 0.3, alpha = 0.6) +
  geom_point(aes(x = X1, y = X2, fill = Phylum, size = igraph.degree),
             data = diet_node, shape = 21, color = "grey30") +
  facet_wrap(. ~ Group) +
  scale_colour_manual(values = c("black", "grey60")) +
  scale_size(range = c(1, 4)) +
  theme_void()

diet_compare <- module.compare.net.pip(
  ps = NULL,
  corg = diet_cortab,
  degree = TRUE,
  zipi = TRUE,
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "pearson",
  padj = FALSE,
  n = 4
)
diet_stats <- diet_compare[[1]]

diet_groups <- unique(diet_node$Group)
diet_props_list <- list()

for (g in diet_groups) {
  nodes_g <- diet_node %>% filter(Group == g)
  edges_g <- diet_edge %>% filter(Group == g)
  
  g_obj <- graph_from_data_frame(
    d = edges_g[, c("OTU_1", "OTU_2", "weight")],
    vertices = nodes_g[, c("elements", "Phylum")],
    directed = FALSE
  )
  
  props <- net_properties.4(g_obj, n.hub = TRUE)
  props_df <- as.data.frame(t(props))
  props_df$Group <- g
  
  diet_props_list[[g]] <- props_df
}
diet_props <- bind_rows(diet_props_list)

diet_plot
diet_stats
diet_props

##format TABLES
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

bind_rows(
  
  # --- MICROBIAL NETWORK ---
  net_props_df %>%
    mutate(Network = "Microbial") %>%
    rename(
      edges = `num.edges(L)`,
      pos_edges = num.pos.edges,
      neg_edges = num.neg.edges,
      vertices = `num.vertices(n)`,
      connectance = `Connectance(edge_density)`,
      avg_degree = `average.degree(Average K)`,
      avg_path = average.path.length,
      diameter = diameter,
      edge_conn = edge.connectivity,
      clustering = `mean.clustering.coefficient(Average.CC)`,
      clusters = no.clusters,
      centr_degree = centralization.degree,
      centr_between = centralization.betweenness,
      centr_close = centralization.closeness,
      rel_modularity = `RM(relative.modularity)`,
      modularity = modularity.net,
      modularity_rand = modularity_random,
      keystone = the.number.of.keystone.nodes
    ),
  
  # --- DIET NETWORK ---
  diet_net_props_df %>%
    mutate(Network = "Diet") %>%
    rename(
      edges = `num.edges(L)`,
      pos_edges = num.pos.edges,
      neg_edges = num.neg.edges,
      vertices = `num.vertices(n)`,
      connectance = `Connectance(edge_density)`,
      avg_degree = `average.degree(Average K)`,
      avg_path = average.path.length,
      diameter = diameter,
      edge_conn = edge.connectivity,
      clustering = `mean.clustering.coefficient(Average.CC)`,
      clusters = no.clusters,
      centr_degree = centralization.degree,
      centr_between = centralization.betweenness,
      centr_close = centralization.closeness,
      rel_modularity = `RM(relative.modularity)`,
      modularity = modularity.net,
      modularity_rand = modularity_random,
      keystone = the.number.of.keystone.nodes
    )
  
) %>%
  
  # --- CLEAN + FORMAT ---
  mutate(
    centr_close = ifelse(is.nan(centr_close), NA, centr_close),
    
    across(
      c(connectance, avg_degree, avg_path, diameter,
        clustering, centr_degree, centr_between,
        rel_modularity, modularity, modularity_rand),
      ~ round(., 3)
    ),
    
    across(
      c(edges, pos_edges, neg_edges, vertices,
        edge_conn, clusters, keystone),
      ~ as.integer(.)
    )
  ) %>%
  
  # --- LONG FORMAT ---
  pivot_longer(
    cols = -c(Group, Network),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  
  # --- WIDE: GROUPS AS COLUMNS ---
  pivot_wider(
    names_from = Group,
    values_from = Value
  ) %>%
  
  # --- ORDER ---
  mutate(
    Metric = factor(Metric, levels = c(
      "edges", "pos_edges", "neg_edges", "vertices",
      "connectance", "avg_degree", "avg_path",
      "diameter", "edge_conn", "clustering",
      "clusters", "centr_degree", "centr_between",
      "rel_modularity", "modularity", "modularity_rand",
      "keystone"
    ))
  ) %>%
  arrange(Network, Metric) %>%
  
  # --- RENDER ---
  kable(
    format = "html",
    caption = "Network properties for microbial and diet networks (C vs H)"
  ) %>%
  kable_styling(full_width = FALSE) %>%
  collapse_rows(columns = 1, valign = "top")

