# Comparing the diversity of Bicuar with that of the other SEOSAW plots
# John Godlee (johngodlee@gmail.com)
# 2019_11_28

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(betapart)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(stargazer)

source("scripts/functions.R")

# Import data ----
trees_list <- readRDS("data/trees_list.rds")
plot_clim_list <- readRDS("data/plot_clim_list.rds")
trees_ab_mat_list <- readRDS("data/trees_ab_mat.rds")

# Compare diversity found among plot groups ----

# Create clean tree abundance matrix of all big plots
trees_ab_mat_big <- bind_rows(trees_ab_mat_list[[1]], 
  trees_ab_mat_list[[3]],
  trees_ab_mat_list[[4]],
  trees_ab_mat_list[[5]])
trees_ab_mat_big[is.na(trees_ab_mat_big)] <- 0

trees_ab_mat_big_clean <- dplyr::select(trees_ab_mat_big, -plotcode)
rownames(trees_ab_mat_big_clean) <- trees_ab_mat_big$plotcode

plot_vec <- rownames(trees_ab_mat_big_clean)
plot_group <- data.frame(plotcode = as.character(plot_vec), 
  group = case_when(
    grepl("ABG", plot_vec) ~ "bicuar",
    grepl("TKW", plot_vec) ~ "kilwa",
    grepl("DKS", plot_vec) ~ "drc",
    grepl("MGR", plot_vec) ~ "nham",
    TRUE ~ "seosaw"))

# Order plot_clim by tree_ab_mat plotcodes
plot_clim_big <- plot_clim_list[c(1, 3:5)]
plot_clim_big_df <- do.call(rbind, plot_clim_big)

plot_clim_s_all <- left_join(plot_group, plot_clim_big_df, by = c("plotcode" = "plotcode"))

# Calculate Sorensen and Bray beta diversity indices for each plot
##' Creates `dist` object
dist_sor_abund <- beta.pair.abund(trees_ab_mat_big_clean, index.family="bray")

# Calc. multivar. homogeneity of group dispersions inside and outside exclosures
bd_sor_abund <- betadisper(dist_sor_abund[[3]], plot_group$group)

# Create plots of community difference
bd_sorr_abund_plot <- plot(bd_sor_abund)

# Run an ANOVA
##' Is the species community composition different between 
##' groups of plots?
##' If p<0.05, significant different between groups of plots
anova(bd_sor_abund)
tukey_betadisper <- TukeyHSD(bd_sor_abund)

# Create table of environmental fits to the ordination
fileConn <- file(paste0("include/tukey_betadisper.tex"))
writeLines(
  stargazer(tukey_betadisper$group, summary = FALSE, label = "tukey_betadisper", 
    digit.separate = 0, rownames = FALSE),
  fileConn)
close(fileConn)

# Recreate the abundance plot above in ggplot2
# Extract data from plot object
bd_sorr_abund_plot[[1]] <- as.data.frame(bd_sorr_abund_plot[[1]])
bd_sorr_abund_plot[[1]]$group <- plot_group$group
bd_sorr_abund_plot[[1]]$plot <- plot_group$plotcode
bd_sorr_abund_plot[[2]] <- as.data.frame(bd_sorr_abund_plot[[2]])
bd_sorr_abund_plot[[2]]$group <- row.names(bd_sorr_abund_plot[[2]]) 
bd_sorr_abund_plot[[3]] <- left_join(bd_sorr_abund_plot[[1]], bd_sorr_abund_plot[[2]], by = "group")
bd_sorr_abund_plot[[3]]$line <- row.names(bd_sorr_abund_plot[[3]])
bd_a <- bd_sorr_abund_plot[[3]][,c(1:3,7)]
bd_b <- bd_sorr_abund_plot[[3]][,c(5,6,3,7)]
names(bd_b) <- names(bd_a)
lines_df <- rbind(bd_a, bd_b)
names(lines_df) <- c("PCoA1", "PCoA2", "group", "line")

hulls_list <- split(bd_sorr_abund_plot[[1]], bd_sorr_abund_plot[[1]]$group)
hulls_list <- lapply(hulls_list, find_hull)
hulls <- do.call(rbind, hulls_list)

# Create plot
s_all_beta_div_plot <- ggplot() + 
  geom_line(data = lines_df, 
    aes(x = PCoA1, y = PCoA2, colour = group, group = line), 
    show.legend = FALSE) + 
  geom_polygon(data = hulls,
    aes(x = PCoA1, y = PCoA2, fill = group), 
    alpha = 0.5) + 
  geom_point(data = bd_sorr_abund_plot[[1]], 
    aes(x = PCoA1, y = PCoA2, fill = group),
    shape = 21, size = 2) + 
  geom_point(data = bd_sorr_abund_plot[[2]], 
    aes(x = PCoA1, y = PCoA2, fill = group),
    shape = 23, size = 4) + 
  # geom_label(data = bd_sorr_abund_plot[[1]],
  #   aes(x = PCoA1, y = PCoA2, fill = group, label = plot),
  #   size = 4, show.legend = FALSE) +
  theme_classic() + 
  scale_fill_manual(name = "", values = big_pal, labels = c("Angola", "DRC", "Tanzania", "Mozambique"))

pdf(file = paste0("img/s_all_beta_div.pdf"), width = 8, height = 6)
s_all_beta_div_plot
dev.off()

# NMDS between groups of plots ----
# Test optimal dimensions
NMDS.scree(trees_ab_mat_big_clean)

##' Use Bray Curtis distance as it's not biased by null values in the matrix
tree_ab_nmds <- metaMDS(trees_ab_mat_big_clean, distance = "bray", try = 500, trymax = 500, 
  k = 4, autotransform = FALSE)

# Stress plot
stressplot(tree_ab_nmds)

# How many tries (Restarts)?
tree_ab_nmds$tries

# Final stress value?
nmdsstress <- sprintf("%.2f", round(tree_ab_nmds$stress, 2))

# Number of axes?
tree_ab_nmds$ndim

# Extract site (plot) scores from NMDS analysis
plot_scores <- as.data.frame(scores(tree_ab_nmds))  

# Extract species scores from NMDS analysis
species_scores <- as.data.frame(scores(tree_ab_nmds, "species")) 
species_scores$species_binomial <- rownames(species_scores)

# Environmental fit with NMDS
tree_ab_envfit <- envfit(tree_ab_nmds, plot_clim_s_all[,c("mat", "map", "mat_sd", "cwd")], permutations = 999)

tree_ab_envfit
plot(tree_ab_nmds, type = "p")
plot(tree_ab_envfit, p.max = 0.05)

# Get arrow vectors
tree_ab_envfit_arrows <- data.frame(tree_ab_envfit$vectors$arrows)
tree_ab_envfit_arrows$var <- rownames(tree_ab_envfit_arrows)
tree_ab_envfit_arrows$r2 <- tree_ab_envfit$vectors$r
tree_ab_envfit_arrows$p <- tree_ab_envfit$vectors$pvals
tree_ab_envfit_arrows$var_plot <- case_when(
  tree_ab_envfit_arrows$var == "map" ~ "MAP", 
  tree_ab_envfit_arrows$var == "mat" ~ "MAT", 
  tree_ab_envfit_arrows$var == "cwd" ~ "CWD", 
  tree_ab_envfit_arrows$var == "mat_sd" ~ "MAT SD"
)

nmds_envfit_plot <- ggplot() + 
  geom_hline(aes(yintercept = 0), linetype = 2) + 
  geom_vline(aes(xintercept = 0), linetype = 2) + 
  # geom_point(data = species_scores,
  #   aes(x = NMDS1, y = NMDS2)) + 
  # geom_label_repel(data = plot_scores,
  #   aes(x = NMDS1, y = NMDS2, label =  plot_group$plotcode, colour = plot_group$group),
  #   label.padding = 0.15, box.padding = 0.5) +
  # geom_label_repel(data = filter(species_scores, NMDS1 < -0.5 & NMDS2 <0),
  #   aes(x = NMDS1, y = NMDS2, label =  species_binomial),
  #   label.padding = 0.15, box.padding = 0.5) +
  geom_point(data = plot_scores,
    aes(x = NMDS1, y = NMDS2, fill = plot_group$group), shape = 23, size = 3) +
  coord_equal() +
  theme_classic() + 
  geom_segment(data = tree_ab_envfit_arrows, 
    aes(xend = NMDS1*r2*2, yend = NMDS2*r2*2, x = 0, y = 0), 
    arrow = arrow(length = unit(0.05, "npc")),
    colour = "blue") + 
  geom_text(data = tree_ab_envfit_arrows,
    aes(x = NMDS1, y = NMDS2, label = var_plot),
    colour = "blue", 
    nudge_y = c(-0.35,0.1,0,0),
    nudge_x = c(0.75,0,0.5,0.25)) + 
  scale_fill_manual(name = "", values = big_pal, labels = c("Angola", "DRC", "Tanzania", "Mozambique"))

nmds_plot <- ggplot() + 
  geom_hline(aes(yintercept = 0), linetype = 2) + 
  geom_vline(aes(xintercept = 0), linetype = 2) + 
  geom_point(data = species_scores,
    aes(x = NMDS1, y = NMDS2)) + 
  geom_label_repel(data = plot_scores,
    aes(x = NMDS1, y = NMDS2, label =  plot_group$plotcode, colour = plot_group$group),
    label.padding = 0.08, box.padding = 0.2, show.legend = FALSE) +
  geom_point(data = plot_scores,
    aes(x = NMDS1, y = NMDS2, fill = plot_group$group), shape = 23, size = 3) +
  coord_equal() +
  theme_classic() + 
  scale_fill_manual(name = "", values = big_pal, labels = c("Angola", "DRC", "Tanzania", "Mozambique"))

nmds_plot_species <- ggplot() + 
  geom_hline(aes(yintercept = 0), linetype = 2) + 
  geom_vline(aes(xintercept = 0), linetype = 2) + 
  geom_point(data = species_scores,
    aes(x = NMDS1, y = NMDS2)) + 
  geom_label_repel(data = filter(species_scores,
      NMDS1 < 0 & NMDS2 > 0),
    aes(x = NMDS1, y = NMDS2, label = species_binomial),
    label.padding = 0.08, box.padding = 0.2, show.legend = FALSE, min.segment.length = 0) + 
  geom_point(data = plot_scores,
    aes(x = NMDS1, y = NMDS2, fill = plot_group$group), shape = 23, size = 3) +
  coord_equal() +
  theme_classic() + 
  scale_fill_manual(name = "", values = big_pal, labels = c("Angola", "DRC", "Tanzania", "Mozambique"))

##' All plot groups are distinct from each other by species composition
##' Bicuar plots 9, 13, and 15 are very very different from the others 

pdf(file = paste0("img/all_nmds_envfit.pdf"), width = 8, height = 6)
nmds_envfit_plot
dev.off()

pdf(file = paste0("img/all_nmds.pdf"), width = 10, height = 8)
nmds_plot
dev.off()

pdf(file = paste0("img/all_nmds_species.pdf"), width = 10, height = 8)
nmds_plot_species
dev.off()

tree_ab_envfit_arrows$p_format <- as.character(p_format(tree_ab_envfit_arrows$p))
tree_ab_envfit_arrows$NMDS1 <- as.character(round(tree_ab_envfit_arrows$NMDS1, digits = 3))
tree_ab_envfit_arrows$NMDS2 <- as.character(round(tree_ab_envfit_arrows$NMDS2, digits = 3))
tree_ab_envfit_arrows$r2 <- as.character(round(tree_ab_envfit_arrows$r2, digits = 2))
tree_ab_envfit_arrows <- tree_ab_envfit_arrows[,c("var_plot", "NMDS1", "NMDS2", "r2", "p_format")]

# Create table of environmental fits to the ordination
fileConn <- file(paste0("include/nmds_envfit.tex"))
writeLines(
  stargazer(tree_ab_envfit_arrows, summary = FALSE, label = "nmds_envfit", 
    digit.separate = 0, rownames = FALSE),
  fileConn)
close(fileConn)

# Output numerical figures ----
nmds_env_extract <- function(x){
  paste0("R\\textsuperscript{2} = ", 
    tree_ab_envfit_arrows[row.names(tree_ab_envfit_arrows) == x,]$r2,  
    ", ",
    tree_ab_envfit_arrows[row.names(tree_ab_envfit_arrows) == x,]$p_format)
}
nmds_mat <- nmds_env_extract("mat")
nmds_map <- nmds_env_extract("map")
nmds_mat_sd <- nmds_env_extract("mat_sd")
nmds_map_sd <- nmds_env_extract("cwd")

fileConn <- file("include/beta_div_figures.tex")
writeLines(
  c(
    paste0("\\newcommand{\\nmdsstress}{", nmdsstress, "}"),
    paste0("\\newcommand{\\nmdsmat}{", nmds_mat, "}"),
    paste0("\\newcommand{\\nmdsmap}{", nmds_map, "}"),
    paste0("\\newcommand{\\nmdsmatsd}{", nmds_mat_sd, "}"),
    paste0("\\newcommand{\\nmdsmapsd}{", nmds_map_sd, "}")),
  fileConn)
close(fileConn)
