# Summarising each stem dataset by base stem (tree)
# John Godlee (johngodlee@gmail.com)
# 2019_11_29

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(dplyr)
library(tidyr)

source("scripts/functions.R")

# Import data ----
stems_list <- readRDS("data/stems_list.rds")

trees_list <- list()

# Collapse stem measurements to trees
trees_list[[1]] <- stems_list$bicuar %>% 
  group_by(base_stem_id) %>%
  summarise(
    plotcode = first(na.omit(plotcode)),
    stem_id = paste0(stem_id, collapse = ","),
    subplot_20_100 = first(na.omit(subplot_20_100)),
    subplot_20_50 = first(na.omit(subplot_20_50)),
    dec_latitude = mean(dec_latitude, na.rm = TRUE),
    dec_longitude = mean(dec_longitude, na.rm = TRUE),
    elevation_m = mean(elevation_m, na.rm = TRUE),
    dbh_cm = sum(na.omit(dbh_cm)),
    pom_m = paste0(pom_m, collapse = ","),
    height_m = mean(height_m, na.rm = TRUE),
    species_binomial = first(na.omit(species_binomial)),
    family = first(na.omit(family)))

trees_list[[2]] <- stems_list$bicuar_degrad %>%
  group_by(base_stem_id) %>%
  summarise(
    plotcode = first(na.omit(plotcode)),
    stem_id = paste0(stem_id, collapse = ","),
    dec_latitude = mean(dec_latitude, na.rm = TRUE),
    dec_longitude = mean(dec_longitude, na.rm = TRUE),
    dbh_cm = sum(na.omit(dbh_cm)),
    pom_m = paste0(pom_m, collapse = ","),
    height_m = mean(height_m, na.rm = TRUE),
    species_binomial = first(na.omit(species_binomial)),
    family = first(na.omit(family)))

trees_list[[3]] <- stems_list$drc %>% 
  group_by(stem_id, plotcode) %>% 
  summarise(
    tag_id = first(na.omit(tag_id)),
    dbh_cm = sum(dbh_cm, na.rm = TRUE),
    pom_m = paste0(pom_m, collapse = ","),
    height_m = mean(height_m, na.rm = TRUE),
    species_binomial = first(na.omit(species_binomial)),
    family = first(na.omit(family))) %>%
  mutate(base_stem_id = stem_id) %>% 
  ungroup() %>%
  dplyr::select(-tag_id, -stem_id)

trees_list[[4]] <- stems_list$kilwa %>%
  group_by(plotcode, multiple) %>% 
  summarise(
    stem_id = first(na.omit(stem_id)),
    dbh_cm = sum(dbh_cm, na.rm = TRUE),
    pom_m = paste0(pom_m, collapse = ","),
    height_m = mean(height_m, na.rm = TRUE),
    species_binomial = first(na.omit(species_binomial)),
    family = first(na.omit(family))) %>%
  mutate(base_stem_id = stem_id) %>%
  ungroup() %>%
  dplyr::select(-multiple, -stem_id)

trees_list[[5]] <- stems_list$nham %>%
  group_by(multiple, plotcode) %>% 
  summarise(
    stem_id = first(na.omit(stem_id)),
    dbh_cm = sum(dbh_cm, na.rm = TRUE),
    pom_m = paste0(pom_m, collapse = ","),
    height_m = mean(height_m, na.rm = TRUE),
    species_binomial = first(na.omit(species_binomial)),
    family = first(na.omit(family))) %>%
  mutate(base_stem_id = stem_id) %>%
  ungroup() %>%
  dplyr::select(-stem_id)

# Save trees list 
names(trees_list) <- names(stems_list)
saveRDS(trees_list, "data/trees_list.rds")

# Create tree abundance matrices
trees_ab_mat_list <- lapply(trees_list, stem_ab, id = plotcode, sp = species_binomial)
names(trees_ab_mat_list) <- names(stems_list)
saveRDS(trees_ab_mat_list, "data/trees_ab_mat.rds")

