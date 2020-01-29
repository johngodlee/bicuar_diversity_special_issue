# Table of species in Bicuar with information
# John Godlee (johngodlee@gmail.com)
# 2019_12_20

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(dplyr)
library(stargazer)

source("scripts/functions.R")

# Import data ----
stems_list <- readRDS("data/stems_list.rds")
trees_list <- readRDS("data/trees_list.rds")
tree_family_lookup <- readRDS("data/tree_family_lookup.rds")
bicuar_unique <- readRDS("data/bicuar_unique.rds")

# Create columns for table ----
# Mean DBH
dbh <- stems_list$bicuar %>%
  group_by(family, species_binomial) %>%
  summarise(median_dbh = median(dbh_cm, na.rm = TRUE),
    samples_dbh = n(),
    sd_dbh = sd(dbh_cm, na.rm = TRUE),
    sem_dbh = 1.2533 * (sd_dbh / sqrt(samples_dbh)))

basal <- stems_list$bicuar %>%
  group_by(species_binomial, plotcode) %>%
  mutate(basal_area = ((dbh_cm / 2)^2) * 0.0001) %>%
  summarise(sum_basal_area = sum(basal_area, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species_binomial) %>%
  summarise(mean_basal_area = mean(sum_basal_area, na.rm = TRUE),
    samples_basal_area = n(),
    sd_basal_area = sd(sum_basal_area, na.rm = TRUE),
    sem_basal_area = sd_basal_area / sqrt(samples_basal_area))

freq <- trees_list$bicuar %>%
  group_by(species_binomial, plotcode) %>%
  summarise(freq = n()) %>%
  ungroup() %>%
  group_by(species_binomial) %>%
  summarise(total_freq = sum(freq, na.rm = TRUE),
    mean_freq = mean(freq, na.rm = TRUE),
    sd_freq = sd(freq, na.rm = TRUE),
    samples_freq = n(),
    sem_freq = sd_freq / sqrt(samples_freq))

# Create dataframe ----
species_df <- left_join(dbh, basal, by = "species_binomial") %>%
  left_join(.,freq, by = "species_binomial") %>%
  arrange(species_binomial)

# Create table for latex ----
ba_format <- function(x){
  dplyr::case_when(x < 0.01 ~ "<0.01",
    TRUE ~ sprintf("%.2f", round(x, digits = 2)))
}

species_df_out <- species_df %>%
  ungroup() %>%
  mutate(dbh = paste0(sprintf("%.1f", round(median_dbh, 1)), 
    "(", sprintf("%.2f", round(sem_dbh, 2)), ")"),
    basal_area =  paste0(ba_format(mean_basal_area), 
      "(", sprintf("%.3f", round(sem_basal_area, 2)), ")"),
    mean_freq = paste0(round(mean_freq, 1),
      "(", sprintf("%.2f", round(sem_freq, 2)), ")"),
    family = as.character(family)) %>%
  dplyr::select(family, species_binomial, dbh, basal_area, total_freq, mean_freq)

species_df_out$species_binomial <- case_when(
  species_df_out$species_binomial %in% bicuar_unique ~ paste0("\\textbf{*", species_df_out$species_binomial, "}"),
  TRUE ~ species_df_out$species_binomial
)

fileConn <- file(paste0("include/bicuar_species.tex"))
writeLines(
  stargazer(species_df_out, summary = FALSE, label = "bicuar_species", 
    digit.separate = 0, rownames = FALSE),
  fileConn)
close(fileConn)

