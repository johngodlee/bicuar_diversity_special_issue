# Descriptive statistics for big plots
# John Godlee (johngodlee@gmail.com)
# 2019_12_20

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(taxize)
library(dplyr)
library(stargazer)

source("scripts/functions.R")

# Import data ----
trees_list <- readRDS("data/trees_list.rds")
plot_clim_list <- readRDS("data/plot_clim_list.rds")
trees_ab_mat_list <- readRDS("data/trees_ab_mat.rds")
plot_div_list <- readRDS("data/plot_div_list.rds")


# Total number of plots sampled ----
trees_list_big <-  trees_list[c("bicuar", "drc", "kilwa", "nham")]
trees_list_big_clean <- lapply(trees_list_big, function(x){
  x %>%
    dplyr::select(plotcode, dbh_cm, species_binomial, family)
})

n_plots <- lapply(trees_list_big_clean, function(x){
  length(unique(x$plotcode))
})
total_n_plots <- do.call(sum, n_plots)
bicuar_n_plots <- length(unique(trees_list_big_clean[['bicuar']]$plotcode))

# Total number of trees sampled ----
n_trees <- lapply(trees_list_big_clean, function(x){
  length(x$plotcode)
})
total_n_trees <- do.call(sum, n_trees)
bicuar_n_trees <- length(trees_list_big_clean[['bicuar']]$plotcode)

# Total number of families ----
group_n_families <- lapply(trees_list_big_clean, function(x){
  length(unique(x$family))
})

total_n_families <- length(unique(do.call(rbind, trees_list_big_clean)$family))
bicuar_n_families <- group_n_families[[1]]

# Total number of species ----
group_n_species <- lapply(trees_list_big_clean, function(x){
  length(unique(x$species_binomial[!is.na(x$species_binomial)]))
})
total_n_species <- length(unique(do.call(rbind, trees_list_big_clean)$species_binomial))
bicuar_n_species <- group_n_species[[1]]

# Most diverse families ----
family_div <- function(x){
  x %>% 
    group_by(family) %>%
    tally() %>%
    arrange(desc(n)) %>%
    pull(family) %>%
    .[1]
}

# Number of species in most diverse family (Fabaceae) ----
most_div_family_group <- lapply(trees_list_big_clean, family_div)

species_div_family <- length(do.call(rbind, trees_list_big_clean) %>%
    group_by(family, species_binomial) %>%
    summarise() %>%
    filter(family == "Fabaceae") %>%
    pull(species_binomial))

# Species found only in Bicuar ----
# Get species in each group
other_species <- do.call(rbind, trees_list_big_clean[c("drc", "kilwa", "nham")])

# Which species are found only in Bicuar? 
bicuar_unique <- setdiff(trees_list_big_clean$bicuar$species_binomial, other_species$species_binomial)
n_bicuar_unique <- length(bicuar_unique)
saveRDS(bicuar_unique, "data/bicuar_unique.rds")

# Which of these unique species are common in Bicuar?
n_bicuar_trees <- trees_list_big_clean$bicuar %>%
  group_by(species_binomial) %>%
  tally() %>%
  filter(species_binomial %in% bicuar_unique) %>%
  arrange(desc(n))

nbg <- c(n_bicuar_trees[n_bicuar_trees$species_binomial == "Brachystegia glaucescens",2])
nbp <- c(n_bicuar_trees[n_bicuar_trees$species_binomial == "Baikiaea plurijuga",2])
nbm <- c(n_bicuar_trees[n_bicuar_trees$species_binomial == "Baphia massaiensis",2])

# Table describing each group of plots ----
plot_clim <- do.call(rbind, plot_clim_list[c(1,3:5)])

plot_clim$group <- case_when(
  grepl("ABG", plot_clim$plotcode) ~ "Angola",
  grepl("DKS", plot_clim$plotcode) ~ "DRC",
  grepl("TKW", plot_clim$plotcode) ~ "Tanzania",
  grepl("MGR", plot_clim$plotcode) ~ "Mozambique"
)

group_descrip <- plot_clim %>%
  group_by(group) %>%
  summarise(mat = round(mean(mat), digits = 1),
    map = round(mean(map), digits = 0),
    cwd = round(mean(cwd), digits = 0),
    dec_latitude = round(mean(dec_latitude), digits = 2),
    dec_longitude = round(mean(dec_longitude), digits = 2),
    n_plots = n())

group_descrip_big <- filter(group_descrip, group != "bicuar_degrad")

trees_ab_mat_list_big <- trees_ab_mat_list[c('bicuar', 'drc', 'kilwa', 'nham')]

group_descrip_big$sp <- sp_big <- sapply(trees_ab_mat_list_big, function(x){
  ncol(x)
})

fileConn <- file(paste0("include/group_descrip.tex"))
writeLines(
  stargazer(group_descrip_big, summary = FALSE, label = "group_descrip", 
    digit.separate = 0, rownames = FALSE),
  fileConn)
close(fileConn)

# Output numerical figures ----
fileConn <- file("include/data_descrip_figures.tex")
writeLines(
  c(
    paste0("\\newcommand{\\nplots}{", total_n_plots, "}"),
    paste0("\\newcommand{\\nplotsbicuar}{", bicuar_n_plots, "}"),
    paste0("\\newcommand{\\nbicuartrees}{", bicuar_n_trees, "}"),
    paste0("\\newcommand{\\ntrees}{", total_n_trees, "}"),
    paste0("\\newcommand{\\nspecies}{", total_n_species, "}"),
    paste0("\\newcommand{\\nfamilies}{", total_n_families, "}"),
    paste0("\\newcommand{\\nfabaceaespecies}{", species_div_family, "}"),
    paste0("\\newcommand{\\nbicuaruniquespecies}{", n_bicuar_unique, "}"),
    paste0("\\newcommand{\\nbicuarspecies}{", bicuar_n_species, "}"),
    paste0("\\newcommand{\\nbicuarfamilies}{", bicuar_n_families, "}"),
    paste0("\\newcommand{\\nbg}{", nbg, "}"),
    paste0("\\newcommand{\\nbp}{", nbp, "}"),
    paste0("\\newcommand{\\nbm}{", nbm, "}")),
  fileConn)
close(fileConn)

