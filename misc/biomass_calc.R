# Calculate biomass per stem for plots in Bicuar National Park
# John Godlee (johngodlee@gmail.com)
# 2019_05_17

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(ggplot2)
library(dplyr)
library(BIOMASS)
library(gridExtra)
library(colortools)

source("scripts/functions.R")

# Import data ----

# Stem measurements
stems_list <- readRDS("data/stems_list.rds")

# Get wood density estimates ----
stems_list <- lapply(stems_list, function(x){
  x$genus <- gsub("([A-z]+).*", "\\1", x$species_binomial)
  x$species <- gsub("^([A-z]+) ", "", x$species_binomial)
  wood_den <- getWoodDensity(x$genus,
    x$species,
    stand = x$plotcode,
    region = "AfricaTrop")
  
  cbind(x, wood_den[,c(1, 4:7)])
})

# Calculate AGB per stem ----
stems_list <- lapply(stems_list, function(x){
  x$agb_stem_eq4 <- computeAGB(D = x$dbh_cm, 
  coord = data.frame(x$dec_longitude, x$dec_latitude),
	WD = x$meanWD)
  return(x)
})

saveRDS(stems_list, "data/stems_list.rds")