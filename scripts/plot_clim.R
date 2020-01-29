# Make a dataset of environmental variables for chosen plots
# John Godlee (johngodlee@gmail.com)
# 2019_12_01

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(raster)
library(rgdal)
library(dplyr)

# Import data ----
load("data/seosaw_plot_summary5Apr2019.Rdata")
plot_loc_list <- readRDS("data/plot_loc_list.rds")

rastlist <- list.files(path = "data/wc2.0_2.5m_bio/", 
  pattern = '.tif$', 
  all.files = TRUE, 
  full.names = TRUE)

allrasters <- lapply(rastlist, raster)

cwd <- raster("data/CWD.tif")

# Plot locations ----
plot_loc_spdf_list <- lapply(plot_loc_list, function(x){
  SpatialPointsDataFrame(
    coords = data.frame(x$dec_longitude, x$dec_latitude),
    data = x)
})

# For each raster in the stack, take pixel value ----
plot_vals_list <- lapply(plot_loc_spdf_list, function(x){
  lapply(allrasters, function(y){
    raster::extract(y, x)
  })
})

cwd_list <- lapply(plot_loc_spdf_list, function(x){
  raster::extract(cwd, x)
})

# Combine into a tidy dataframe ----
plot_vals_df_list <- lapply(1:length(plot_loc_spdf_list), function(x){
  df <- cbind(plot_loc_spdf_list[[x]]@data, as.data.frame(do.call(cbind, plot_vals_list[[x]])))
  df$cwd <- cwd_list[[x]]
  return(df)
  })

clim_df_names <- c("plotcode", "dec_longitude", "dec_latitude", 
  "mat", "temp_diurnal_range", "isothermality", "mat_sd", "max_temp", 
  "min_temp", "temp_annual_range", "mean_temp_wet_q", "mean_temp_dry_q",
  "mean_temp_warm_q", "mean_temp_cold_q", "map", "precip_wet_m", "precip_dry_m",
  "map_sd", "precip_wet_q", "precip_dry_q", "precip_warm_q", "precip_cold_q", "cwd")

plot_vals_df_clean_list <- lapply(plot_vals_df_list, setNames, nm = clim_df_names)

# Write combined dataset  ----
saveRDS(plot_vals_df_clean_list, "data/plot_clim_list.rds")
