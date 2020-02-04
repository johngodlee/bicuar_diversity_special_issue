# Generate a dataset of precipitation and temperature means across the Miombo region
# John Godlee (johngodlee@gmail.com)
# 2019_10_02

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(raster)
library(rgdal)

# Import data ----
# Get Vegetation type shapefile
miombo <- readOGR(dsn="data/teow_miombo/miombo", layer="miombo")

# Temperature raster
# Get list of files
rastlist_t <- list.files(path = "data/wc2.0_30s_tavg", 
  pattern = '.tif$', 
  all.files = TRUE, 
  full.names = TRUE)

# Import
allrasters_t <- lapply(rastlist_t, raster)

# Crop to miombo extent
allrasters_t_crop <- lapply(allrasters_t, function(x){
  crop(x, miombo)
})

# Stack
allrasters_t_crop_stack <- raster::stack(allrasters_t_crop)

# Take mean of all in stack
allrasters_t_mean <- calc(allrasters_t_crop_stack, mean, na.rm = TRUE)

# Crop to miombo outline
allrasters_t_mean_crop_mask <- mask(allrasters_t_mean, miombo)

# Extract all values
t_vals <- values(allrasters_t_mean_crop_mask)

# Precipitation raster
# Get list of files
rastlist_p <- list.files(path = "data/wc2.0_30s_prec", 
  pattern = '.tif$', 
  all.files = TRUE, 
  full.names = TRUE)

# Import
allrasters_p <- lapply(rastlist_p, raster)

# Crop to plot extent
allrasters_p_crop <- lapply(allrasters_p, function(x){
  crop(x, miombo)
})

# Stack
allrasters_p_crop_stack <- raster::stack(allrasters_p_crop)

# Get sum of all in stack
allrasters_p_mean <- calc(allrasters_p_crop_stack, sum, na.rm = TRUE)

# Crop to miombo outline
allrasters_p_mean_crop_mask <- mask(allrasters_p_mean, miombo)

# Extract all values
p_vals <- values(allrasters_p_mean_crop_mask)

# Combine to dataframe
t_p <- data.frame(t_vals, p_vals)

saveRDS(t_p, "data/region_temp_precip.rds")

# Save precipitation raster
saveRDS(allrasters_p_mean_crop_mask, "data/precip_raster.rds")

