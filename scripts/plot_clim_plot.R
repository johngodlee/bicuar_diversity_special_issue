# Make maps of locations of plots in Africa and climate space plot
# John Godlee (johgodlee@gmail.com)
# 2019_12_20

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(ggplot2)
library(ggnewscale)
library(maps)
library(rgdal)
library(dplyr)

source("scripts/functions.R")

# Import data ----
plot_loc_list <- readRDS("data/plot_loc_list.rds")
plot_clim_list <- readRDS("data/plot_clim_list.rds")
miombo_clim <- readRDS("data/region_temp_precip.rds")
precip_raster <- readRDS("data/precip_raster.rds")

# Whites veg map
white_veg <- readOGR(dsn="data/whitesveg", 
  layer="Whites vegetation")

white_veg_fort <- fortify(white_veg, region = "DESCRIPTIO")
names(white_veg_fort)
length(unique(white_veg_fort$id))

white_veg_miombo <- white_veg_fort %>%
  filter(id %in% c("Moist-infertile savanna"),
    lat < -2)

# Plot of MAP~MAT and plots ----
# Subset to big plots
plot_clim_list_big <- plot_clim_list[c(1, 3:5)]

# Create dataframes
plot_clim_names <- mapply(function(x,y){ 
  rep(y, length(x[,1])) 
  }, x = plot_clim_list_big, y = c("Angola", "DRC", "Tanzania", "Mozambique"))

plot_clim_list_big_group <- lapply(1:length(plot_clim_list_big), function(x){
  cbind(plot_clim_list_big[[x]], group = plot_clim_names[[x]])
})

plot_clim_df <- do.call(rbind, plot_clim_list_big_group)

# Create plot
pdf(file = "img/temp_precip.pdf", width = 6, height = 6)
ggplot() + 
  stat_binhex(data = miombo_clim, 
    mapping = aes(x = t_vals, y = p_vals, colour = ..count.., fill = ..count..),
    bins = 500) +
  scale_fill_gradient(name = "Density", low = "#D9D9D9", high = "black",
    trans = "log", breaks = c(1, 5, 10, 50, 100, 200, 400)) + 
  scale_colour_gradient(name = "Density", low = "#D9D9D9", high = "black",
    trans = "log", breaks = c(1, 5, 10, 50, 100, 200, 400)) + 
  new_scale_fill() +
  geom_point(data = plot_clim_df,
    mapping = aes(x = mat, y = map, fill = group), 
    colour = "black", shape = 21, size = 4) + 
  scale_fill_manual(name = "", values = big_pal) +
  theme_classic() + 
  labs(x= expression("MAT" ~ (degree*C)), 
    y = expression("MAP" ~ (mm ~ y^-1)))
dev.off()

# Map of plots in miombo ----
s_af <- iso.expand(c("ZAF", "COD", "NAM", "ZMB", "BWA", "ZWE", "MOZ", 
  "MWI", "AGO", "TZA", "COG", "RWA", "BDI", "UGA", "KEN"))

# Create map of country outlines
map_africa <- borders(database = "world", regions = s_af, fill = "grey90", colour = "black")
map_africa_fill <- borders(database = "world", regions = s_af, fill = "grey90")
map_africa_colour <- borders(database = "world", regions = s_af, colour = "black")

plot_loc_list_big <- plot_loc_list[c("bicuar", "drc", "kilwa", "nham")]
plot_loc_list_big_group <- lapply(1:length(plot_loc_list_big), function(x){
  cbind(plot_loc_list_big[[x]], group = plot_clim_names[[x]])
})

plot_loc_df <- do.call(rbind, plot_loc_list_big_group)

# Create simplified spatial object from precipitation raster
# precip_raster_agg <- aggregate(precip_raster, fact = 10)
# precip_spdf <- as(precip_raster_agg, "SpatialPixelsDataFrame")
# precip_df <- as.data.frame(precip_spdf)
# colnames(precip_df) <- c("value", "x", "y")
# saveRDS(precip_df, file = "data/precip_df.rds")
precip_df <- readRDS("data/precip_df.rds")

# Make plot
plot_map <- ggplot() + 
  geom_tile(data = precip_df, aes(x = x, y = y, fill = value, colour = value)) + 
  scale_colour_gradient(name = expression("MAP"~(mm~y^-1)),
    low = "#FAE987", high = "#332B00") + 
  scale_fill_gradient(name = expression("MAP"~(mm~y^-1)),
    low = "#FAE987", high = "#332B00") + 
  new_scale_fill() + 
  map_africa_colour +
  geom_point(data = plot_loc_df, 
    aes(x = dec_longitude, y = dec_latitude, fill = group), 
    colour = "black", shape = 21, size = 4, alpha = 1, show.legend = FALSE) +
  scale_fill_manual(name = "", values = big_pal, 
    labels = c("Angola", "DRC", "Tanzania", "Mozambique")) +
  coord_map() + 
  ylim(-35.5, 10) + 
  labs(x = "Longitude", y = "Latitude") + 
  theme_classic() + 
  theme(legend.position = c(0.89, 0.15), 
    legend.background = element_blank())

pdf(file = "img/plot_map.pdf", plot_map, width = 6, height = 6)
plot_map
dev.off()

