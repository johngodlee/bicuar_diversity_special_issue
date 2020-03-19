# Map of Bicuar plots
# John Godlee (johngodlee@gmail.com)
# 2019_12_20

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(maps)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(rgdal)
library(ggmap)
library(ggsn)

source("scripts/functions.R")

# Import data ----
bicuar_poly <- readOGR(dsn = "data/bicuar_shp", 
  layer = "WDPA_Mar2018_protected_area_350-shapefile-polygons")
plot_loc_list <- readRDS("data/plot_loc_list.rds")

bicuar_df <- fortify(bicuar_poly)

plot_loc_df <- do.call(rbind,  plot_loc_list[c("bicuar", "bicuar_degrad")])
plot_loc_df$group <- case_when(grepl("ABGD", plot_loc_df$plotcode) ~ "degrad",
  TRUE ~ "intact")

plot_bicuar_df <- plot_loc_df %>%
  filter(grepl("ABG", plotcode)) %>%
  dplyr::select(plotcode, group, dec_longitude, dec_latitude) %>%
  mutate(plotcode = gsub("ABG", "", .$plotcode) %>%
      gsub("-0*", "", .))

bicuar_centre <- c(14.8, -15.25)

bicuar_plot_map <- ggplot() +
  geom_polygon(data = bicuar_df, aes(x = long, y = lat), 
    fill = NA, colour = "black", show.legend = FALSE) + 
  geom_point(data = plot_bicuar_df, 
    aes(x = dec_longitude, y = dec_latitude, colour = group)) +
  geom_label_repel(data = plot_bicuar_df, 
    aes(x = dec_longitude, y = dec_latitude, label = plotcode),
    label.padding = 0.1) +
  coord_equal() +
  theme_bw() + 
  labs(x = "Longitude", y = "Latitude") + 
  geom_point(aes(x = bicuar_centre[1], y = bicuar_centre[2]))

# register_google(key = "", write = TRUE)
# bicuar_ggmap <- get_googlemap(bicuar_centre, zoom = 9, maptype = "satellite")
# saveRDS(bicuar_ggmap, "data/bicuar_ggtiles.rds")
bicuar_ggmap <- readRDS("data/bicuar_ggtiles.rds")

bicuar_ggmap_plot <- ggmap(bicuar_ggmap) +
  geom_polygon(data = bicuar_df, aes(x = long, y = lat), 
    fill = NA, colour = "#CF44D4", size = 1, show.legend = FALSE) + 
  geom_point(data = plot_bicuar_df, 
    aes(x = dec_longitude, y = dec_latitude, fill = group), 
    colour = "black", shape = 21, size = 3) +
  # geom_label_repel(data = plot_bicuar_df, 
  #   aes(x = dec_longitude, y = dec_latitude, label = plotcode, colour = group),
  #   size = 4,
  #   segment.colour = "black",  min.segment.length = 0,
  #   label.padding = 0.075, point.padding = 0.5, box.padding = 0.6) +
  theme_bw() + 
  theme(legend.position = "right") +
  labs(x = "Longitude", y = "Latitude")  + 
  lims(x = c(14.2, 15.5), y = c(-15.75, -14.75)) + 
  scale_fill_manual(name = "", labels = c("Disturbed", "Not\ndisturbed"), values = degrad_pal) + 
  scale_colour_manual(name = "", labels = c("Disturbed", "Not\ndisturbed"), values = degrad_pal, 
    guide = FALSE) + 
  scalebar(x.min = 14.3, x.max = 14.7,
    y.min = -15.725,  y.max = -13, st.dist = 0.009, 
    transform = TRUE, model = "WGS84", dist = 20, dist_unit = "km",
    height = 0.01,
    box.fill = c("black", "white"), st.color = "black", box.color = "grey")

pdf("img/bicuar_map.pdf", width = 8, height = 6)
bicuar_ggmap_plot
dev.off()

