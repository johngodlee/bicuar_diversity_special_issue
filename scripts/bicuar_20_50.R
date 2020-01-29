# Create 20x50 m subplots in Bicuar 2018 plot data
# John Godlee (johngodlee@gmail.com)
# 2019_11_29

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(rgdal)
library(geosphere)
library(dplyr)

# Define a CRS
epsg <- make_EPSG()  # Create list of EPSG dataset to search for CRS
wgs84 <- epsg[grep("WGS 84", epsg$note, ignore.case=TRUE),]  # Search for wgs84
wgs84[grep("longlat", wgs84$prj4, ignore.case=TRUE),]  # grep proj4string to check
wgs84_crs <- CRS(wgs84[grep("longlat", wgs84$prj4, ignore.case=TRUE),]$prj4[2])  # Store string as vector

# Import data ----
stems_list <- readRDS("data/stems_list.rds")
corners_bicuar <- read.csv("data/bicuar_plot_corner_loc.csv")

# Subset data to 2018 plots ----
s_bicuar_2018 <- stems_list[[1]][stems_list[[1]]$plotcode %in% c("ABG-001", "ABG-002", "ABG-003", "ABG-004"),]
corners_bicuar_2018 <- corners_bicuar[corners_bicuar$plot %in% c("ABG-001", "ABG-002", "ABG-003", "ABG-004"),]

# Convert to stems spatialpoints dataframe ----
s_bicuar_nona <- s_bicuar_2018[!is.na(s_bicuar_2018$dec_longitude) & 
    !is.na(s_bicuar_2018$dec_latitude),]

s_bicuar_points <- SpatialPointsDataFrame(s_bicuar_nona[,c("dec_longitude", "dec_latitude")],  # extract only long:lat coords
  proj4string = wgs84_crs,
  data = data.frame(s_bicuar_nona[,"plotcode"]))

# Find midpoint along north and south edge to construct a polygon to split the plot
mid_poly_list <- list()
mid_points_list <- list()

# Create polygons which split the plot in half  perpendicular to the 20x100 strips
for(i in 1:4){
  corners <- data.frame(corners_bicuar_2018[corners_bicuar_2018$plot == paste0("ABG-00", i), c("plotcode", "dec_longitude", "dec_latitude")])
  colnames(corners) <- c("id", "x", "y")
  
  corners_points <- SpatialPointsDataFrame(corners[,2:3], 
    proj4string = wgs84_crs, 
    data = data.frame(corners[,1]))
  
  mid_n <- midPoint(corners_points@coords[1,], corners_points@coords[2,])
  mid_s <- midPoint(corners_points@coords[3,], corners_points@coords[4,])
  
  mid_poly_w <- Polygon(rbind(corners_points@coords[1,], mid_n, 
    mid_s, corners_points@coords[4,], corners_points@coords[1,]))
  mid_poly_w_sp <- Polygons(list(mid_poly_w), "w")
  
  mid_poly_e <- Polygon(rbind(mid_n, corners_points@coords[2,], 
    corners_points@coords[3,], mid_s, mid_n))
  mid_poly_e_sp <- Polygons(list(mid_poly_e), "e")
  
  mid_poly_sps = SpatialPolygons(list(mid_poly_w_sp, mid_poly_e_sp), proj4string = wgs84_crs)

  mid_poly_list[[i]] <- mid_poly_sps
  
  mid_points_list[[i]] <- over(s_bicuar_points , mid_poly_list[[i]] , fn = NULL) 
}

# Combine matches from `over()` into one vector
s_bicuar_2018$w_e <- pmin(mid_points_list[[1]], mid_points_list[[2]], 
  mid_points_list[[3]], mid_points_list[[4]],
  na.rm = TRUE)

# Format the subplot matches 
s_bicuar_2018_subs <- s_bicuar_2018 %>%
  filter(!is.na(w_e)) %>%
  mutate(subplot_20_50 = paste0(subplot_20_100, "_", w_e)) %>%
  mutate(subplot_20_50 = case_when(
    subplot_20_50 == "1_1" ~ 1,
    subplot_20_50 == "1_2" ~ 2,
    subplot_20_50 == "2_1" ~ 3,
    subplot_20_50 == "2_2" ~ 4,
    subplot_20_50 == "3_1" ~ 5,
    subplot_20_50 == "3_2" ~ 6,
    subplot_20_50 == "4_1" ~ 7,
    subplot_20_50 == "4_2" ~ 8,
    subplot_20_50 == "5_1" ~ 9,
    subplot_20_50 == "5_2" ~ 10)) %>%
  dplyr::select(stem_id, subplot_20_50)

# Join back into main dataset
stems_list[[1]] <- left_join(stems_list[[1]], s_bicuar_2018_subs, by = "stem_id")

stems_list[[1]]$subplot_20_50 <- pmin(stems_list[[1]]$subplot_20_50.x,
  stems_list[[1]]$subplot_20_50.y, na.rm = TRUE)

stems_list[[1]] <- stems_list[[1]] %>%
  dplyr::select(-subplot_20_50.x, -subplot_20_50.y)

# Write to csv
saveRDS(stems_list, "data/stems_list.rds")

