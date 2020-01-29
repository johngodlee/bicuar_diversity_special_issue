# Get plot locations for each plot
# John Godlee (johngodle#gmail.com)
# 2019_12_20

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages 
library(dplyr)

# Import data ----
stems_list <- readRDS("data/stems_list.rds")

# Isolate plot locations ----
plot_loc_list <- lapply(stems_list, function(x){
  x %>% 
    group_by(plotcode) %>%
    summarise(dec_longitude = mean(dec_longitude, na.rm = TRUE),
      dec_latitude = mean(dec_latitude, na.rm = TRUE))
})

# Save
saveRDS(plot_loc_list, "data/plot_loc_list.rds")
