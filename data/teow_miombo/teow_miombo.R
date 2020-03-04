library(rgdal)
library(dplyr)
library(ggplot2)

test <- readOGR(dsn = ".", layer = "wwf_terr_ecos")

miombo <- test[grepl("miombo", test$ECO_NAME, ignore.case = TRUE),]

writeOGR(miombo, dsn = "miombo", layer = "miombo", driver = "ESRI Shapefile")

