# Clean stem data before analysis
# John Godlee (johngodlee@gmail.com)
# 2019_11_28

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages 
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
library(tibble)

source("scripts/functions.R")

# Import data ----

# Raw stems data
bicuar_stems <- read.csv("data/bicuar_stems.csv")
bicuar_degrad_stems <- read.csv("data/bicuar_degrad_stems.csv")
load("data/seosaw_plot_summary5Apr2019.Rdata")
load("data/seosaw_stems.Rdata")

# Clean bicuar data ----

# Remove stems with no DBH value, no species, add plotcode, remove dead stems
bicuar_stems_clean <- bicuar_stems %>%
  filter(!is.na(species_binomial),
    !is.na(dbh_cm),
    alive_dead == "A") %>%
  mutate(plotcode = paste0("ABG-", str_pad(plot, 3, pad = "0")),
    dbh_cm = POMadj(.$dbh_cm, .$pom_m))
  
bicuar_degrad_stems_clean <- bicuar_degrad_stems %>%
  filter(!is.na(species_binomial),
    !is.na(dbh_cm),
    alive_dead == "A") %>%
  mutate(plotcode = paste0("ABGD-", str_pad(plot, 3, pad = "0")),
    dbh_cm = POMadj(.$dbh_cm, .$pom))

# Clean seosaw_data ----

# Flatten lists
ssaw8_flat <- c(ssaw8[sapply(ssaw8, is.data.frame)], ssaw8$floristicsSz)

# Exclude our Bicuar plots from 2018
ssaw8_nobnp <- lapply(ssaw8_flat, function(x){
  x %>% filter(!plotcode %in% c("ABG-001", "ABG-002", "ABG-003", "ABG-004"))
})

# Which plots should we use?
ssaw8_nobnp$plotInfoFull %>%
  arrange(desc(area_of_plot)) %>%
  print(n = 50)

ssaw8_nobnp$plotInfoFull %>%
  filter(area_of_plot == 1) %>%
  print(n = 50)

# Setup of DRC big plot
pdf("img/dks001_subplot_map.pdf", width = 10, height = 10)
s %>%
  filter(plotcode == "DKS001") %>%
  mutate(subplot1 =  factor(subplot1, levels = 1:40)) %>%
  ggplot(., aes(x = x , y = y, colour = subplot1, shape = subplot1)) + 
  geom_point() + 
  scale_shape_manual(values = rep(0:9, times = 4)) +
  coord_equal() + 
  theme_bw()
dev.off()

##' Kilwa - 25
##' Nhambita - 15
##' DRC - Subdivide into 50x50 m parcels (1:4, 5:9, 10:13 ...) in clean_stem_data and 
##'   combine four of them to a 1 ha plot, maximum of 10 plots + 
##'   2 separate 1 ha plots = 12

# Assemble DRC plot stem data ----
s_drc_subs <- s %>% 
  filter(plotcode == "DKS001") %>%
  mutate(subplot = cut(as.integer(.$subplot1), 
    breaks = seq(from = 1, to = 45, by = 4), 
    labels = FALSE, right = FALSE)) %>%
  mutate(plotcode = paste0(plotcode, "_", subplot)) %>%
  dplyr::select(-subplot)

s_drc <- s %>% 
  filter(plotcode %in% c("DKS002", "DKS003")) %>%
  bind_rows(., s_drc_subs) %>%
  group_by(stem_id, plotcode) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  mutate(plotcode = gsub("DKS", "DKS-", .$plotcode)) %>%
  mutate(multiple = gsub("^.*_", "", .$stem_id)) %>%
  mutate(stem_id = gsub("^(.*)[_].*", "\\1", .$stem_id))

# Assemble Kilwa plot stem data ----

# Exclude kilwa plots which are weird in species composition 
# Create stem abundance matrix
kilwa_ab_mat <- s %>%
  filter(grepl("TKW", plotcode)) %>%
  mutate(species_binomial = paste0(genus, "_", species)) %>%
  filter(!is.na(species_binomial), !is.na(plotcode)) %>%
  group_by(plotcode, species_binomial, .drop = FALSE) %>%
  tally() %>%
  spread(species_binomial, n, fill = 0) %>%
  ungroup() %>%
  as.data.frame() %>%
  dplyr::select(-plotcode) %>%
  filter_all(any_vars(. != 0))

rownames(kilwa_ab_mat) <- s %>% 
  filter(grepl("TKW", plotcode)) %>% 
  pull(plotcode) %>% 
  unique()

# Perform NMDS
##' Use Bray Curtis distance as it's not affected by null values in the matrix
kilwa_ab_nmds <- metaMDS(kilwa_ab_mat, distance = "bray", 
  try = 100)

# Extract site (plot) scores from NMDS analysis
kilwa_plot_scores <- as.data.frame(scores(kilwa_ab_nmds))  
kilwa_plot_scores$plot <- rownames(kilwa_plot_scores)

# Extract species scores from NMDS analysis
kilwa_species_scores <- as.data.frame(scores(kilwa_ab_nmds, "species")) 
kilwa_species_scores$species_binomial <- rownames(kilwa_species_scores)

# Plot extracted scores in ggplot2
ggplot() + 
  geom_point(data = kilwa_species_scores,
    aes(x = NMDS1, y = NMDS2)) + 
  geom_label(data = kilwa_plot_scores,
    aes(x = NMDS1, y = NMDS2, label = plot, colour = plot)) + 
  coord_equal() +
  theme_classic() + 
  theme(legend.position = "none") 

##' Exclude Kilwa TKW-002, TKW-003, TKW-004
kilwa_chosen <- c(
  "TKW-001", "TKW-005", "TKW-006", "TKW-007", 
  "TKW-008", "TKW-009", "TKW-010", "TKW-011", 
  "TKW-012", "TKW-013", "TKW-014", "TKW-015", 
  "TKW-016", "TKW-017", "TKW-018", "TKW-019", 
  "TKW-020", "TKW-021", "TKW-022", "TKW-023", "TKW-024", "TKW-025") 

s_kilwa <- s %>%
  filter(plotcode %in% kilwa_chosen) %>%
  group_by(plotcode, stem_id) %>%
  filter(year == max(year))

# Assemble Nhambita plot stem data ----

s_nham <- s %>%
  filter(grepl("MGR", plotcode)) %>%
  group_by(plotcode, tag_id) %>%
  filter(year == max(year))

# Combine all stem data dataframes ----

s_seosaw <- bind_rows(s_drc, s_kilwa, s_nham)

s_seosaw$species_binomial <- paste0(s_seosaw$genus, " ", s_seosaw$species)

s_seosaw$group <- case_when(
  grepl("TKW", s_seosaw$plotcode) ~ "kilwa",
  grepl("DKS", s_seosaw$plotcode) ~ "drc",
  grepl("MGR", s_seosaw$plotcode) ~ "nham",
  TRUE ~ "seosaw")

seosaw_stems_clean <- s_seosaw %>%
  mutate(species_binomial = case_when(
    is.na(species_binomial) ~ paste(spp_local, "local", sep = " "),
    species_binomial == "Indet indet" ~ paste(spp_local, "local", sep = " "),
    TRUE ~ species_binomial)) %>%
  filter(!is.na(species_binomial),
    !is.na(diam),
    alive == "A" | is.na(alive))

seosaw_stems_clean$pom_m <- seosaw_stems_clean$pom / 100

# Add approximate stem locations
seosaw_plot_loc <- data.frame(plotcode = ssaw8$plotInfoFull$plotcode, 
  dec_latitude = ssaw8$plotInfoFull$latitude_of_centre,
  dec_longitude = ssaw8$plotInfoFull$longitude_of_centre)

seosaw_plot_loc$plotcode <- gsub("DKS", "DKS-", seosaw_plot_loc$plotcode)
dks001_loc <- seosaw_plot_loc[rep(grep("DKS-001", seosaw_plot_loc$plotcode), each = 10),] 
dks001_loc$plotcode <- paste0(dks001_loc$plotcode, "_", seq(from = 1, to = 10))

seosaw_plot_loc <- seosaw_plot_loc %>% filter(plotcode != "DKS-001") %>%
  bind_rows(.,dks001_loc)

seosaw_stems_clean <- left_join(seosaw_stems_clean, seosaw_plot_loc, by = "plotcode")

# Split by group 
seosaw_stems_clean_list <- split(seosaw_stems_clean, seosaw_stems_clean$group)

# Compile data to one list of tidy dataframes ----
bicuar_stems_clean <- bicuar_stems_clean %>%
  dplyr::select(plotcode, stem_id, base_stem_id, subplot_20_50, subplot_20_100, 
    dec_latitude, dec_longitude, elevation_m,
    dbh_cm, pom_m, height_m, species_binomial)

bicuar_degrad_stems_clean <- bicuar_degrad_stems_clean %>%
  dplyr::select(plotcode, stem_id, base_stem_id, dec_latitude, dec_longitude,
    dbh_cm, pom_m = pom, height_m, species_binomial)

seosaw_stems_clean_list <- lapply(seosaw_stems_clean_list, function(x){
  x %>%
    dplyr::select(plotcode, stem_id, tag_id, multiple, dbh_cm = diam130, pom_m, 
      height_m = height, dec_latitude, dec_longitude, species_binomial)
})

stems_list <- c(list(bicuar = bicuar_stems_clean,
  bicuar_degrad = bicuar_degrad_stems_clean), 
  seosaw_stems_clean_list)

# Save data ----
saveRDS(stems_list, "data/stems_list.rds")


