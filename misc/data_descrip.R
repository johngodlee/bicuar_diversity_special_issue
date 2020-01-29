# Descriptive statistics of diversity across all plots
# John Godlee (johngodlee@gmail.com)
# 2019_11_29

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(dplyr)
library(taxize)
library(vegan)
library(ggplot2)
library(ggnewscale)
library(stargazer)
library(maps)
library(rgdal)
library(tidyr)
library(ggrepel)
library(ggmap)

source("scripts/functions.R")

# Import data ----
trees_list <- readRDS("data/trees_list.rds")
stems_list <- readRDS("data/stems_list.rds")
trees_ab_mat_list <- readRDS("data/trees_ab_mat.rds")
plot_clim <- readRDS("data/plot_clim_list.rds")
miombo_clim <- readRDS("data/region_temp_precip.rds")
bicuar_poly <- readOGR(dsn = "data/bicuar_shp", 
  layer = "WDPA_Mar2018_protected_area_350-shapefile-polygons")

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
# Add plot group column
plot_clim$group <- case_when(
  grepl("ABGD", plot_clim$plotcode) ~ "bicuar_degrad",
  grepl("ABG", plot_clim$plotcode) ~ "Angola",
  grepl("TKW", plot_clim$plotcode) ~ "Tanzania",
  grepl("DKS", plot_clim$plotcode) ~ "DRC",
  grepl("MGR", plot_clim$plotcode) ~ "Mozambique",
  TRUE ~ "seosaw")

# Create plot of MAP~MAT
pdf(file = "img/temp_precip.pdf", width = 6, height = 6)
ggplot() + 
  stat_binhex(data = miombo_clim, 
    mapping = aes(x = t_vals, y = p_vals, colour = ..count.., fill = ..count..),
    bins = 500) +
  scale_fill_gradient(name = "Density", low = "#D9D9D9", high = "black",
    trans = "log", breaks = c(1, 5, 10, 50, 100, 200, 400)) + 
  scale_colour_gradient(name = "Density", low = "#D9D9D9", high = "black",
    trans = "log", breaks = c(1, 5, 10, 50, 100, 200, 400)) + 
  # scale_fill_continuous(name = "Density", type = "viridis", trans = "log",
  #   breaks = c(1, 5, 10, 50, 100, 200, 400), alpha = 0.8) + 
 # scale_colour_continuous(name = "Density", type = "viridis", trans = "log", 
 #   breaks = c(1, 5, 10, 50, 100, 200, 400)) + 
  new_scale_fill() +
  geom_point(data = filter(plot_clim, group != "bicuar_degrad"),
    mapping = aes(x = mat, y = map, fill = group), 
    colour = "black", shape = 21, size = 4) + 
  scale_fill_discrete(name = "") +
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

plot_map <- ggplot() + 
  map_africa_fill + 
  geom_polygon(aes(x = long, y = lat, group = group), 
    fill = "#A37803", colour = "#A37803",
    data = white_veg_miombo, alpha = 1) +
  # geom_point(data = ssaw8$plotInfoFull, aes(x = longitude_of_centre, y = latitude_of_centre), alpha = 0.6) + 
  map_africa_colour +
  geom_point(data = filter(plot_clim, group != "bicuar_degrad"), 
    aes(x = dec_longitude, y = dec_latitude, fill = group), 
    colour = "black", shape = 21, size = 4, alpha = 1) +
  coord_map() + 
  ylim(-35.5, 10) + 
  labs(x = "Longitude", y = "Latitude") + 
  theme_classic() + 
  scale_fill_discrete(name = "", labels = c("Angola", "DRC", "Tanzania", "Mozambique")) +
  theme(legend.position = "none")

pdf(file = "img/plot_map.pdf", plot_map, width = 6, height = 6)
plot_map
dev.off()

# Total number of plots sampled ----
n_plots <- lapply(trees_list, function(x){
  length(unique(x$plotcode))
})

sum(n_plots[[1]], n_plots[[3]], n_plots[[4]], n_plots[[5]])

# Total number of trees sampled ----
n_trees <- lapply(trees_list, function(x){
 length(x$plotcode)
})

# Total number of families ----
# usethis::edit_r_environ()
# ENTREZ_KEY="d9cd1e651a746d82b62e4726838109259107"
# httr::set_config(httr::config(http_version = 0))
# families <- lapply(trees_list, function(x){
#   classification(unique(as.character(x$species_binomial)),
#     db = "gbif", accepted = TRUE)})
#saveRDS(families, "data/tree_family_lookup.rds")
families <- readRDS("data/tree_family_lookup.rds")

families_list <- lapply(families, taxa_df_gen)

families_n_group <- lapply(families_list, function(x){
  length(unique(x$family[!is.na(x$family)]))
})

families_total_bicuar <- rbind(families_list[[1]], families_list[[2]])
families_n_total_bicuar <- length(unique(families_total_bicuar$family[!is.na(families_total_bicuar$family)]))

families_total <- do.call(rbind, families_list)
families_n_total <- length(unique(families_total$family[!is.na(families_total$family)]))


# Total number of species ----
species_n_group <- lapply(families_list, function(x){
  length(unique(x$species[!is.na(x$species)]))
})

species_total_bicuar <- rbind(families_list[[1]], families_list[[2]])
species_n_total_bicuar <- length(unique(species_total_bicuar$species[!is.na(species_total_bicuar$species)]))

species_total <- do.call(rbind, families_list)
species_n_total <- length(unique(species_total$species[!is.na(species_total$species)]))

# Table of species in Bicuar
dbh <- stems_list$bicuar %>%
  group_by(family, species_binomial) %>%
  summarise(mean_dbh = mean(dbh_cm, na.rm = TRUE),
    sd_dbh = sd(dbh_cm, na.rm = TRUE)) 

basal <- stems_list$bicuar %>%
  group_by(species_binomial, plotcode) %>%
  mutate(basal_area = ((dbh_cm / 2)^2) * 0.0001) %>%
  summarise(sum_basal_area = sum(basal_area, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species_binomial) %>%
  summarise(mean_basal_area = mean(sum_basal_area, na.rm = TRUE),
    sd_basal_area = sd(sum_basal_area, na.rm = TRUE))

freq <- trees_list$bicuar %>%
  group_by(species_binomial, plotcode) %>%
  summarise(freq = n()) %>%
  ungroup() %>%
  group_by(species_binomial) %>%
  summarise(total_freq = sum(freq, na.rm = TRUE),
    mean_freq = mean(freq, na.rm = TRUE))

View(left_join(dbh, basal, by = "species_binomial") %>%
  left_join(.,freq, by = "species_binomial") %>%
  arrange(species_binomial))

# Most diverse families ----
family_div <- function(x){
  x %>% 
    group_by(family) %>%
    tally() %>%
    arrange(desc(n)) %>%
    pull(family) %>%
    .[1]
}

most_div_family_group <- lapply(families_list, family_div)

# Number of species in most diverse family (Fabaceae) ----
length(families_total %>%
  group_by(family, species) %>%
  summarise() %>%
  filter(family == "Fabaceae") %>%
    pull(species))

# Total basal area of sampled stems ----
stems_list <- lapply(stems_list, function(x){
  x$ba <- (x$dbh_cm / 2)^2
  return(x)
})

##' in m^2
ba_list <- lapply(stems_list, function(x){
  sum(x$ba) / 10000
})

# Species found only in Bicuar ----
# Get species in each group
species_group_list <- lapply(
  list(bicuar = trees_ab_mat_list$bicuar, drc = trees_ab_mat_list$drc, kilwa = trees_ab_mat_list$kilwa, nham = trees_ab_mat_list$nham), function(x){
    unique(unlist(apply(x[,-1], 1, function(y){
      names(y[y > 0])
    })))
  })

species_seosaw_list <- unique(do.call(c, list(species_group_list$drc, species_group_list$kilwa, species_group_list$nham)))

# Which species are found only in Bicuar? 
bicuar_unique <- setdiff(species_group_list$bicuar, species_seosaw_list)

# Which of these unique species are common in Bicuar?
sort(colSums(trees_ab_mat_list$bicuar[,bicuar_unique]), decreasing = TRUE)

# Diversity estimates ----

# Alpha diversity
div_df <- function(ab_mat){
  ab_mat_clean <- ab_mat[,-1]
  plotcode <- ab_mat[,1]
  
  alpha <- data.frame(plotcode, shannon = diversity(ab_mat_clean), 
    simpson = diversity(ab_mat_clean, index = "simpson"))
  return(alpha)
}

alpha_list <- lapply(trees_ab_mat_list, div_df)

alpha_df <- do.call(rbind, alpha_list)

alpha_df$group <- gsub("\\..*$", "", rownames(alpha_df))

alpha_df <- filter(alpha_df, group != "bicuar_degrad")

alpha_df$group <- case_when(
  grepl("bicuar", alpha_df$group) ~ "Angola",
  grepl("drc", alpha_df$group) ~ "DRC",
  grepl("kilwa", alpha_df$group) ~ "Tanzania",
  grepl("nham", alpha_df$group) ~ "Mozambique"
)

##' ISO3 country codes
alpha_df$plotcode <- gsub("ABGD-", "ANGD-", alpha_df$plotcode) %>%
        gsub("ABG-", "ANG-", .) %>%
        gsub("DKS-", "DRC-", .) %>%
        gsub("TKW-", "TZA-", .) %>%
        gsub("MGR-", "MOZ-", .) %>%
        gsub("-0*", "-", .) %>%
        gsub("DRC-2", "DRC-11", .) %>%
        gsub("DRC-3", "DRC-12", .) %>%
        gsub("DRC-1_", "DRC-", .)


alpha_df$plotcode_shannon = factor(alpha_df$plotcode, levels = alpha_df$plotcode[order(desc(alpha_df$shannon))])
alpha_df$plotcode_simpson = factor(alpha_df$plotcode, levels = alpha_df$plotcode[order(desc(alpha_df$simpson))])

shannon_by_plot <- ggplot(alpha_df, aes(x = plotcode_shannon, y = shannon)) + 
  geom_bar(stat = "identity", colour = "black", aes(fill = "grey")) + 
  labs(x = "Plot", y = "Shannon diversity index") + 
  facet_grid(. ~ group, scales = "free_x", space = "free") + 
  theme_bw() + 
  theme(axis.text.x=element_text(
    angle=45, 
    vjust=1, 
    hjust=1)) + 
  theme(legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")) + 
  scale_fill_manual(name = "DBH bin (cm)", values = "grey",
    guide = guide_legend(override.aes = list(colour = "white", fill = "white")))

simpson_by_plot <- ggplot(alpha_df, aes(x = plotcode_simpson, y = simpson)) + 
  geom_bar(stat = "identity", colour = "black", aes(fill = "grey")) + 
  labs(x = "Plot", y = "Simpson diversity index") + 
  facet_grid(. ~ group, scales = "free_x", space = "free") + 
  theme_bw() + 
  theme(axis.text.x=element_text(
    angle=45, 
    vjust=1, 
    hjust=1)) + 
  theme(legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")) + 
  scale_fill_manual(name = "DBH bin (cm)", values = "grey",
    guide = guide_legend(override.aes = list(colour = "white", fill = "white"))) + 
  scale_y_continuous(limits = c(0,1))

pdf(file = "img/shannon_plot.pdf", width = 15, height = 6)
shannon_by_plot
dev.off()

pdf(file = "img/simpson_plot.pdf", width = 15, height = 6)
simpson_by_plot
dev.off()

shannon_group_lm <- lm(shannon ~ group, data = alpha_df)
summary(shannon_group_lm)
TukeyHSD(aov(shannon_group_lm))

shannon_summ <- alpha_df %>%
  group_by(group) %>%
  summarise(shannon_mean = mean(shannon),
    shannon_sd = sd(shannon))


# Calculate gamma diversity for each group of plots
ab_mat_pool <- lapply(trees_ab_mat_list, function(x){
  colSums(x[,-1])})

gamma <- lapply(ab_mat_pool, function(x){
  diversity(x, index = "shannon")
})

# Proportional species turnover across all sites
bmt_1 <- lapply(1:5, function(x){
  (gamma[[x]] - mean(alpha_list[[x]]$shannon)) / gamma[[x]]
})
names(bmt_1) <- names(ab_mat_pool)

# Pairwise Jaccard-Sorensen beta diversity ----
# Generate pairwise list for each group of plots
species_group_list_plot <- lapply(
  list(bicuar = trees_ab_mat_list$bicuar, drc = trees_ab_mat_list$drc, kilwa = trees_ab_mat_list$kilwa, nham = trees_ab_mat_list$nham), function(x){
    apply(x[,-1], 1, function(y){
      names(y[y > 0])
    })
  })

sample_mat <- combn(species_group_list_plot, 2, simplify = FALSE)

sample_pairs <- lapply(sample_mat, function(x){
  expand.grid(x[[1]], x[[2]])
})

sample_pairs_clean <- lapply(sample_pairs, function(x){
  lapply(1:length(x$Var1), function(y){
    list(x$Var1[[y]], x$Var2[[y]])
  })
})

# Find number of species shared (a) in each pairwise list
shared_n_sp <- lapply(sample_pairs_clean, function(x){
  lapply(x, function(y){
    length(intersect(y[[1]], y[[2]]))
  })
})

# Find number of species unique to community 1
comm_1_n_sp <- lapply(sample_pairs_clean, function(x){
  lapply(x, function(y){
    length(y[[1]][!y[[1]] %in% y[[2]]])
  })
})

# Find number of species unique to community 2
comm_2_n_sp <- lapply(sample_pairs_clean, function(x){
  lapply(x, function(y){
    length(y[[2]][!y[[2]] %in% y[[1]]])
  })
})

# Calculate Jaccard-Sorensen for each pair
js <- lapply(1:6, function(x){
  sapply(1:length(shared_n_sp[[x]]), function(y){
    (2*shared_n_sp[[x]][[y]]) / ((2*shared_n_sp[[x]][[y]]) + comm_1_n_sp[[x]][[y]] + comm_2_n_sp[[x]][[y]]) * 100
  })
})

# Get mean Jaccard-Sorensen for each pairwise group
mean_js <- sapply(js, mean)

pair_names <- t(sapply(sample_mat, names))
pair_names <- data.frame(pair_names)
pair_names$mean_js <- mean_js
names(pair_names) <- c("col1", "col2", "mean_js")

# Create plot
js_pairs <- ggplot(pair_names, aes(x = col1, y = col2, label = round(mean_js, digits = 1))) + 
  geom_point(aes(fill = mean_js), colour = "black", shape = 21, size = 15) +
  geom_text(colour = "white") + 
  scale_fill_gradient(name = expression(S[S]), low = "blue", high = "red") + 
  scale_y_discrete(limits = rev(levels(as.factor(pair_names$col2))), labels = c("Mozambique", "Tanzania", "DRC")) + 
  scale_x_discrete(labels = c("Angola", "DRC", "Tanzania")) + 
  theme_classic() + 
  theme(panel.grid.major = element_line(colour = "grey")) + 
  labs(x = "", y = "")

# Save plot
pdf("img/js_pairs.pdf", width = 8, height = 6)
js_pairs
dev.off()

# Express as a sparse matrix
pairs_js_mat <- pivot_wider(pair_names, names_from = col1, values_from = mean_js)
pairs_js_mat$col2 <- as.character(pairs_js_mat$col2)
pairs_js_mat$bicuar <- round(pairs_js_mat$bicuar, digits = 1)
pairs_js_mat$drc <- round(pairs_js_mat$drc, digits = 1)
pairs_js_mat$kilwa <- round(pairs_js_mat$kilwa, digits = 1)

fileConn <- file(paste0("include/pairs_js.tex"))
writeLines(
  stargazer(pairs_js_mat, summary = FALSE, label = "pairs_js", 
    digit.separate = 0, rownames = FALSE),
  fileConn)
close(fileConn)

# Clean data for structure measures ----

stems_list_fil <- lapply(stems_list, function(x){
  x[,c("plotcode", "species_binomial", "dbh_cm", "agb_stem_eq4")]
})

stems_all <- rbind(stems_list_fil$bicuar, stems_list_fil$drc, 
  stems_list_fil$kilwa, stems_list_fil$nham)

stems_all$group <- case_when(
  grepl("ABG", stems_all$plotcode) ~ "Angola",
  grepl("DKS", stems_all$plotcode) ~ "DRC",
  grepl("TKW", stems_all$plotcode) ~ "Tanzania",
  grepl("MGR", stems_all$plotcode) ~ "Mozambique"
)

# Adjust plotcodes and group names
##' ISO3 country codes
stems_all$plotcode  <- gsub("ABGD-", "ANGD-", stems_all$plotcode) %>%
        gsub("ABG-", "ANG-", .) %>%
        gsub("DKS-", "DRC-", .) %>%
        gsub("TKW-", "TZA-", .) %>%
        gsub("MGR-", "MOZ-", .) %>%
        gsub("-0*", "-", .) %>%
        gsub("DRC-2", "DRC-11", .) %>%
        gsub("DRC-3", "DRC-12", .) %>%
        gsub("DRC-1_", "DRC-", .)

# Calculate variation in basal area across plots by group ----
ba <- stems_all %>%
  mutate(basal_area = ((dbh_cm / 2)^2) * 0.0001) %>%
  group_by(plotcode) %>%
  summarise(basal_area = sum(basal_area, na.rm = TRUE),
    group = first(group))

ba_summ <- ba %>%
  group_by(group, plotcode) %>%
  summarise(total_ba = sum(basal_area, na.rm = TRUE)) %>%
  summarise(mean_ba = mean(total_ba, na.rm = TRUE),
    sd_ba = sd(total_ba, na.rm = TRUE)) %>%
  mutate(cov_ba = sd_ba / mean_ba * 100)

ba_plot <- ggplot() + 
  geom_bar(data = ba_summ,
    aes(x = group, y = mean_ba),
    stat = "identity", colour = "black", fill = "grey") + 
  geom_errorbar(data = ba_summ, 
    aes(x = group, ymax = mean_ba + sd_ba, ymin = mean_ba - sd_ba),
    width = 0.4) +
  theme_classic() + 
  labs(x = "Site", y = expression("Mean plot basal area ("*m^2*")"))

pdf(file = "img/basal_area_group.pdf", height = 6, width = 8)
ba_plot
dev.off()

ba_group_lm <- lm(basal_area ~ group, data = ba)
summary(ba_group_lm)
TukeyHSD(aov(ba_group_lm))

ba$plotcode_ba = factor(ba$plotcode, levels = ba$plotcode[order(desc(ba$basal_area))])

ba_by_plot <- ggplot(ba, aes(x = plotcode_ba, y = basal_area)) + 
  geom_bar(stat = "identity", colour = "black", aes(fill = "grey")) + 
  labs(x = "Site", y = expression("Basal area ("*m^2~ha^-1*")")) + 
facet_grid(. ~ group, scales = "free_x", space = "free") + 
  theme_bw() + 
  theme(axis.text.x=element_text(
    angle=45, 
    vjust=1, 
    hjust=1)) + 
  theme(legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")) + 
  scale_fill_manual(name = "DBH bin (cm)", values = "grey",
    guide = guide_legend(override.aes = list(colour = "white", fill = "white")))

pdf(file = "img/basal_area.pdf", width = 15, height = 6)
ba_by_plot
dev.off()

# Stem size class bins by DBH ----

stems_all <- stems_all %>%
  mutate(dbh_bin = factor(case_when(dbh_cm <= 10.0 ~ "5-10",
    dbh_cm >= 10.0 & dbh_cm <= 20.0 ~ "10-20",
    dbh_cm >= 20.0 & dbh_cm <= 30.0 ~ "20-30",
    dbh_cm >= 30.0 & dbh_cm <= 40.0 ~ "30-40",
    dbh_cm >= 40.0 ~ "40+"), 
    levels = c("5-10", "10-20", "20-30", "30-40", "40+")))

abund_dbh <- stems_all %>%
  group_by(group, plotcode, dbh_bin) %>%
  tally() %>%
  filter(!is.na(dbh_bin)) %>%
  ungroup()

abund_dbh_order_group <- abund_dbh %>%
  group_by(plotcode) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total)) %>%
  mutate(ord = seq(1:length(.$total))) %>%
  dplyr::select(plotcode, ord) %>%
  ungroup()

abund_dbh_ord <- abund_dbh %>%
  left_join(., abund_dbh_order_group, by = "plotcode") %>%
  mutate(plotcode_ord = factor(plotcode, levels = unique(plotcode[order(.$ord)])))

stem_ab_dbh_bin_plot <- ggplot(abund_dbh_ord, 
  aes(x = plotcode_ord, y = n, fill = dbh_bin)) +
  geom_bar(stat = "identity", colour = "black", 
    position = position_stack(reverse = TRUE)) + 
  labs(x = "Plot", y = "Number of stems") + 
  guides(fill=guide_legend(title="DBH bin (cm)")) + 
  scale_fill_brewer(type = "seq", palette = "Greens") + 
  facet_grid(. ~ group, scales = "free_x", space = "free") + 
  theme_bw() + 
  theme(axis.text.x=element_text(
    angle=45, 
    vjust=1, 
    hjust=1)) + 
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500))

pdf(file = "img/stem_ab_dbh_bin.pdf", width = 15, height = 6)
stem_ab_dbh_bin_plot
dev.off()

# Proportion of stems in largest classes
abund_prop <- abund_dbh_ord %>%
  group_by(plotcode, group) %>%
  summarise(n = sum(n))

abund_big <- abund_dbh_ord %>%
  group_by(plotcode, group) %>%
  filter(dbh_bin %in% c("30-40", "40+")) %>%
  summarise(n_big = sum(n)) 

abund_prop_big <- left_join(abund_prop, abund_big, by = "plotcode") %>%
  dplyr::select(plotcode, n, n_big, group = group.x) %>%
  mutate(big_prop = n_big / n * 100)

prop_big_group_lm <- lm(big_prop ~ group, data = abund_prop_big)
summary(prop_big_group_lm)
TukeyHSD(aov(prop_big_group_lm))


# Plot biomass ----

agb_summ <- stems_all %>%
  group_by(group, plotcode) %>%
  summarise(agb_total = sum(agb_stem_eq4, na.rm = TRUE)) %>%
  mutate(plotcode = factor(plotcode, levels = .$plotcode[order(desc(.$agb_total))]))

agb_by_plot <- ggplot(agb_summ, aes(x = plotcode, y = agb_total)) + 
  geom_bar(stat = "identity", colour = "black", aes(fill = "grey")) + 
  labs(x = "Plot", y = expression(paste("Aboveground biomass (Mg ha"^-1,")"))) + 
  facet_grid(. ~ group, scales = "free_x", space = "free") + 
  theme_bw() + 
  theme(axis.text.x=element_text(
    angle=45, 
    vjust=1, 
    hjust=1)) + 
  theme(legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")) + 
  scale_fill_manual(name = "DBH bin (cm)", values = "grey",
    guide = guide_legend(override.aes = list(colour = "white", fill = "white")))

pdf(file = paste0("img/agb_plot.pdf"), width = 15, height = 6)
agb_by_plot
dev.off()

agb_group_summ <- agb_summ %>%
  group_by(group) %>%
  summarise(agb_mean = mean(agb_total, na.rm = TRUE),
    agb_sd = sd(agb_total, na.rm = TRUE))

agb_group <- ggplot(agb_group_summ, aes(x = group, y = agb_mean)) + 
  geom_bar(stat = "identity", colour = "black", fill = "grey") + 
  geom_errorbar(aes(ymin = agb_mean - agb_sd, ymax = agb_mean + agb_sd),
    width = 0.4) + 
  labs(x = "Site", y = expression(paste("Aboveground biomass (Mg ha"^-1,")"))) + 
  theme_classic()

pdf(file = paste0("img/agb_group.pdf"), width = 8, height = 6)
agb_group
dev.off()

# Table describing each group of plots ----
group_descrip <- plot_clim %>%
  group_by(group) %>%
  summarise(mat = round(mean(mat), digits = 1),
    map = round(mean(map), digits = 0),
    #dec_latitude = round(mean(dec_latitude), digits = 2),
    #dec_longitude = round(mean(dec_longitude), digits = 2),
    n_plots = n())

group_descrip_big <- filter(group_descrip, group != "bicuar_degrad")

sp_big <- sapply(trees_ab_mat_list, function(x){
  ncol(x)
})

group_descrip_big$sp <- c(sp_big[[1]], sp_big[[3]], sp_big[[4]], sp_big[[5]])

group_descrip_big$shannon <- paste0(round(shannon_summ$shannon_mean, digits = 2), 
  "(", round(shannon_summ$shannon_sd , digits = 3), ")")

group_descrip_big$ba <- paste0(round(ba_summ$mean_ba, digits = 2), 
  "(", round(ba_summ$sd_ba , digits = 3), ")")

fileConn <- file(paste0("include/group_descrip.tex"))
writeLines(
  stargazer(group_descrip_big, summary = FALSE, label = "group_descrip", 
    digit.separate = 0, rownames = FALSE),
  fileConn)
close(fileConn)

# Map of Bicuar plot locations ----

bicuar_df <- fortify(bicuar_poly)

plot_bicuar_df <- plot_clim %>%
  filter(grepl("ABG", plotcode)) %>%
  dplyr::select(plotcode, group, dec_longitude, dec_latitude) %>%
  mutate(plotcode = gsub("ABG", "", .$plotcode) %>%
      gsub("-0*", "", .))

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
  labs(x = "Longitude", y = "Latitude")

bic_map <- openmap(upperLeft = c(max(bicuar_df$lat) + 0.2, min(bicuar_df$long) - 0.2),
  lowerRight = c(min(bicuar_df$lat) - 0.2, max(bicuar_df$long) + 0.2),
  type = "bing",
  zoom = 10) %>% openproj(., 
    projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
