# Degrad. vs. non-degrad. plots 
# John Godlee (johngodlee@gmail.com)
# 2019_12_20

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(betapart)
library(vegan)
library(ggrepel)

source("scripts/functions.R")

# Import data ----
trees_ab_mat_list <- readRDS("data/trees_ab_mat.rds")
stems_list <- readRDS("data/stems_list.rds")
trees_list <- readRDS("data/trees_list.rds")

# Variation in species composition between Bicuar degrad. plots and big plots ----
# From each of the 15 big plots, randomly sample a 20x50 m subplot
bicuar_trees_list <- split(trees_list$bicuar, trees_list$bicuar$plotcode)

bicuar_trees_sub_list <- lapply(bicuar_trees_list, function(x){
  rand <- sample(1:10, 1)
  
  return(x %>% filter(!is.na(subplot_20_50), subplot_20_50 == rand)) 
})

bicuar_trees_sub_df <- do.call(rbind, bicuar_trees_sub_list) 

# Do the same with stem level data 
bicuar_stems_list <- split(stems_list$bicuar, stems_list$bicuar$plotcode)

bicuar_stems_sub_list <- lapply(bicuar_stems_list, function(x){
	rand <- sample(1:10, 1)
	
	return(x %>% filter(!is.na(subplot_20_50), subplot_20_50 == rand)) 
})

bicuar_stems_sub_df <- do.call(rbind, bicuar_stems_sub_list) 


# Join degrad plots with subsampled big plots
bicuar_trees_all_df <- bicuar_trees_sub_df %>%
  dplyr::select(plotcode, species_binomial, dbh_cm) %>%
  bind_rows(., dplyr::select(trees_list[[2]], plotcode, species_binomial, dbh_cm)) %>%
  mutate(plotcode = gsub("ABG", "", .$plotcode) %>%
      gsub("-0*", "", .),
    group = case_when(grepl("D", plotcode) ~ "degrad",
      TRUE ~ "intact")
  )

degrad_stems <- dplyr::select(stems_list[[2]], plotcode, species_binomial, dbh_cm, base_stem_id) %>%
	mutate(base_stem_id = paste0("D", base_stem_id))

bicuar_stems_all_df <- bicuar_stems_sub_df %>%
	dplyr::select(plotcode, species_binomial, dbh_cm, base_stem_id) %>%
	bind_rows(., degrad_stems) %>%
	mutate(plotcode = gsub("ABG", "", .$plotcode) %>%
			gsub("-0*", "", .),
		group = case_when(grepl("D", plotcode) ~ "degrad",
			TRUE ~ "intact")
	)

# Tree abundance matrix for Bicuar degradation plots and big plots
tree_ab_mat_bicuar_degrad <- bicuar_trees_all_df %>%
  filter(!is.na(species_binomial)) %>%
  group_by(plotcode, species_binomial, .drop = FALSE) %>%
  tally() %>%
  spread(species_binomial, n, fill = 0) %>%
  ungroup()

plot_vec_bicuar_degrad <- tree_ab_mat_bicuar_degrad$plotcode
plot_group_degrad <- data.frame(plotcode = plot_vec_bicuar_degrad, 
  group = case_when(
    grepl("D", plot_vec_bicuar_degrad) ~ "degrad",
    TRUE ~ "bicuar"))

tree_ab_mat_bicuar_degrad_clean <- tree_ab_mat_bicuar_degrad %>% 
  dplyr::select(-plotcode)

# Run NMDS 
NMDS.scree(tree_ab_mat_bicuar_degrad_clean)

bicuar_degrad_tree_ab_nmds <- metaMDS(tree_ab_mat_bicuar_degrad_clean, 
  distance = "bray", try = 500, k = 5)

stressplot(bicuar_degrad_tree_ab_nmds)

# How many tries (Restarts)?
bicuar_degrad_tree_ab_nmds$tries

# Final stress value?
bicuar_degrad_tree_ab_nmds$stress

# Number of axes?
bicuar_degrad_tree_ab_nmds$ndim

# Extract site (plot) scores from NMDS analysis
bicuar_degrad_plot_scores <- data.frame(scores(bicuar_degrad_tree_ab_nmds))  
bicuar_degrad_plot_scores$plot <- plot_vec_bicuar_degrad
bicuar_degrad_plot_scores$group <- case_when(
  grepl("D", bicuar_degrad_plot_scores$plot) ~ "degrad",
  TRUE ~ "inside")

# Extract species scores from NMDS analysis
bicuar_degrad_species_scores <- as.data.frame(scores(bicuar_degrad_tree_ab_nmds, "species")) 
bicuar_degrad_species_scores$species_binomial <- rownames(bicuar_degrad_species_scores)

bicuar_degrad_species_scores_text <- bicuar_degrad_species_scores %>%
	filter(species_binomial %in% c(
		"Brachystegia tamarindoides",
		"Burkea africana", 
		"Baikiaea plurijuga",
		"Baphia massaiensis",
		"Julbernardia paniculata"))

# Plot extracted scores in ggplot2
bicuar_degrad_nmds_plot <- ggplot() + 
  geom_point(data = bicuar_degrad_plot_scores,
    aes(x = NMDS1, y = NMDS2, fill = group), shape = 23, size = 3) +
  geom_point(data = bicuar_degrad_species_scores_text,
    aes(x = NMDS1, y = NMDS2)) + 
  geom_label_repel(data = bicuar_degrad_species_scores_text,
    aes(x = NMDS1, y = NMDS2, label = species_binomial),
    label.padding = 0.08, box.padding = 1, point.padding = 0.15, show.legend = FALSE) +
  coord_equal() +
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_fill_manual(name = "", labels = c("Previously\ndisturbed", "Not disturbed"), values = degrad_pal) +
  scale_colour_manual(name = "", labels = c("Previously\ndisturbed", "Not disturbed"), values = degrad_pal)

pdf(file = paste0("img/bicuar_degrad_nmds.pdf"), width = 5, height = 5)
bicuar_degrad_nmds_plot
dev.off()

# Calculate Beta diversity homoegeneity of variances
dist_sor_abund_degrad <- beta.pair.abund(tree_ab_mat_bicuar_degrad_clean, index.family="bray")
dist_sor_abund_degrad <- betadisper(dist_sor_abund_degrad[[3]], plot_group_degrad$group)
bd_sorr_abund_degrad_plot <- plot(dist_sor_abund_degrad)

bd_sorr_abund_degrad_plot[[1]] <- as.data.frame(bd_sorr_abund_degrad_plot[[1]])
bd_sorr_abund_degrad_plot[[1]]$group <- plot_group_degrad$group
bd_sorr_abund_degrad_plot[[1]]$plot <- plot_group_degrad$plotcode
bd_sorr_abund_degrad_plot[[2]] <- as.data.frame(bd_sorr_abund_degrad_plot[[2]])
bd_sorr_abund_degrad_plot[[2]]$group <- row.names(bd_sorr_abund_degrad_plot[[2]]) 
bd_sorr_abund_degrad_plot[[3]] <- left_join(bd_sorr_abund_degrad_plot[[1]], bd_sorr_abund_degrad_plot[[2]], by = "group")
bd_sorr_abund_degrad_plot[[3]]$line <- row.names(bd_sorr_abund_degrad_plot[[3]])
bd_a <- bd_sorr_abund_degrad_plot[[3]][,c(1:3,7)]
bd_b <- bd_sorr_abund_degrad_plot[[3]][,c(5,6,3,7)]
names(bd_b) <- names(bd_a)
lines_df <- rbind(bd_a, bd_b)
names(lines_df) <- c("PCoA1", "PCoA2", "group", "line")

hulls_list <- split(bd_sorr_abund_degrad_plot[[1]], bd_sorr_abund_degrad_plot[[1]]$group)
hulls_list <- lapply(hulls_list, find_hull)
hulls <- do.call(rbind, hulls_list)

# Create plot in ggplot2
bicuar_trees_degrad_beta_div_plot <- ggplot() + 
  geom_line(data = lines_df, 
    aes(x = PCoA1, y = PCoA2, colour = group, group = line), 
    show.legend = FALSE) + 
  geom_polygon(data = hulls,
    aes(x = PCoA1, y = PCoA2, fill = group), 
    alpha = 0.5) + 
  geom_point(data = bd_sorr_abund_degrad_plot[[1]], 
    aes(x = PCoA1, y = PCoA2, fill = group),
    shape = 21, size = 2) + 
  geom_point(data = bd_sorr_abund_degrad_plot[[2]], 
    aes(x = PCoA1, y = PCoA2, fill = group),
    shape = 23, size = 4) + 
  # geom_label(data = bd_sorr_abund_degrad_plot[[1]],
  #   aes(x = PCoA1, y = PCoA2, fill = group, label = plot),
  #   size = 4, show.legend = FALSE) +
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_fill_manual(name = "", values = degrad_pal, labels = c("Degraded", "Non-degraded")) + 
  scale_colour_manual(name = "", values = degrad_pal, labels = c("Degraded", "Non-degraded"))

pdf(file = paste0("img/bicuar_degrad_beta_div.pdf"), width = 8, height = 6)
bicuar_trees_degrad_beta_div_plot
dev.off()

# Which species most common in degraded plots? ----
degrad_species <- trees_list[[2]] %>%
  group_by(species_binomial) %>%
  tally() %>%
  arrange(desc(n))

big_species <- bicuar_trees_sub_df %>%
  group_by(species_binomial) %>%
  tally() %>%
  arrange(desc(n))

# Which species were only found in degraded plots? ----
degrad_only_species <- setdiff(degrad_species$species_binomial, big_species$species_binomial)
big_only_species <- setdiff(big_species$species_binomial, degrad_species$species_binomial)

degrad_species %>% 
  filter(species_binomial %in% degrad_only_species)

bigonlyspeciesn <- big_species %>%
  filter(species_binomial %in% big_only_species)

# What are the mean dimensions of the most common species? ----
common_dim_degrad <- stems_list$bicuar_degrad %>%
  filter(species_binomial == "Baphia massaiensis") %>%
  summarise(mean_dbh = mean(dbh_cm, na.rm = TRUE),
    mean_height = mean(height_m, na.rm = TRUE),
    sd_dbh = sd(dbh_cm, na.rm = TRUE),
    sd_height = sd(height_m, na.rm = TRUE))

common_dim_bicuar <- stems_list$bicuar %>%
  filter(species_binomial == "Julbernardia paniculata") %>%
  summarise(mean_dbh = mean(dbh_cm, na.rm = TRUE),
    mean_height = mean(height_m, na.rm = TRUE),
    sd_dbh = sd(dbh_cm, na.rm = TRUE),
    sd_height = sd(height_m, na.rm = TRUE))

# What is the mean stem density? ----
stem_dens_bicuar <- stems_list$bicuar %>%
	group_by(plotcode) %>%
	tally()

stem_dens_degrad <- stems_list$bicuar_degrad %>%
	group_by(plotcode) %>%
	tally()

# Which plots have more multi stemmed trees? ----
multi_stem <- bicuar_stems_all_df %>%
	group_by(plotcode, base_stem_id) %>%
	tally() %>%
	filter(n > 1) %>%
	mutate(group = case_when(
		grepl("D", plotcode) ~ "degrad",
		TRUE ~ "intact")) %>%
	ungroup() %>%
	group_by(group) %>%
	summarise(mean_n = mean(n, na.rm = TRUE),
		sd_n = sd(n, na.rm = TRUE))

bicuar_stems_all_df_bin <- bicuar_stems_all_df %>%
	mutate(dbh_bin = factor(case_when(dbh_cm <= 10.0 ~ "5-10",
		dbh_cm >= 10.0 & dbh_cm <= 20.0 ~ "10-20",
		dbh_cm >= 20.0 & dbh_cm <= 30.0 ~ "20-30",
		dbh_cm >= 30.0 & dbh_cm <= 40.0 ~ "30-40",
		dbh_cm >= 40.0 ~ "40+"), 
		levels = c("5-10", "10-20", "20-30", "30-40", "40+")))

bicuar_stems_all_df_bin_list <- split(bicuar_stems_all_df_bin, bicuar_stems_all_df_bin$dbh_bin)

lapply(bicuar_stems_all_df_bin_list, function(x){
	df <- x %>%
		group_by(plotcode, group) %>%
		tally()
	
	mod <- glm(n ~ group, data = df, family = "poisson")
	
	summary(mod)
})

bicuar_stems_all_df_bin_mean <- bicuar_stems_all_df_bin %>%
	group_by(group, plotcode, dbh_bin) %>%
	tally() %>%
	filter(!is.na(dbh_bin)) %>% 
	group_by(group, dbh_bin) %>%
	summarise(mean_n = mean(n, na.rm = TRUE),
		samples = n(),
		sem = sd(n, na.rm = TRUE) / sqrt(samples)) %>%
	ungroup() %>%
	mutate(group = factor(group,
		levels = c("degrad", "intact"), labels = c("Disturbed", "Not disturbed")))

degrad_dbh_bin <- ggplot(bicuar_stems_all_df_bin_mean,
	aes(x = dbh_bin, y = mean_n)) + 
	geom_bar(aes(fill = group), 
		stat = "identity", colour = "black", position = position_dodge(width=0.9)) +
	geom_errorbar(aes(ymin = mean_n - sem, ymax = mean_n + sem, group = group), 
		width = 0.5, position = position_dodge(width=0.9)) + 
	theme_classic() + 
	theme(axis.text.x=element_text(
		angle=45, 
		vjust=1, 
		hjust=1)) + 
	labs(x = "Stem diameter bin (cm)", y = "Number of stems") +
	theme(legend.position = c(0.9, 0.5)) + 
	scale_fill_manual(name = "", values = degrad_pal, labels = c("Disturbed", "Not disturbed")) + 
	geom_text(x = 1, y = 85, label = "***", size = 8) + 
	geom_text(x = 2, y = 25, label = "***", size = 8)
	

pdf("img/degrad_dbh_bin.pdf", width = 6, height = 4)
degrad_dbh_bin
dev.off()

# Estimate shannon equitability for evenness in degrad and non-degrad. ----

plot_group_degrad$shannon <- diversity(tree_ab_mat_bicuar_degrad_clean)

plot_group_degrad$rich <- rowSums(tree_ab_mat_bicuar_degrad_clean != 0)

plot_group_degrad$equit <- plot_group_degrad$shannon / log(plot_group_degrad$rich)

plot_group_degrad_summ <- plot_group_degrad %>%
  group_by(group) %>%
  summarise(shannon_mean = mean(shannon),
    shannon_sd = sd(shannon),
    rich_mean = mean(rich),
    rich_sd = sd(rich, na.rm = TRUE),
    equit_mean = mean(equit),
    equit_sd = sd(equit))

bicuar_div_gather <- bicuar_trees_all_df %>%
  mutate(ba = ((dbh_cm / 2)^2) * 0.0001) %>%
  group_by(plotcode) %>%
  summarise(ba = sum(ba)) %>%
  left_join(.,plot_group_degrad, by = "plotcode") %>%
  gather("index", "value", -group, -plotcode) %>%
  mutate(index_plot = case_when(
    index == "ba" ~ "Basal area",
    index == "equit" ~ "Shannon equitability",
    index == "rich" ~ "Species richness",
    index == "shannon" ~ "Shannon-Wiener index"
  )) %>%
  mutate(index_plot = factor(index_plot, 
    levels = c("Species richness", "Basal area", 
      "Shannon-Wiener index", "Shannon equitability"),
    labels = c(expression("Species"~"richness"), 
      expression("Basal"~"area"~(m^2~ha^-1)), 
      expression("Shannon-Wiener"~"index"), 
      expression("Shannon"~"equitability")))) %>%
  mutate(group = case_when(
    group == "bicuar" ~ "Not disturbed",
    group == "degrad" ~ "Disturbed"
  ))

sp_degrad_lm <- lm(rich ~ group, data = plot_group_degrad)
summary(sp_degrad_lm)
TukeyHSD(aov(sp_degrad_lm))

ba_degrad_lm <- bicuar_div_gather %>% filter(index == "ba") %>%
  mutate(group = case_when(
    group == "Not disturbed" ~ "bicuar",
    group == "Disturbed" ~ "degrad"
  )) %>%
  lm(value ~ group, data = .)

summary(ba_degrad_lm)
TukeyHSD(aov(ba_degrad_lm))

sh_degrad_lm <- lm(shannon ~ group, data = plot_group_degrad)
summary(sh_degrad_lm)
TukeyHSD(aov(sh_degrad_lm))

se_degrad_lm <- lm(equit ~ group, data = plot_group_degrad)
summary(se_degrad_lm)
TukeyHSD(aov(se_degrad_lm))

fileConn <- file(paste0("include/degrad_anova_table.tex"))
writeLines(
  stargazer(sp_degrad_lm, ba_degrad_lm, sh_degrad_lm, se_degrad_lm, 
    notes.label	= ""),
  fileConn)
close(fileConn)

# Boxplots ----

box_degrad <- ggplot() + 
  geom_boxplot(data = bicuar_div_gather,
    aes(x = group, y = value, fill = group)) + 
  facet_wrap(~index_plot, scales = "free_y", labeller = label_parsed) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_fill_manual(name = "", values = degrad_pal) + 
  labs(x = "", y = "")

pdf("img/degrad_box.pdf", width = 6, height = 6)
box_degrad
dev.off()


  


# Output numerical figures ----
bmdbhdegrad <- paste0(round(common_dim_degrad$mean_dbh, digits = 1), 
  "$\\pm$",
  round(common_dim_degrad$sd_dbh, digits = 2))
bmheightdegrad <-  paste0(round(common_dim_degrad$mean_height, digits = 1), 
  "$\\pm$",
  round(common_dim_degrad$sd_height, digits = 2))
jpdbhbicuar <-  paste0(round(common_dim_bicuar$mean_dbh, digits = 1), 
  "$\\pm$",
  round(common_dim_bicuar$sd_dbh, digits = 2))
jpheightbicuar <-  paste0(round(common_dim_bicuar$mean_height, digits = 1), 
  "$\\pm$",
  round(common_dim_bicuar$sd_height, digits = 2))

ndegradplots <- length(unique(trees_list[[2]]$plotcode))
nbmdegrad <- c(degrad_species[degrad_species$species_binomial == "Baphia massaiensis", "n"])
njpdegrad <- c(big_species[big_species$species_binomial == "Julbernardia paniculata", "n"])
 
degradshannon <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "shannon" & bicuar_div_gather$group == "Disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "shannon" & bicuar_div_gather$group == "Disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "shannon" & bicuar_div_gather$group == "Disturbed", "value"]$value)), 2))

bicuarsubshannon <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "shannon" & bicuar_div_gather$group == "Not disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "shannon" & bicuar_div_gather$group == "Not disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "shannon" & bicuar_div_gather$group == "Not disturbed", "value"]$value)), 2))

degradequit <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "equit" & bicuar_div_gather$group == "Disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "equit" & bicuar_div_gather$group == "Disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "equit" & bicuar_div_gather$group == "Disturbed", "value"]$value)), 1))

bicuarsubequit <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "equit" & bicuar_div_gather$group == "Not disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "equit" & bicuar_div_gather$group == "Not disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "equit" & bicuar_div_gather$group == "Not disturbed", "value"]$value)), 2))

bicuarsubrich <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "rich" & bicuar_div_gather$group == "Not disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "rich" & bicuar_div_gather$group == "Not disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "rich" & bicuar_div_gather$group == "Not disturbed", "value"]$value)), 2))

degradrich <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "rich" & bicuar_div_gather$group == "Disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "rich" & bicuar_div_gather$group == "Disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "rich" & bicuar_div_gather$group == "Disturbed", "value"]$value)), 2))

bicuarsubba <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "ba" & bicuar_div_gather$group == "Not disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "ba" & bicuar_div_gather$group == "Not disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "ba" & bicuar_div_gather$group == "Not disturbed", "value"]$value)), 2))

degradba <- paste0(
  round(mean(bicuar_div_gather[bicuar_div_gather$index == "ba" & bicuar_div_gather$group == "Disturbed", "value"]$value), 1),
  "$\\pm$",
  round(sd(bicuar_div_gather[bicuar_div_gather$index == "ba" & bicuar_div_gather$group == "Disturbed", "value"]$value) / 
      sqrt(length(bicuar_div_gather[bicuar_div_gather$index == "ba" & bicuar_div_gather$group == "Disturbed", "value"]$value)), 2))

stemdensbicuar <- paste0(
	round(mean(stem_dens_bicuar$n), 1),
	"$\\pm$",
	round(sd(stem_dens_bicuar$n), 2)
)

stemdensdegrad <- paste0(
	round(mean(stem_dens_degrad$n * 10), 1),
	"$\\pm$",
	round(sd(stem_dens_degrad$n * 10), 2)
)

multistemdegrad <- paste0(
	round(multi_stem$mean_n[1], 1),
	"$\\pm$",
	round(multi_stem$sd_n[1], 2))

multistembicuar <- paste0(
	round(multi_stem$mean_n[2], 1),
	"$\\pm$",
	round(multi_stem$sd_n[2], 2))

ndegradonlyspecies <- length(degrad_only_species)
nbigonlyspecies <- length(big_only_species)
nccdegrad <- degrad_species[degrad_species$species_binomial == "Combretum celastroides",]$n
nvrdegrad <- degrad_species[degrad_species$species_binomial == "Acacia reficiens",]$n
ngtdegrad <- degrad_species[degrad_species$species_binomial == "Gardenia ternifolia",]$n

nbsbig <- bigonlyspeciesn[bigonlyspeciesn$species_binomial == "Brachystegia spiciformis",]$n
nbpbig <- bigonlyspeciesn[bigonlyspeciesn$species_binomial == "Baikiaea plurijuga",]$n
ncabig <- bigonlyspeciesn[bigonlyspeciesn$species_binomial == "Combretum apiculatum",]$n


fileConn <- file("include/degrad_figures.tex")
writeLines(
  c(
    paste0("\\newcommand{\\ndegradplots}{", ndegradplots, "}"),
    paste0("\\newcommand{\\nbmdegrad}{", nbmdegrad, "}"),
    paste0("\\newcommand{\\njpdegrad}{", njpdegrad, "}"),
    paste0("\\newcommand{\\degradshannon}{", degradshannon, "}"),
    paste0("\\newcommand{\\bicuarsubshannon}{", bicuarsubshannon, "}"),
    paste0("\\newcommand{\\degradrich}{", degradrich, "}"),
    paste0("\\newcommand{\\bicuarsubrich}{", bicuarsubrich, "}"),
    paste0("\\newcommand{\\degradba}{", degradba, "}"),
    paste0("\\newcommand{\\bicuarsubba}{", bicuarsubba, "}"),
    paste0("\\newcommand{\\degradequit}{", degradequit, "}"),
    paste0("\\newcommand{\\bicuarsubequit}{", bicuarsubequit, "}"),
    paste0("\\newcommand{\\lmshannondegrad}{", lm_format(shannon_degrad_lm), "}"),
    paste0("\\newcommand{\\lmequitdegrad}{", lm_format(equit_degrad_lm), "}"),
    paste0("\\newcommand{\\ndegradonlyspecies}{", ndegradonlyspecies, "}"),
    paste0("\\newcommand{\\nbigonlyspecies}{", nbigonlyspecies, "}"),
    paste0("\\newcommand{\\nccdegrad}{", nccdegrad, "}"),
    paste0("\\newcommand{\\nvrdegrad}{", nvrdegrad, "}"),
    paste0("\\newcommand{\\ngtdegrad}{", ngtdegrad, "}"),
    paste0("\\newcommand{\\nbsbig}{", nbsbig, "}"),
    paste0("\\newcommand{\\nbpbig}{", nbpbig, "}"),    
    paste0("\\newcommand{\\ncabig}{", ncabig, "}"),
    paste0("\\newcommand{\\bmdbhdegrad}{", bmdbhdegrad, "}"),
    paste0("\\newcommand{\\bmheightdegrad}{", bmheightdegrad, "}"),
    paste0("\\newcommand{\\jpdbhbicuar}{", jpdbhbicuar, "}"),
    paste0("\\newcommand{\\jpheightbicuar}{", jpheightbicuar, "}"),
    paste0("\\newcommand{\\stemdensbicuar}{", stemdensbicuar, "}"),
    paste0("\\newcommand{\\stemdensdegrad}{", stemdensdegrad, "}"),
  	paste0("\\newcommand{\\multistemdegrad}{", multistemdegrad, "}"),
  	paste0("\\newcommand{\\multistembicuar}{", multistembicuar, "}")
  ),
  fileConn)
close(fileConn)