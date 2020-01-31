# Alpha and beta diversity estimates 
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
library(vegan)
library(stargazer)

source("scripts/functions.R")

# Import data 
trees_ab_mat_list <- readRDS("data/trees_ab_mat.rds")
trees_ab_mat_list_big <- trees_ab_mat_list[c("bicuar", "drc", "kilwa", "nham")]
stems_list <- readRDS("data/stems_list.rds")

# Calculate shannon, evenness richness for each plot ----
div_df <- function(ab_mat){
  ab_mat_clean <- ab_mat[,-1]
  plotcode <- ab_mat[,1]
  
  alpha <- data.frame(
    plotcode, 
    shannon = diversity(ab_mat_clean), 
    rich = rowSums(ab_mat_clean != 0))
  alpha$shannon_equit <- alpha$shannon / log(alpha$rich)

  return(alpha)
}

alpha_list <- lapply(trees_ab_mat_list, div_df)
saveRDS(alpha_list, "data/plot_div_list.rds")

alpha_list_big <- alpha_list[c(1,3:5)]

alpha_df <- do.call(rbind, alpha_list_big)
alpha_df$group <- gsub("\\..*$", "", rownames(alpha_df))

shannon_group_lm <- lm(shannon ~ group, data = alpha_df)
summary(shannon_group_lm)
shannon_group_tukey <- TukeyHSD(aov(shannon_group_lm))

# Basal area
stems_list_bin <- lapply(stems_list, function(x){
  x %>%
    mutate(dbh_bin = factor(case_when(dbh_cm <= 10.0 ~ "5-10",
      dbh_cm >= 10.0 & dbh_cm <= 20.0 ~ "10-20",
      dbh_cm >= 20.0 & dbh_cm <= 30.0 ~ "20-30",
      dbh_cm >= 30.0 & dbh_cm <= 40.0 ~ "30-40",
      dbh_cm >= 40.0 ~ "40+"), 
      levels = c("5-10", "10-20", "20-30", "30-40", "40+")))
})

stems_list_fil <- lapply(stems_list_bin, function(x){
  x[,c("plotcode", "species_binomial", "dbh_cm", "dbh_bin")]
})

stems_list_big <- stems_list_fil[c("bicuar", "drc", "kilwa", "nham")]

stems <- do.call(rbind, stems_list_big)

# Find 5-10 cm basal area proportion in other plots
est_ba <- stems %>%
  mutate(basal_area = ((dbh_cm / 2)^2) * 0.0001) %>%
  group_by(plotcode, dbh_bin) %>%
  summarise(basal_area = sum(basal_area, na.rm = TRUE)) %>%
  mutate(prop_ba = basal_area / sum(basal_area)) %>%
  filter(dbh_bin == "5-10") %>%
  mutate(group =  case_when(
    grepl("ABG", plotcode) ~ "bicuar",
    grepl("TKW", plotcode) ~ "kilwa",
    grepl("MGR", plotcode) ~ "nham"
  )) %>%
  filter(!is.na(group))

# Run a linear model
##' does the proportion of basal area in small stems differ between plot groups?
lm_ba_small_stems <- lm(prop_ba ~ group, data = est_ba)
summary(lm_ba_small_stems)
##' Nope

# Mean basal area of 5-10 cm stems as a proportion of total ba
mean_small_ba <- mean(est_ba$prop_ba)

ba <- stems %>%
  mutate(basal_area = ((dbh_cm / 2)^2) * 0.0001) %>%
  group_by(plotcode) %>%
  summarise(basal_area = sum(basal_area, na.rm = TRUE))

# Create dataframe with estimated ba for DRC
alpha_df <- left_join(alpha_df, ba, by = "plotcode") %>%
  mutate(basal_area = case_when(
    group == "drc" ~ (basal_area + (basal_area * mean_small_ba)),
    TRUE ~ basal_area
  ))

# Linear model, groups vs. basal area
ba_group_lm <- lm(basal_area ~ group, data = alpha_df)
summary(ba_group_lm)
ba_group_tukey <- TukeyHSD(aov(ba_group_lm))

alpha_df$plotcode <-
  gsub("ABG-", "ANG-", alpha_df$plotcode) %>%
  gsub("DKS-", "DRC-", .) %>%
  gsub("TKW-", "TZA-", .) %>%
  gsub("MGR-", "MOZ-", .) %>%
  gsub("-0*", "-", .) %>%
  gsub("DRC-2", "DRC-11", .) %>%
  gsub("DRC-3", "DRC-12", .) %>%
  gsub("DRC-1_", "DRC-", .)

# Three panel boxplot of diversity by group ----
alpha_df_gather <- alpha_df %>%
  group_by(group) %>%
  gather("index", "value", -group, -plotcode) %>%
  mutate(index = case_when(
    index == "rich" ~ "Species richness",
    index == "shannon" ~ "Shannon-Wiener index", 
    index == "shannon_equit" ~ "Shannon equitability",
    index == "basal_area" ~ "Basal area",
    TRUE ~ index)) %>%
  mutate(index = factor(index, 
    levels = c("Species richness", "Basal area", 
      "Shannon-Wiener index", "Shannon equitability"),
    labels = c(expression("Species"~"richness"), 
      expression("Basal"~"area"~(m^2~ha^-1)), 
      expression("Shannon-Wiener"~"index"), 
      expression("Shannon"~"equitability")))) %>%
  ungroup() %>%
  mutate(group = case_when(
    group == "bicuar" ~ "Angola",
    group == "drc" ~ "DRC", 
    group == "kilwa" ~ "Tanzania",
    group == "nham" ~ "Mozambique"))

pdf("img/div_box.pdf", width = 12, height = 6)
ggplot() +
  geom_boxplot(data = alpha_df_gather, aes(x = group, y = value, fill = group)) + 
  facet_wrap(~index, scales = "free_y", labeller = label_parsed) +
  theme_classic() + 
  scale_fill_manual(name = "", values = big_pal) + 
  theme(legend.position = "none") + 
  labs(x = "Site", y = "")
dev.off()

# Pairwise Jaccard-Sorensen beta diversity ----
# Generate pairwise list for each group of plots
species_group_list_plot <- lapply(
  list(bicuar = trees_ab_mat_list$bicuar, drc = trees_ab_mat_list$drc, 
    kilwa = trees_ab_mat_list$kilwa, nham = trees_ab_mat_list$nham), function(x){
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
    (2*shared_n_sp[[x]][[y]]) / ((2*shared_n_sp[[x]][[y]]) + 
      comm_1_n_sp[[x]][[y]] + comm_2_n_sp[[x]][[y]]) * 100
  })
})

# Get mean Jaccard-Sorensen for each pairwise group
mean_js <- sapply(js, mean)

pair_names <- t(sapply(sample_mat, names))
pair_names <- data.frame(pair_names)
pair_names$mean_js <- mean_js
names(pair_names) <- c("col1", "col2", "mean_js")

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

# Site level pairwise Sorensen index ----
species_site_list <- lapply(species_group_list_plot, function(x){
  unique(do.call(c, x))
  })

site_sample_mat <- combn(species_site_list, 2, simplify = FALSE)

# Find number of species shared (a) in each pairwise list
site_shared_n_sp <- lapply(site_sample_mat, function(x){
  length(intersect(x[[1]], x[[2]]))
})

# Find number of species unique to community 1
site_comm_1_n_sp <- lapply(site_sample_mat, function(x){
  length(setdiff(x[[1]], x[[2]]))
})

# Find number of species unique to community 2
site_comm_2_n_sp <- lapply(site_sample_mat, function(x){
  length(setdiff(x[[2]], x[[1]]))
})

# Calculate Jaccard-Sorensen for each pair
site_js <- sapply(1:6, function(x){
  sapply(1:length(site_shared_n_sp[[x]]), function(y){
    (2*site_shared_n_sp[[x]][[y]]) / ((2*site_shared_n_sp[[x]][[y]]) + 
      site_comm_1_n_sp[[x]][[y]] + site_comm_2_n_sp[[x]][[y]]) * 100
  })
})

# Create clean dataframe for printing
site_pair_names <- pair_names[,1:2]
site_pair_names_print <- site_pair_names %>%
  mutate(js = sprintf("%.1f", round(site_js, digits = 1)),
    sp_share = unlist(site_shared_n_sp),
    col1 = paste0(col1, "(", unlist(site_comm_1_n_sp), ")"),
    col2 = paste0(col2, "(", unlist(site_comm_2_n_sp), ")")) 

# Save to table
fileConn <- file(paste0("include/site_pairs_js.tex"))
writeLines(
  stargazer(site_pair_names_print, summary = FALSE, label = "site_pairs_js", 
    digit.separate = 0, rownames = FALSE),
  fileConn)
close(fileConn)

# Numerical outputs ----
bicuarbamin <- alpha_df %>% 
  filter(group == "bicuar") %>%
  filter(basal_area == min(basal_area)) %>%
  pull(basal_area) %>%
  round(., 2)

bicuarbamax <- alpha_df %>% 
  filter(basal_area == max(basal_area)) %>%
  pull(basal_area) %>%
  round(., 2)

ba_mean <- alpha_df %>% 
  group_by(group) %>%
  summarise(mean_ba = mean(basal_area),
    sem_ba = sd(basal_area) / sqrt(n()))

babicuar <- paste0(round(ba_mean[ba_mean$group == "bicuar",]$mean_ba, 2), 
  "$\\pm$", round(ba_mean[ba_mean$group == "bicuar",]$sem_ba, 3))
badrc <- paste0(round(ba_mean[ba_mean$group == "drc",]$mean_ba, 2), 
  "$\\pm$", round(ba_mean[ba_mean$group == "drc",]$sem_ba, 3))
banham <- paste0(round(ba_mean[ba_mean$group == "nham",]$mean_ba, 2), 
    "$\\pm$", round(ba_mean[ba_mean$group == "nham",]$sem_ba, 3))
bakilwa <- paste0(round(ba_mean[ba_mean$group == "kilwa",]$mean_ba, 2), 
  "$\\pm$", round(ba_mean[ba_mean$group == "kilwa",]$sem_ba, 3))

bicuarshannon <- paste0(
  round(mean(alpha_df[alpha_df$group == "bicuar",]$shannon), 1), 
  "$\\pm$",
  round(sd(alpha_df[alpha_df$group == "bicuar",]$shannon) / sqrt(length(alpha_df[alpha_df$group == "bicuar",]$shannon)), 2)
  )
nhamshannon <- paste0(
  round(mean(alpha_df[alpha_df$group == "nham",]$shannon), 1), 
  "$\\pm$",
  round(sd(alpha_df[alpha_df$group == "nham",]$shannon) / sqrt(length(alpha_df[alpha_df$group == "nham",]$shannon)), 2)
  )
kilwashannon <- paste0(
  round(mean(alpha_df[alpha_df$group == "kilwa",]$shannon), 1), 
  "$\\pm$",
  round(sd(alpha_df[alpha_df$group == "kilwa",]$shannon) / sqrt(length(alpha_df[alpha_df$group == "kilwa",]$shannon)), 2)
  )
drcshannon <- paste0(
  round(mean(alpha_df[alpha_df$group == "drc",]$shannon), 1), 
  "$\\pm$",
  round(sd(alpha_df[alpha_df$group == "drc",]$shannon) / sqrt(length(alpha_df[alpha_df$group == "drc",]$shannon)), 2)
  )
bicuarminshannon <- round(min(alpha_list$bicuar$shannon), 2)
bicuarminshannonplot <- alpha_list$bicuar[
  alpha_list$bicuar$shannon == min(alpha_list$bicuar$shannon),]$plotcode
bicuarmaxshannon <- round(max(alpha_list$bicuar$shannon), 2)
bicuarmaxshannonplot <- alpha_list$bicuar[
  alpha_list$bicuar$shannon == max(alpha_list$bicuar$shannon),]$plotcode
threshbabicuar <- ba %>% 
  filter(grepl("ABG-", plotcode)) %>%
  pull(basal_area) %>% 
  max() %>% 
  round(., 1)
tukeyshannonbicuardrc <- p_format(shannon_group_tukey$group[1,4])
tukeyshannonbicuarkilwa <- p_format(shannon_group_tukey$group[2,4])
tukeyshannonbicuarnham <- p_format(shannon_group_tukey$group[3,4])
tukeybabicuardrc <- p_format(ba_group_tukey$group[1,4])
tukeybabicuarkilwa <- p_format(ba_group_tukey$group[2,4])
tukeybabicuarnham <- p_format(ba_group_tukey$group[3,4])
lmshannon <- lm_format(shannon_group_lm)
lmba <- lm_format(ba_group_lm)

fileConn <- file("include/plot_div_figures.tex")
writeLines(
  c(
    paste0("\\newcommand{\\bicuarshannon}{", bicuarshannon, "}"),
    paste0("\\newcommand{\\nhamshannon}{", nhamshannon, "}"),
    paste0("\\newcommand{\\kilwashannon}{", kilwashannon, "}"),
    paste0("\\newcommand{\\drcshannon}{", drcshannon, "}"),
    paste0("\\newcommand{\\bicuarminshannon}{", bicuarminshannon, "}"),
    paste0("\\newcommand{\\bicuarminshannonplot}{", bicuarminshannonplot, "}"),
    paste0("\\newcommand{\\bicuarmaxshannon}{", bicuarmaxshannon, "}"),
    paste0("\\newcommand{\\bicuarmaxshannonplot}{", bicuarmaxshannonplot, "}"),
    paste0("\\newcommand{\\threshbabicuar}{", threshbabicuar, "}"),
    paste0("\\newcommand{\\tukeyshannonbicuardrc}{", tukeyshannonbicuardrc, "}"),
    paste0("\\newcommand{\\tukeyshannonbicuarkilwa}{", tukeyshannonbicuarkilwa, "}"),
    paste0("\\newcommand{\\tukeyshannonbicuarnham}{", tukeyshannonbicuarnham, "}"),
    paste0("\\newcommand{\\tukeybabicuardrc}{", tukeybabicuardrc, "}"),
    paste0("\\newcommand{\\tukeybabicuarkilwa}{", tukeybabicuarkilwa, "}"),
    paste0("\\newcommand{\\tukeybabicuarnham}{", tukeybabicuarnham, "}"),
    paste0("\\newcommand{\\lmshannon}{", lmshannon, "}"),
    paste0("\\newcommand{\\lmba}{", lmba, "}"),
    paste0("\\newcommand{\\babicuar}{", babicuar, "}"),
    paste0("\\newcommand{\\badrc}{", badrc, "}"),
    paste0("\\newcommand{\\banham}{", banham, "}"),
    paste0("\\newcommand{\\bakilwa}{", bakilwa, "}"),
    paste0("\\newcommand{\\bicuarbamin}{", bicuarbamin, "}"),
    paste0("\\newcommand{\\bicuarbamax}{", bicuarbamax, "}"),
    paste0("\\newcommand{\\lmbasmallstems}{", lm_format(lm_ba_small_stems), "}")
    ),
  fileConn)
close(fileConn)

