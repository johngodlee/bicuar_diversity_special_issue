# Plot of stem abundance by DBH bin
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
library(broom)

source("scripts/functions.R")

# Import data ----
stems_list <- readRDS("data/stems_list.rds")

# Create DBH bins in all plots
stems_list_bin <- lapply(stems_list, function(x){
  x %>%
    mutate(dbh_bin = factor(case_when(dbh_cm <= 10.0 ~ "5-10",
      dbh_cm >= 10.0 & dbh_cm <= 20.0 ~ "10-20",
      dbh_cm >= 20.0 & dbh_cm <= 30.0 ~ "20-30",
      dbh_cm >= 30.0 & dbh_cm <= 40.0 ~ "30-40",
      dbh_cm >= 40.0 ~ "40+"), 
      levels = c("5-10", "10-20", "20-30", "30-40", "40+")))
})

# Estimate 5-10 cm stems in DRC ----
# Find 5-10 cm stem proportion in other plots
est_stems <- lapply(stems_list_bin, function(x){
  x %>%
  group_by(plotcode, dbh_bin) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(dbh_bin == "5-10")
})

# Collapse to dataframe
est_stems_other_df <- do.call(rbind, est_stems[c("bicuar", "kilwa", "nham")])

est_stems_other_df$group <- case_when(
  grepl("ABG", est_stems_other_df$plotcode) ~ "bicuar",
  grepl("TKW", est_stems_other_df$plotcode) ~ "kilwa",
  grepl("MGR", est_stems_other_df$plotcode) ~ "nham"
)

# Run a linear model
##' does the proportion of small stems differ between plot groups?
lm_small_stems <- lm(freq ~ group, data = est_stems_other_df)
summary(lm_small_stems)
##' Nope

# Get mean frequency across plot groups
mean_small_freq <- mean(est_stems_other_df$freq)

# Get number of stems in DRC, estimate 5-10 cm stems from regional mean proportion
small_stems <- stems_list$drc %>%
  filter(dbh_cm >= 10) %>%
  group_by(plotcode) %>%
  tally() %>%
  mutate(n_small = n / (1 - mean_small_freq))

# Generate small stems rows
small_stems_rows <- small_stems[rep(row.names(small_stems), small_stems$n_small), 1]
small_stems_rows$dbh_bin <- "5-10"

# Stem abundance by DBH bin ----
# Isolate big plots
stems_big <- lapply(stems_list_bin, function(x){
  dplyr::select(x, plotcode, dbh_bin)
})

# Collapse to dataframe
stems <- do.call(rbind, stems_big[c("bicuar", "drc", "kilwa", "nham")])

# Add small stem estimates from drc
stems_est_df <- rbind(stems, small_stems_rows)

# Add plot group id
stems_est_df$group <- case_when(
  grepl("ABG", stems_est_df$plotcode) ~ "Angola",
  grepl("DKS", stems_est_df$plotcode) ~ "DRC",
  grepl("TKW", stems_est_df$plotcode) ~ "Tanzania",
  grepl("MGR", stems_est_df$plotcode) ~ "Mozambique",
  TRUE ~ "seosaw")

# Summarise dataframe, add column for estimated
abund_dbh <- stems_est_df %>%
  group_by(group, plotcode, dbh_bin) %>%
  tally() %>%
  filter(!is.na(dbh_bin)) %>%
  ungroup() %>%
  mutate(est = case_when(
    group == "DRC" & dbh_bin == "5-10" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(est = factor(est, levels = c(TRUE, FALSE)))

# Order bars within a facet in descending total stem freq.
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

# Create plot
stem_ab_dbh_bin_plot <- ggplot(abund_dbh_ord, 
  aes(x = plotcode_ord, y = n, fill = dbh_bin)) +
  geom_bar(aes(linetype = est), stat = "identity", colour = "black", 
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
  scale_linetype_manual(values = c(5,1), guide = FALSE)

# Save plot
pdf(file = "img/stem_ab_dbh_bin.pdf", width = 15, height = 6)
stem_ab_dbh_bin_plot
dev.off()

# Estimate abundance of 5-10 cm in DRC using extrapolation of log abundance ----
stem_ab_lm_df <- abund_dbh_ord %>% 
  filter(est == FALSE) %>%
  group_by(plotcode, group) %>% 
  do(model = lm(log(n) ~ as.numeric(.$dbh_bin), data = .))

stem_ab_lm_df_drc <- stem_ab_lm_df %>%
  filter(group == "DRC")

extrap_drc_5_10 <- unname(sapply(stem_ab_lm_df_drc$model, function(x){
  exp(predict(x, data.frame(y=1))[1])
}))

extrap_drc_5_10_df <- data.frame(group = "DRC", 
  plotcode = stem_ab_lm_df_drc$plotcode,
  dbh_bin = "5-10", n = extrap_drc_5_10, est = TRUE)

# Bar chart of mean per site ----
abund_dbh_ord_mean <- abund_dbh_ord %>%
  dplyr::select(-ord, -plotcode_ord) %>%
  filter(est == FALSE) %>%
  rbind(., extrap_drc_5_10_df) %>%
  group_by(group, dbh_bin) %>%
  summarise(mean_n = mean(n, na.rm = TRUE),
    samples = n(),
    sem = sd(n, na.rm = TRUE) / sqrt(samples)) %>%
  mutate(est = case_when(
    group == "DRC" & dbh_bin == "5-10" ~ TRUE,
    TRUE ~ FALSE))

stem_ab_dbh_bin_plot_group <- ggplot(abund_dbh_ord_mean,
  aes(x = dbh_bin, y = mean_n, fill = group)) + 
  geom_bar(aes(linetype = est), stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = mean_n - sem, ymax = mean_n + sem), width = 0.5) + 
  facet_grid(. ~ group, scales = "free_x", space = "free") + 
  theme_classic() + 
  theme(axis.text.x=element_text(
    angle=45, 
    vjust=1, 
    hjust=1)) + 
  labs(x = "Stem diameter bin (cm)", y = "Number of stems") +
  theme(legend.position = "none") + 
  scale_fill_manual(name = "", values = big_pal, labels = c("Angola", "DRC", "Tanzania", "Mozambique")) + 
  scale_linetype_manual(values = c(1,5), guide = FALSE)

pdf(file = "img/stem_ab_dbh_bin_group.pdf", width = 8, height = 4)
stem_ab_dbh_bin_plot_group
dev.off()



# Line plot of DBH abundance distribution for each plot ----
stem_ab_lm_group_df <- abund_dbh_ord %>% 
  filter(est == FALSE) %>%
  group_by(group) %>% 
  do(model = lm(log(n) ~ as.numeric(.$dbh_bin), data = .)) %>% 
  tidy(model) %>%
  filter(term == "as.numeric(.$dbh_bin)") %>%
  dplyr::select(group, estimate, std.error)

stem_ab_lm_plot <- ggplot(filter(abund_dbh_ord, est == FALSE), 
  aes(x = as.numeric(dbh_bin), y = n, group = group, colour = group)) + 
  geom_point(position = position_jitter(width = 0.2)) + 
  stat_smooth(method = "lm") + 
  theme_classic() +
  scale_colour_manual(name = "", values = big_pal, 
    labels = c("Angola", "DRC", "Tanzania", "Mozambique")) + 
  labs(x = "Stem diameter bin (cm)", y = "Number of stems")

pdf(file = "img/stem_ab_lm.pdf", width = 8, height = 5)
stem_ab_lm_plot
dev.off()

# Proportion of stems in largest classes ----
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

# Numerical outputs ----
dbhslopebicuar <- paste0(round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "Angola",]$estimate, 2),
  "$\\pm$",
  round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "Angola",]$std.error, 3))

dbhslopedrc <- paste0(round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "DRC",]$estimate, 2),
  "$\\pm$",
  round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "DRC",]$std.error, 3))

dbhslopenham <- paste0(round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "Mozambique",]$estimate, 2),
  "$\\pm$",
  round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "Mozambique",]$std.error, 3))

dbhslopekilwa <- paste0(round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "Tanzania",]$estimate, 2),
  "$\\pm$",
  round(stem_ab_lm_group_df[stem_ab_lm_group_df$group == "Tanzania",]$std.error, 3))

fileConn <- file("include/dbh_bin_figures.tex")
writeLines(
  c(    
    paste0("\\newcommand{\\dbhslopebicuar}{", dbhslopebicuar, "}"),
    paste0("\\newcommand{\\dbhslopedrc}{", dbhslopedrc, "}"),
    paste0("\\newcommand{\\dbhslopenham}{", dbhslopenham, "}"),
    paste0("\\newcommand{\\dbhslopekilwa}{", dbhslopekilwa, "}"),
    paste0("\\newcommand{\\lmsmallstems}{", lm_format(lm_small_stems), "}"),
    paste0("\\newcommand{\\lmbigstems}{", lm_format(prop_big_group_lm), "}")
  ),
  fileConn)
close(fileConn)


