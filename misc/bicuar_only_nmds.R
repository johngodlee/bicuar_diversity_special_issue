# Compare diversity within the Bicuar plots only ----
stem_ab_mat_bicuar <- stem_ab(trees_list_fil$bicuar, trees_list_fil$bicuar$plotcode, 
  trees_list_fil$bicuar$species_binomial)

plot_vec_bicuar <- stem_ab_mat_bicuar$plotcode
stem_ab_mat_bicuar_clean <- stem_ab_mat_bicuar %>% 
  dplyr::select(-plotcode)
NMDS.scree(stem_ab_mat_bicuar_clean)

# NMDS of species composition across plots
##' Use Bray Curtis distance as it's not biased by null values in the matrix
bicuar_stem_ab_nmds <- metaMDS(stem_ab_mat_bicuar_clean, distance = "bray", 
  try = 500, trymax = 500, k = 3, autotransform = FALSE)

stressplot(bicuar_stem_ab_nmds)
##' original dissimilarities are  well preserved in the reduced dimensions

# Extract site (plot) scores from NMDS analysis
bicuar_plot_scores <- data.frame(scores(bicuar_stem_ab_nmds))  
bicuar_plot_scores$plot <- plot_vec_bicuar

# Extract species scores from NMDS analysis
bicuar_species_scores <- as.data.frame(scores(bicuar_stem_ab_nmds, "species")) 
bicuar_species_scores$species_binomial <- rownames(bicuar_species_scores)

# Plot extracted scores in ggplot2
bicuar_nmds_plot <- ggplot() + 
  geom_point(data = bicuar_plot_scores,
    aes(x = NMDS1, y = NMDS2, fill = plot), shape = 23, size = 3) +
  geom_point(data = bicuar_species_scores,
    aes(x = NMDS1, y = NMDS2)) + 
  geom_label_repel(data = bicuar_plot_scores,
    aes(x = NMDS1, y = NMDS2, label =  plot, colour = plot),
    label.padding = 0.15, box.padding = 0.5) +
  # geom_label_repel(data = filter(bicuar_species_scores, NMDS1 > 0.45),
  #   aes(x = NMDS1, y = NMDS2, label = species_binomial),
  #   label.padding = 0.15, box.padding = 0.5, alpha = 0.6)
  coord_equal() +
  theme_classic() + 
  theme(legend.position = "none")

pdf(file = paste0("img/bicuar_nmds.pdf"), width = 8, height = 6)
bicuar_nmds_plot
dev.off()

##' Plots 15, 13, and 9 are different from those in the rest of Bicuar
##' These are Baikiaea plurijuga / Baphia massaiensis woodlands.
