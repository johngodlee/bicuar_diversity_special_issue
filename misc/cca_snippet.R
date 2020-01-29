stem_ab_env_cca <- cca(stem_ab_mat ~ mat + map + mat_cov + map_cov, plot_clim_s_all[,5:9])

stem_ab_env_cca_plot <- plot(stem_ab_env_cca, type="p")
stem_ab_env_cca_plot$biplot <- as.data.frame(stem_ab_env_cca_plot$biplot * attributes(stem_ab_env_cca_plot$biplot)$arrow.mul)
stem_ab_env_cca_plot$species <- as.data.frame(stem_ab_env_cca_plot$species)
stem_ab_env_cca_plot$sites <- as.data.frame(stem_ab_env_cca_plot$sites)
stem_ab_env_cca_plot$biplot$var <- row.names(stem_ab_env_cca_plot$biplot)
stem_ab_env_cca_plot$species$species <- row.names(stem_ab_env_cca_plot$species)
stem_ab_env_cca_plot$sites$plotcode <- row.names(stem_ab_env_cca_plot$sites)
stem_ab_env_cca_plot$sites$plotcode <- plot_group$plotcode
stem_ab_env_cca_plot$sites$group <- plot_group$group

env_cca_plot <- ggplot() + 
  geom_hline(aes(yintercept = 0), linetype = 2) + 
  geom_vline(aes(xintercept = 0), linetype = 2) +  
  geom_point(data = stem_ab_env_cca_plot$species,
    aes(x = CCA1, y = CCA2), size = 1) +
  geom_point(data = stem_ab_env_cca_plot$sites,
    aes(x = CCA1, y = CCA2, fill = group), 
    shape = 23, size = 3) + 
  geom_segment(data = stem_ab_env_cca_plot$biplot, 
    aes(xend = CCA1, yend = CCA2, x = 0, y = 0), 
    arrow = arrow(length = unit(0.03, "npc"))) + 
  geom_text(data = stem_ab_env_cca_plot$biplot,
    aes(x = CCA1, y = CCA2, label = var),
    nudge_y = 0.2) + 
  coord_equal() +
  theme_classic()