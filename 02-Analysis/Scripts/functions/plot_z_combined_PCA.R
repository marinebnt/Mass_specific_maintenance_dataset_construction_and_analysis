# beneat marine 
# 14/10/24 
#* Combine plots to shape one unique plot
#* with traits related to time or traits of different kinds
#* 

library(ggpubr)


load(file=paste0(path_plots, "/ELASMO_PCA_time.RData"))
elasmo <- plot
load(file=paste0(path_plots, "/TOT_PCA_time.RData"))
tot <- plot
load(file=paste0(path_plots, "/TELEO_PCA_time.RData"))
teleo <- plot


# Extract the legend from one of the plots
tot_plot <- tot + guides(
  shape = "none",     # remove shape legend
  fill = "none",      # remove fill legend
)
teleo_plot <- teleo + guides(
  shape = "none",     # remove shape legend
  fill = "none",      # remove fill legend
)
elasmo_plot <- elasmo + guides(
  shape = "none",     # remove shape legend
  fill = "none",      # remove fill legend
)

teleo_legend <- tot + guides(
  linetype = "none"
)

legend <- ggfun::get_legend(teleo_legend)  # from cowplot, or use get_legend() from ggpubr if preferred
legend$grobs[[2]] <- legend$grobs[[3]]

# Now arrange them in a 2x2 matrix
plotall <- ggarrange(
  tot_plot, teleo_plot,
  legend, elasmo_plot,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "", "C")
)

ggsave(plotall, file = paste0(path_plots, "/PCA_combined_time.pdf"), width=13, height=9.5)

