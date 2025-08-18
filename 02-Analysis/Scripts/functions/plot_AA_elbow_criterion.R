# beneat marine 
# 14/10/24 
# Plot elbow criterion 

library(archetypes)

#############
model <- "TLstdmeca"
OUTPUT <- paste0("Outputs/", model)
path_plots <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")
#############

# produce the data to plot
# beta_iv1 = dataplot[, which(colnames(dataplot) %in%  traits)]
# repAA_tot<-stepArchetypes(data=beta_iv1, k=1:6, verbose=TRUE, nrep=30) # this is long to run

# save.image(file=paste0(path_plots, "/ELBOW_tot.RData"))

# plot
load(paste0(path_plots, "/ELBOW_tot.RData"))

rss_data <- rss(repAA_tot)

num_components <- 1:nrow(rss_data)
rss_mean <- rowMeans(rss_data, na.rm = T)
rss_sd <- apply(rss_data, 1, sd)  # or se = sd / sqrt(ncol(rss_data))


elbow_df <- data.frame(
  Components = num_components,
  RSS_Mean = rss_mean,
  RSS_SD = rss_sd
)

elbow <- ggplot(elbow_df, aes(x = Components, y = RSS_Mean)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = RSS_Mean - RSS_SD, ymax = RSS_Mean + RSS_SD), width = 0.2) +
  theme_minimal() +
  labs(title = "Elbow Plot of RSS vs. Number of Components",
       x = "Number of Components",
       y = "Mean RSS ± SD")

ggsave(elbow, file=paste0(path_plots, "/plot_elbow.pdf"), width=5, height=7)


#
# a <- rss(repAA_tot)^2
# plot(a/sum(a), type = "b",
#      xlab = "Principal Component",
#      ylab = "Percentage of Variance Explained")
#
#
#
#
#
# pdf(file = paste0(path_plots, "/elbow_tot.pdf"))
# print(screenplot(rss(repAA_tot)))
# dev.off()
#
#
#
# beta_iv1 = dataplot[, which(colnames(dataplot_noteleo) %in%  traits)]
# repAA_elasmo<-stepArchetypes(data=beta_iv1, k=1:kmax, verbose=TRUE, nrep=30)
#
# pdf(file = paste0(pathoutput, "/elbow_elasmo"))
# print(screenplot(repAA_elasmo))
# dev.off()
#
#
#
# beta_iv1 = dataplot[, which(colnames(dataplot_noelasmo) %in%  traits)]
# repAA_teleo<-stepArchetypes(data=beta_iv1, k=1:kmax, verbose=TRUE, nrep=30)
#
# pdf(file = paste0(pathoutput, "/elbow_teleo"))
# print(screenplot(repAA_teleo))
# dev.off()
