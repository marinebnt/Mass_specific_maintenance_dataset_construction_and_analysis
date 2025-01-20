# beneat marine 
# 14/10/24 
#* Combine plots to shape one unique plot
#* with traits related to time or traits of different kinds
#* 
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(getwd(), "/02-Analysis/Outputs/IMAGE_AA_FOR_ANALYSIS.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))


load(file=paste0(path_plots, "/ELASMO_PCA_time.RData"))
elasmo <- plot
load(file=paste0(path_plots, "/TELEO_PCA_time.RData"))
teleo <- plot
load(file=paste0(path_plots, "/TOT_PCA_time.RData"))
tot <- plot

plotall <- ggarrange(tot, teleo, NA, elasmo, labels = c("A","B","", "C"))

ggsave(file = paste0(path_plots, "/PCA_combined_time.pdf"), width=13, height=9.5)

