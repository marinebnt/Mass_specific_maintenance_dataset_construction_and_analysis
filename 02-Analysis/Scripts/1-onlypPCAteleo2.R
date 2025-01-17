
setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/tests/repository_for_zenodo")
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(getwd(), "/02-Analysis/Outputs/IMAGETOT_STD_Log_BodyDepth.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))

treeteleo <- treeforpPCA(dataphylo_noelasmo_2)
rownames(dataplot_noelasmo_2) <- stringr::str_replace(dataphylo_noelasmo_2$Species, " ", "_")
pPCAteleo_2 <- runpPCA(dataplot_noelasmo_2, traits_2, tree=treeteleo)

save.image(paste0(getwd(), "/02-Analysis/Outputs/IMAGETOT_STD_Log_BodyDepth.RData"))