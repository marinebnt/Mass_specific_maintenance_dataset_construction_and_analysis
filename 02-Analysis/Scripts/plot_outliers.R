# beneat marine 
# 14/10/24 
# PCA of the Teleostei species

library(tidyverse)
setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/tests/repository_for_zenodo")
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(getwd(), "/02-Analysis/Outputs/IMAGETOT_STD_Log_BodyDepth.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))



######
# TRAIT TIME RELATED
#####

AXESTOREPRESENT=c(1,2)

# creating the PCA and the data frames associated
PCAtot <- runPCA(dataplot_2, traits_2)
plot(PCAtot$x[,1], PCAtot$x[,2])
plot(PCAtot)
res.pca <- PCAtot
# res.pca <- PCA(beta_iv1_tot)
dataacpPLOT <- dataacp_add_colorvector(dataphylo_2, kclusters=7, dataacp_2)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacpPLOT, AA=AAtot_2, PCA=PCAtot)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))] # species data and position in PCA and clusters

# loading the data for the pictures
# periodic_uuid <- get_uuid(name = "mugilidae")
# periodic_pic <- get_phylopic(uuid = periodic_uuid)
# elasmo_uuid <- get_uuid(name = "Squalidae", n=2)
# elasmo_pic <- get_phylopic(uuid = elasmo_uuid[2])
# opp_uuid <- get_uuid(name = "Gobiidae", n=5)
# opp_pic <- get_phylopic(uuid = opp_uuid[3])

# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,1:2]
rownames(datatoadd)<- NULL
labels <- c("Periodic", "Equilibrium", "Opportunistic")
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Amat", "Amax", "M", "K", "Tlvl", "Hab")

# Pseudopleuronectes_americanus
# Pseudopleuronectes_herzensteini
# Pseudopleuronectes_yokohamae
# Gilchristella_aestuaria
# Glossogobius_kokius
# Mugil_cephalus
# Mugil_curema
# Buglossidium_luteum
# Chaenocephalus_aceratus
# Gillichthys_seta
# Gillichthys_mirabilis
# Glossogobius_giuris
# Glossogobius_minutus


select.ind = list(name = dataphylo$Species[which(dataphylo$Species %in% c(#"Gillichthys seta/mirabilis", removed because tm<tmax or Lm<Loo
                                                                                            "Pseudopleuronectes herzensteini",
                                                                                            "Lipophrys pholis",
                                                                                            "Buglossidium luteum",
                                                                                            "Gilchristella aestuaria", 
                                                                                            "Chaenocephalus aceratus",
                                                                                            "Glossogobius kokius", 
                                                                                            "Brevoortia tyrannus",
                                                                                            "Mugil cephalus",
                                                                                            "Glossogobius giuris"))])
custom_points <- data.frame(res.pca$x)
plotnames <- dataphylo$Species %>% as.character()
mylabels <- sapply(plotnames, function(x) ifelse(dataphylo$Species %in% select.ind, "", x)) %>% as.vector()

p<- fviz_pca_biplot(res.pca,
                    label = c("var"),
                    pointshape=19,
                    arrowsize = 0.25, 
                    geom.ind = "none",
                    # alpha.ind=ifelse((dataphylo$Species %in% select.ind), 1,0),
                    # select.ind = list(name = select.ind),
                    col.var = "darkblue", repel= T
                    # col.ind=c(dataphylo$Species %in% c("Gillichthys_mirabilis", 
                    #                                           "Pseudopleuronectes herzensteini",
                    #                                           "Lipophrys pholis",
                    #                                           "Buglossidium luteum",
                    #                                           "Cyclothone braueri",
                    #                                           "Chaenocephalus aceratus",
                    #                                           "Glossogobius kokius", 
                    #                                           "Brevoortia tyrannus"))
)+
  scale_colour_manual(values = mycol, name="MSM cluster \ncentroid")
p <- p + geom_point(data = custom_points[which(dataphylo$Species %in% select.ind$name),], 
                    aes(x = PC1, y = PC2), 
                    size = 3, col=c("darkorchid4", "#4575b4", "#91bfdb", "#fc8d59", "gold", "red", "forestgreen", "pink", "black")) #, "forestgreen", "pink"

p <- p + geom_text_repel(data = custom_points[which(dataphylo$Species %in% select.ind$name),], 
              aes(x = PC1, y = PC2), 
              label=rownames(custom_points[which(dataphylo$Species %in% select.ind$name),]))


pdf(file = paste0(path_plots, "/outliers_teleo.pdf"), height = 8, width=6.8)
print(p)
dev.off()


# 
# 
# selected_species <- select.ind
# 
# fviz_pca_biplot(res.pca,
#                 label = "var",
#                 # habillage = dataphylo$Species,
#                 fill.ind = as.factor(dataphylo$Species),
#                 col.ind = as.factor(dataphylo$Species),  # Coloration des points sélectionnés
#                 # pointshape = ifelse(dataphylo$Species %in% selected_species, 17, 19),  # Formes différentes pour les points sélectionnés
#                 arrowsize = 0.25,
#                 col.var = "darkblue",
#                 alpha = 0.5) +
#   scale_colour_manual(values = c("red", "gray"), name = "Selected Points") +
#   theme_minimal() +
#   theme(
#     text = element_text(size = 12),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 12))
