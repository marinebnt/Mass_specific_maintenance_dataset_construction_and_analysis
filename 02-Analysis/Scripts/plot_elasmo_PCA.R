# beneat marine 
# 14/10/24 
# PCA for elasmobranchii species

setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/tests/repository_for_zenodo")
path_phylosem_out <- paste0(getwd(), "/02-Analysis/Outputs")
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")

load(paste0(path_phylosem_out, "/IMAGETOT_STD_Log_BodyDepth.RData")) # archetypal analysis and pca outputs
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))

######
# TRAITS NOT TIME RELATED
######
  
# creating the PCA and the data frames associated
beta_iv1_tot = dataplot_noteleo[, which(colnames(dataplot_noteleo) %in%  traits)]
PCAelasmo <- runPCA(dataplot_noteleo, traits)
plot(PCAelasmo$x[,1], PCAelasmo$x[,2])
plot(PCAelasmo)
res.pca <- PCAelasmo
# res.pca <- PCA(beta_iv1_tot
dataacp_noteleoPLOT <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=4, dataacp_noteleo)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacp_noteleoPLOT, AA=AAelasmo, PCA=PCAelasmo)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))] # species data and position in PCA and clusters

# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,1:2]
rownames(datatoadd)<- NULL
labels <- c("Benthopelagic", "Pelagic", "Demersal")
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Amax", "K", "Winf", "Tlvl", "Fec", "Hab")

# plot
mycol =  c("darkorchid4", "#4575b4", "#fee090", "#fc8d59")
plot <- fviz_pca_biplot(res.pca, axes = c(1,3),
                        label = c("var"),# label = "none", # hide individual labels
                        habillage = , # color by groups
                        # palette = c("darkblue", "royalblue", "#00AFBB", "lightblue", "#E7B800", "#FC4E07", "darkred"),
                        col.ind = as.factor(round(dataacpPLOT$colACP, 2)), #dataacpPLOT$fecundity, #as.factor(dataacpPLOT$Order), # #as.factor(dataacpPLOT$Trophic.lvl), #
                        # col.ind= dataacpPLOT$c_m, 
                        arrowsize = 1.5, pointshape=19,labelsize = 5,
                        col.var = "grey14", alpha = 0.5, repel= T)+
  # scale_colour_viridis_c(direction = -1, name="Continuous\nMSM color\ngradient") +  # Adjust the color scale
  scale_colour_manual(values = mycol, name="MSM cluster \ncentroid")+
  # geom_text(data = datatoadd, aes(x = x + 0.1, y = y - 0.4, label = z),                 
  #           color = "black", size = 4.5, fontface = "bold") +
  # geom_point(data = datatoadd, aes(x = x, y = y),                 
  #            color = c("#00b159", "tomato", "royalblue"), size = 6,
  #            stroke = 1, fill = "white") +
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank(),
        # legend.position = 'none'
        ) +
  # Manually adjust line types to match ellipses
  scale_linetype_manual(values = c("solid", "dotted")) +
  # Manually adjust the guide to show the ellipse shapes and lines
  guides(fill = "none",
         linetype = guide_legend("Class"),
         color = guide_legend(override.aes = list(fill = NA)))

plot


pdf(file = paste0(path_plots, "/PCA_elasmo.pdf"), width=3.14961, height=3.14961)
plot
dev.off()

ggsave(plot, file = paste0(path_plots, "/PCA_elasmo.pdf"), height = 8, width=6.8) #, width = 80, height = 180, dpi = 600, units = "mm"


save(plot, file=paste0(path_plots, "/ELASMO_PCA.RData"))




#MORPHO






# creating the PCA and the data frames associated
beta_iv1_tot = dataplot_noteleo[, which(colnames(dataplot_noteleo) %in%  traits_morpho)]
PCAelasmo <- runPCA(dataplot_noteleo, traits_morpho)
plot(PCAelasmo$x[,1], PCAelasmo$x[,2])
plot(PCAelasmo)
res.pca <- PCAelasmo
# res.pca <- PCA(beta_iv1_tot)
dataacp_noteleoPLOT <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=4, dataacp_noteleo)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacp_noteleoPLOT, AA=AAelasmo_morpho, PCA=PCAelasmo)
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
# labels <- c("Benthopelagic", "Pelagic", "Demersal")
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Peddepth", "JawL", "Bdepth", "Bwidth")

# plot
mycol =  c("darkorchid4", "#4575b4", "#fee090", "#fc8d59")
plot <- fviz_pca_biplot(res.pca, axes = c(1,3),
                        label = c("var"),# label = "none", # hide individual labels
                        habillage = , # color by groups
                        # palette = c("darkblue", "royalblue", "#00AFBB", "lightblue", "#E7B800", "#FC4E07", "darkred"),
                        col.ind = as.factor(round(dataacpPLOT$colACP, 2)), #dataacpPLOT$fecundity, #as.factor(dataacpPLOT$Order), # #as.factor(dataacpPLOT$Trophic.lvl), #
                        # col.ind= dataacpPLOT$c_m, 
                        arrowsize = 1.5, pointshape=19,labelsize = 5,
                        col.var = "grey14", alpha = 0.5, repel= T)+
  # scale_colour_viridis_c(direction = -1, name="Continuous\nMSM color\ngradient") +  # Adjust the color scale
  scale_colour_manual(values = mycol, name="MSM cluster \ncentroid")+
  # geom_text(data = datatoadd, aes(x = x + 0.1, y = y - 0.4, label = z),                 
  #           color = "black", size = 4.5, fontface = "bold") +
  # geom_point(data = datatoadd, aes(x = x, y = y),                 
  #            color = c("#00b159", "tomato", "royalblue"), size = 6,
  #            stroke = 1, fill = "white") +
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank(),
        # legend.position = 'none'
  ) +
  # Manually adjust line types to match ellipses
  scale_linetype_manual(values = c("solid", "dotted")) +
  # Manually adjust the guide to show the ellipse shapes and lines
  guides(fill = "none",
         linetype = guide_legend("Class"),
         color = guide_legend(override.aes = list(fill = NA)))

plot


pdf(file = paste0(path_plots, "/PCA_elasmo.pdf"), width=3.14961, height=3.14961)
plot
dev.off()





######
# TRAITS TIME RELATED
######

# creating the PCA and the data frames associated
PCAelasmo <- runPCA(dataplot_noteleo_2, traits_2)
#PCAelasmo <- runPCA(dataplot_noteleo_2[, traits_2], scale.unit = TRUE)
plot(PCAelasmo$x[,1], PCAelasmo$x[,2])
plot(PCAelasmo)
res.pca <- PCAelasmo
# res.pca <- PCA(beta_iv1_tot)
dataacp_noteleoPLOT <- dataacp_add_colorvector(dataphylo_noteleo_2, kclusters=10, dataacp_noteleo_2)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacp_noteleoPLOT, AA=AAelasmo_2, PCA=PCAelasmo)
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
# labels <- c("Benthopelagic", "Pelagic", "Demersal")
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Amat", "Amax", "M", "K", "Tlvl", "Hab")

# plot
mycol =  c("darkorchid4", "#4575b4", "#fee090", "#fc8d59") #c("darkorchid4", "#4575b4", "#fee090", "#fc8d59") #c("darkorchid4", "#4575b4", "#91bfdb", "#fc8d59","#d73027") 
plot <- fviz_pca_biplot(res.pca, axes = c(1,3),
                        label = c("var"),# label = "none", # hide individual labels
                        habillage = , # color by groups
                        # palette = c("darkblue", "royalblue", "#00AFBB", "lightblue", "#E7B800", "#FC4E07", "darkred"),
                        # col.ind = dataacpPLOT$c_m, #as.factor(round(dataacpPLOT$colACP, 2)),
                        #* gradient.cols = mycol,#
                        # col.ind = as.factor(round(dataacpPLOT$colACP, 2)),      
                        col.ind= dataacpPLOT$c_m,
                        arrowsize = 1.5, pointshape=19,labelsize = 5,
                        col.var = "grey14", alpha = 0.5, repel= T)+
  # scale_colour_viridis_c(direction = -1, name="Continuous\nMSM color\ngradient") +  # Adjust the color scale
  # scale_colour_manual(values = mycol, name="MSM cluster \ncentroid")+ #*
  scale_colour_gradient(
    low = "#4575b4",
    high = "red", 
    name="MSM color gradient"
  )+
  # geom_text(data = datatoadd, aes(x = x + 0.1, y = y - 0.4, label = z),                 
  #           color = "black", size = 4.5, fontface = "bold") +
  # geom_point(data = datatoadd, aes(x = x, y = y),                 
  #            color = c("#00b159", "tomato", "royalblue"), size = 6,
  #            stroke = 1, fill = "white") +
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank(),
        # legend.position = 'none'
  ) +
  # xlim(-6,2)+
  # Manually adjust line types to match ellipses
  scale_linetype_manual(values = c("solid", "dotted")) +
  # Manually adjust the guide to show the ellipse shapes and lines
  guides(fill = "none",
         linetype = guide_legend("Class"))
         # color = guide_legend(override.aes = list(fill = NA)))

plot


pdf(file = paste0(path_plots, "/PCA_elasmo_time.pdf"), width=3.14961, height=3.14961)
plot
dev.off()

ggsave(plot, file = paste0(path_plots, "/PCA_elasmo_time.pdf")) #, width = 80, height = 180, dpi = 600, units = "mm"


save(plot, file=paste0(path_plots, "/ELASMO_PCA_time.RData"))




#MORPHO






# creating the PCA and the data frames associated
beta_iv1_tot = dataplot_noteleo[, which(colnames(dataplot_noteleo) %in%  traits_morpho)]
PCAelasmo <- runPCA(dataplot_noteleo, traits_morpho)
plot(PCAelasmo$x[,1], PCAelasmo$x[,2])
plot(PCAelasmo)
res.pca <- PCAelasmo
# res.pca <- PCA(beta_iv1_tot)
dataacp_noteleoPLOT <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=4, dataacp_noteleo)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacp_noteleoPLOT, AA=AAelasmo_morpho, PCA=PCAelasmo)
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
# labels <- c("Benthopelagic", "Pelagic", "Demersal")
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Peddepth", "JawL", "Bdepth", "Bwidth"   )

# plot
mycol = c("darkorchid4", "#4575b4", "#fee090", "#fc8d59")  #c("darkorchid4", "#4575b4", "#91bfdb", "#fc8d59","#fc8d59") 
plot <- fviz_pca_biplot(res.pca, axes = c(1,3),
                        label = c("var"),# label = "none", # hide individual labels
                        habillage = , # color by groups
                        # palette = c("darkblue", "royalblue", "#00AFBB", "lightblue", "#E7B800", "#FC4E07", "darkred"),
                        col.ind = as.factor(round(dataacpPLOT$colACP, 2)), #dataacpPLOT$fecundity, #as.factor(dataacpPLOT$Order), # #as.factor(dataacpPLOT$Trophic.lvl), #
                        arrowsize = 1.5, pointshape=19,labelsize = 6,
                        col.var = "grey14", alpha = 1, repel= T)+
  # scale_colour_viridis_d(direction = -1, name="MSM cluster \ncentroid") +  # Adjust the color scale
  scale_colour_manual(values = mycol, name="MSM cluster \ncentroid")+
  # geom_text(data = datatoadd, aes(x = x + 0.1, y = y - 0.4, label = z),                 
  #           color = "black", size = 4.5, fontface = "bold") +
  # geom_point(data = datatoadd, aes(x = x, y = y),                 
  #            color = c("#00b159", "tomato", "royalblue"), size = 6,
  #            stroke = 1, fill = "white") +
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank(),
        # legend.position = 'none'
  ) +
  # Manually adjust line types to match ellipses
  scale_linetype_manual(values = c("solid", "dotted")) +
  # Manually adjust the guide to show the ellipse shapes and lines
  guides(fill = "none",
         linetype = guide_legend("Class"),
         color = guide_legend(override.aes = list(fill = NA)))

plot


pdf(file = paste0(path_plots, "/PCA_elasmo_morpho.pdf"), width=3.14961, height=3.14961)
plot
dev.off()

