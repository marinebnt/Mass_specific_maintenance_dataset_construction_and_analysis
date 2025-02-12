# beneat marine 
# 14/10/24 
# PCA of the complete dataset



OUTPUT = "Outputs_WITHOUTUNITCVLOGNOTKM"
labelsAA <- c("Equilibrium", "Periodic", "Opportunistic") 
labelsAA <- c("Opportunistic", "Equilibrium", "Periodic")
colAA <- c("royalblue", "#00b159", "tomato")
colAA <- c("tomato", "royalblue", "#00b159")

load(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_FOR_ANALYSIS.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
path_phylosem_out <- paste0(getwd(), "/02-Analysis/", OUTPUT)
path_plots <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")


AXESTOREPRESENT = c(1,2)


# creating the PCA and the data frames associated
PCAtot <- runPCA(dataplot, traits)
plot(PCAtot$x[,AXESTOREPRESENT[1]], PCAtot$x[,AXESTOREPRESENT[2]])
plot(PCAtot)
res.pca <- PCAtot
plot(PCAtot)
dataacp_totPLOT <- dataacp_add_colorvector(dataphylo = dataphylo, kclusters=7, dataacp = dataacp)
clusterscentroids <- dataacp_totPLOT[2][[1]]
listforplot <- preparedataforplot(numbPCA1=AXESTOREPRESENT[1], numbPCA2=AXESTOREPRESENT[2], 
                                  dataacp=dataacp_totPLOT, AA=AAtot, PCA=PCAtot)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]

# loading the data for the pictures
periodic_uuid <- get_uuid(name = "Sarpa salpa")
periodic_pic <- get_phylopic(uuid = periodic_uuid)
elasmo_uuid <- get_uuid(name = "Squalus suckleyi")
elasmo_pic <- get_phylopic(uuid = elasmo_uuid)
opp_uuid <- get_uuid(name = "Amblyeleotris guttata") # gobiidae
opp_pic <- get_phylopic(uuid = opp_uuid)
plot(opp_pic)

# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,AXESTOREPRESENT]
rownames(datatoadd)<- NULL
labels <- labelsAA
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Amat", "Amax", "M", "K", "Tlvl", "Hab")

# plot
mycol = c("darkorchid4", "cyan3", "#4575b4", "#91bfdb", "#fee090", "#fc8d59","#d73027", "white", "midnightblue")#
plot <- fviz_pca_biplot(res.pca, axes = AXESTOREPRESENT,
                        label = c("var"), # label = "none" for individual labels
                        habillage = , # color by groups
                        col.ind = as.factor(round(dataacpPLOT$colACP, 2)),
                        pointshape = 19, alpha = 0.5,
                        pointsize = 1.5, labelsize = 5,
                        arrowsize = 1.5, 
                        col.var = "midnightblue", 
                        repel = T
                        ) +  # Suppression de la légende des flèches
  geom_mark_ellipse(aes(fill = as.factor(dataacpPLOT$Class), 
                        linetype = as.factor(dataacpPLOT$Class)), 
                    expand = unit(0.5, "mm"), alpha = 0) +
  scale_colour_manual(values = mycol, name = "MSM cluster \ncentroid") +  # Ajustement de la couleur
  geom_point(data = datatoadd, aes(x = x, y = y), pch=c(23,21,22),                 
             fill = colAA, size = 6, 
             stroke = 1) +
  geom_point(data = res.pca$x[c(which(dataacp$Species %in% c("Amblyeleotris guttata", "Squalus suckleyi", "Sarpa salpa"))),],
             aes(x = PC1, y = PC2),
             pch = 8, color = "red", size = 5) +
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank()) +
  add_phylopic(periodic_pic, alpha = 1, x = as.numeric(datatoadd[which(labelsAA=="Periodic"), 1]*1.3), 
               y = as.numeric(datatoadd[which(labelsAA=="Periodic"), 2]*1.3), ysize = 0.4) +
  add_phylopic(elasmo_pic, alpha = 1, x = as.numeric(datatoadd[which(labelsAA=="Equilibrium"), 1]*1.3), 
               y = as.numeric(datatoadd[which(labelsAA=="Equilibrium"), 2]*1.3), ysize = 0.4) +
  add_phylopic(opp_pic, alpha = 1, x = as.numeric(datatoadd[which(labelsAA=="Opportunistic"), 1]*1.3), 
               y = as.numeric(datatoadd[which(labelsAA=="Opportunistic"), 2]*1.3), ysize = 0.4) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  guides(linetype = guide_legend("Class"), fill="none", 
         point = guide_legend(override.aes=list(fill=c(mycol[1:4], "white", "white"))))
plot

pdf(file = paste0(path_plots, "/PCA_tot_time.pdf"))
plot
dev.off()



save(plot, file=paste0(path_plots, "/TOT_PCA_time.RData"))       
