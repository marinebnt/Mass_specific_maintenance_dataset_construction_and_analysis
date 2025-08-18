# beneat marine 
# 14/10/24 
# PCA of the complete dataset



i_elasmo = 2 
i_opp = 6

# veccol veccnames
colAA <- c()
labelAA <- c()
pchvec <- c()
iopp <- which.max(as.data.frame(AAtot$archetypes)$K) 
colAA[iopp] = "#F8766D"
labelAA[iopp] = "Fast"
pchvec[iopp] = c(24)
ieq <- which.min(as.data.frame(AAtot$archetypes)$K) 
colAA[ieq] = "#00BFC4"
labelAA[ieq] = "Slow"
pchvec[ieq] = c(22)
iper <- c(1:3)[-c(ieq, iopp)]
colAA[iper] = "#7CAE00"
labelAA[iper] ="Intermediate"
pchvec[iper] = c(21)

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
periodic_uuid <- get_uuid(name = "Mugilidae")
periodic_pic <- get_phylopic(uuid = periodic_uuid)
elasmo_uuid <- get_uuid(name = "Squalus suckleyi")
elasmo_pic <- get_phylopic(uuid = elasmo_uuid)
opp_uuid <- get_uuid(name = "Amblygobius sphynx") # gobiidae
opp_pic <- get_phylopic(uuid = opp_uuid)
# plot(opp_pic)

# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,AXESTOREPRESENT]
rownames(datatoadd)<- NULL
labels <- labelAA
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels, pchvec=pchvec)
datatoadd$Archetypes <- labelAA
datatoadd$Archetypes <- as.factor(labelAA)
datatoadd$colAA <- colAA
order_table <- order(datatoadd$Archetypes)
order_table <- c(which((factor(labelAA)) == "Fast"),which((factor(labelAA)) == "Slow"),which((factor(labelAA)) == "Intermediate"))
datatoadd <- datatoadd[order_table,]
order_table <- c(which(levels(factor(labelAA)) == "Fast"),which(levels(factor(labelAA)) == "Slow"),which(levels(factor(labelAA)) == "Intermediate"))
datatoadd$Archetypes <- factor(datatoadd$Archetypes, levels = levels(datatoadd$Archetypes)[order_table])


# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- traits

# plot
mycol = c("darkorchid4", "cyan3", "#4575b4", "#91bfdb", "#fee090", "#fc8d59","#d73027", "white", "midnightblue")#
plot <- fviz_pca_biplot(res.pca, axes = AXESTOREPRESENT,
                        label = c("none"), # label = "none" for individual labels
                        habillage = , # color by groups
                        col.ind = as.factor(round(dataacpPLOT$colACP, 2)),
                        pointshape = 19, alpha = 0.5,
                        pointsize = 1.5, labelsize = 5,
                        arrowsize = 1.5, 
                        col.var = "midnightblue", 
                        repel = T
                        ) +  # Suppression de la légende des flèches
  geom_text_repel(data = as.data.frame(res.pca$rotation), 
                  aes(x = PC1*c(rep(14, length(PC1))), y = PC2*c(rep(7, length(PC1))), label = rownames(res.pca$rotation)),
                  size = 5, color = "midnightblue")+
  scale_color_manual(values = mycol, name = "RMR0 cluster \ncentroid") +  # Ajustement de la couleur
  geom_point(data = res.pca$x[c(which(dataacp$Species %in% c("Amblygobius sphynx", "Squalus suckleyi", "Planiliza carinata"))),],
             aes(x = PC1, y = PC2),
             pch = 8, color = "red", size = 5) +
  geom_mark_ellipse(aes(linetype = as.factor(dataacpPLOT$Class)), 
                    expand = unit(0.5, "mm"), alpha = 1, inherit.aes = T) +
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank()) +
  add_phylopic(periodic_pic, alpha = 1, x = as.numeric(datatoadd[which(labelAA=="Periodic"), 1]*1.3), 
               y = as.numeric(datatoadd[which(labelAA=="Periodic"), 2]*1.3), ysize = 0.4) +
  add_phylopic(elasmo_pic, alpha = 1, x = as.numeric(datatoadd[which(labelAA=="Equilibrium"), 1]*1.3), 
               y = as.numeric(datatoadd[which(labelAA=="Equilibrium"), 2]*1.3), ysize = 0.4) +
  add_phylopic(opp_pic, alpha = 1, x = as.numeric(datatoadd[which(labelAA=="Opportunistic"), 1]*1.3), 
               y = as.numeric(datatoadd[which(labelAA=="Opportunistic"), 2]*1.3), ysize = 0.4) +
  guides(linetype = guide_legend("Class"))

plot <- plot+ new_scale_fill()+
  geom_point(
    data = datatoadd,
    aes(
      x = x, 
      y = y,
      shape = Archetypes,       
      fill = Archetypes         
    ),
    colour = "black",          
    size = 6,
    stroke = 1
  ) +
  scale_shape_manual(values = datatoadd$pchvec) +
  scale_fill_manual(values = datatoadd$colAA)

plot



# save 
pdf(file = paste0(path_plots, "/PCA_tot_time.pdf"))
print(plot)
dev.off()



save(plot, file=paste0(path_plots, "/TOT_PCA_time.RData"))       
