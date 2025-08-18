# beneat marine 
# 14/10/24 
# PCA for elasmobranch species


# veccol veccnames
colAA <- c()
labelAA <- c()
pchvec <- c()
iopp <- which.max(as.data.frame(AAelasmo$archetypes)$K) 
colAA[iopp] = "#F8766D"
labelAA[iopp] = "Fast"
pchvec[iopp] = c(24)
ieq <- which.min(as.data.frame(AAelasmo$archetypes)$K) 
colAA[ieq] = "#00BFC4"
labelAA[ieq] = "Slow"
pchvec[ieq] = c(22)
iper <- c(1:3)[-c(ieq, iopp)]
colAA[iper] = "#7CAE00"
labelAA[iper] ="Intermediate"
pchvec[iper] = c(21)

AXESTOREPRESENT = c(c(-1),2)
dataacp_noteleoPLOT <- NULL
# creating the PCA and the data frames associated
PCAelasmo <- runPCA(dataplot_noteleo, traits)

# if the axes to represent are specific (ex : change an axis orientation), it will be adjusted
PCAelasmo$rotation[,abs(AXESTOREPRESENT[1])] <- sign(AXESTOREPRESENT[1])*PCAelasmo$rotation[,abs(AXESTOREPRESENT[1])]
PCAelasmo$rotation[,abs(AXESTOREPRESENT[2])] <- sign(AXESTOREPRESENT[2])*PCAelasmo$rotation[,abs(AXESTOREPRESENT[2])]
PCAelasmo$x[,abs(AXESTOREPRESENT[1])] <- sign(AXESTOREPRESENT[1])*PCAelasmo$x[,abs(AXESTOREPRESENT[1])]
PCAelasmo$x[,abs(AXESTOREPRESENT[2])] <- sign(AXESTOREPRESENT[2])*PCAelasmo$x[,abs(AXESTOREPRESENT[2])]

plot(PCAelasmo$x[,abs(AXESTOREPRESENT[1])], PCAelasmo$x[,abs(AXESTOREPRESENT[2])])
plot(PCAelasmo)
res.pca <- PCAelasmo
dataacp_noteleoPLOT <- dataacp_add_colorvector(dataphylo = dataphylo_noteleo, kclusters=7, dataacp = dataacp_noteleo)
listforplot <- preparedataforplot(numbPCA1=abs(AXESTOREPRESENT[1]), numbPCA2=abs(AXESTOREPRESENT[1]), 
                                  dataacp=dataacp_noteleoPLOT, AA=AAelasmo, PCA=PCAelasmo)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]


# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,abs(AXESTOREPRESENT)]
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
mycol =  c("darkorchid4", "cyan3", "#4575b4", "#91bfdb", "#fee090", "#fc8d59","#d73027") 
plot <- fviz_pca_biplot(res.pca, axes = abs(AXESTOREPRESENT),
                        label = c("none"),  
                        col.ind= as.factor(round(dataacpPLOT$colACP, 2)), #col_temp, #
                        arrowsize = 1.5, pointshape=19,labelsize = 5,
                        col.var = "darkblue", alpha = 0.5, repel= T)+
  scale_color_manual(values = mycol, name = "RMR0 cluster \ncentroid") +  # Ajustement de la couleur
  geom_text_repel(data = as.data.frame(res.pca$rotation), 
                  aes(x = PC1*c(rep(7, length(PC1))), y = PC2*c(rep(3, length(PC1))), label = rownames(res.pca$rotation)),
                  size = 5, color = "midnightblue")+
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank()
  ) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  # Manually adjust the guide to show the ellipse shapes and lines
  guides(fill = "none",
         linetype = guide_legend("Class"))
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


# save
pdf(file = paste0(path_plots, "/PCA_elasmo_time.pdf"), width=3.14961, height=3.14961)
print(plot)
dev.off()


save(plot, file=paste0(path_plots, "/ELASMO_PCA_time.RData"))




