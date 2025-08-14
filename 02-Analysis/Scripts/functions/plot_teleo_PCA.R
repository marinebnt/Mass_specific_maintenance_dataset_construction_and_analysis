# beneat marine 
# 14/10/24 
# PCA of the Teleostei species




# veccol veccnames
colAA <- c()
labelAA <- c()
pchvec <- c()
iopp <- which.max(as.data.frame(AAteleo$archetypes)$K) 
colAA[iopp] = "tomato"
labelAA[iopp] = "Fast"
pchvec[iopp] = c(23)
ieq <- which.min(as.data.frame(AAteleo$archetypes)$K) 
colAA[ieq] = "royalblue"
labelAA[ieq] = "Slow"
pchvec[ieq] = c(22)
iper <- c(1:3)[-c(ieq, iopp)]
colAA[iper] = "darkgreen"
labelAA[iper] ="Intermediate"
pchvec[iper] = c(21)


AXESTOREPRESENT = c(1,2)

# creating the PCA and the data frames associated
PCAteleo <- runPCA(dataplot_noelasmo, traits)
plot(PCAteleo$x[,1], PCAteleo$x[,2])
plot(PCAteleo)
res.pca <- PCAteleo
dataacp_noelasmoPLOT <- dataacp_add_colorvector(dataphylo = dataphylo_noelasmo, kclusters=7, dataacp = dataacp_noelasmo)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacp_noelasmoPLOT, AA=AAteleo, PCA=PCAteleo)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]


# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,AXESTOREPRESENT]
rownames(datatoadd)<- NULL
labels <- labelAA
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels, pchvec=pchvec)
datatoadd$Archetypes <- labelAA
datatoadd$colAA <- colAA
order_table <- order(datatoadd$Archetypes)
order_table <- c(2,3,1)
datatoadd <- datatoadd[order_table, ]

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- traits

# plot
mycol =  c("darkorchid4", "cyan3", "#4575b4", "#91bfdb", "#fee090", "#fc8d59","#d73027") 
plot <- fviz_pca_biplot(res.pca, axes = AXESTOREPRESENT,
                        label = c("none"),
                        habillage = , 
                        col.ind = as.factor(round(dataacpPLOT$colACP, 2)), 
                        arrowsize = 1.5, pointshape=19, labelsize = 5, alpha = 0.5,
                        col.var = "darkblue", repel= T)+
  geom_text_repel(data = as.data.frame(res.pca$rotation), 
                  aes(x = PC1*c(rep(13, length(PC1)-1),10), y = PC2*c(6.5,rep(5, length(PC1)-1)), label = rownames(res.pca$rotation)),
                  size = 5, color = "midnightblue")+
  scale_colour_manual(values = mycol, name="RMR0 cluster \ncentroid") + #*  # Adjust the color scale
  ggtitle(NULL) +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key = element_blank()) +
  # Manually adjust line types to match ellipses
  scale_linetype_manual(values = c("solid", "dotted")) +
  # Manually adjust the guide to show the ellipse shapes and lines
  guides(fill = "none", 
         linetype = guide_legend("Class"))
plot <- plot+ new_scale_color()+
  geom_point(data = datatoadd,
             aes(x = x, y = y,
                 shape = Archetypes, 
                 colour = Archetypes),
             size = 6, stroke = 1) +
  scale_color_manual(
    values = datatoadd$colAA, guide = "legend")

print(plot)
pdf(file = paste0(path_plots, "/PCA_teleo_time.pdf"))
plot
dev.off()



# plot
fviz_pca_biplot(res.pca,
                label = c("var"),
                habillage = , 
                fill.ind = as.factor(dataphylo_noelasmo$Species),
                col.ind=as.factor(dataphylo_noelasmo$Species),
                
                select.ind = list(name = which(dataphylo$Species %in% c("Gadus macrocephalus", "Chaenocephalus aceratus", "Lipophrys pholis",
                                                           "Engraulis australis", "Brevoortia tyrannus"))), 
                pointshape=19,
                arrowsize = 0.25, 
                col.var = "darkblue", alpha = 0.2, repel= T)+
  scale_colour_manual(values = c("darkorchid4", "#4575b4", "#91bfdb", "#fc8d59", "gold"), name="Outliers of the \ncross-validation") +
  ggtitle(NULL)+
  labs()+ 
  theme_minimal()+
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  guides(shape="none", fill="none")

save(plot, file=paste0(path_plots, "/TELEO_PCA_time.RData"))

