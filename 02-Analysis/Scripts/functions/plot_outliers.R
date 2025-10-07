# beneat marine 
# 14/10/24 
# PCA of the Teleostei species



dataseterrorE <- read.csv( paste0(path_CV,  "dataseterror", "SEM", "all", "c_m", ".csv"))

AXESTOREPRESENT=c(1,-2)

# creating the PCA and the data frames associated
PCAtot <- runPCA(dataplot, traits)
plot(PCAtot$x[,1], PCAtot$x[,2])
plot(PCAtot)
res.pca <- PCAtot

res.pca$rotation[, abs(AXESTOREPRESENT)] <- cbind(sign(AXESTOREPRESENT[1])*res.pca$rotation[, abs(AXESTOREPRESENT)[1]], 
                                                  sign(AXESTOREPRESENT[2])*res.pca$rotation[, abs(AXESTOREPRESENT)[2]])
res.pca$x[, abs(AXESTOREPRESENT)] <- cbind(sign(AXESTOREPRESENT[1])*res.pca$x[, abs(AXESTOREPRESENT)[1]], 
                                           sign(AXESTOREPRESENT[2])*res.pca$x[, abs(AXESTOREPRESENT)[2]])

# res.pca <- PCA(beta_iv1_tot)
dataacpPLOT <- dataacp_add_colorvector(dataphylo, kclusters=7, dataacp)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacpPLOT, AA=AAtot, PCA=PCAtot)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]
dataacpPLOT$boolean <- ifelse(stringr::str_replace(dataacp$Species, " ", "_") %in% dataseterrorE$lab, 1, 0)


dataacpPLOT$Inference <- dataacpPLOT$boolean
dataacpPLOT$Name <- dataacpPLOT$boolean
k=0
for (i in stringr::str_replace(dataacp$Species, " ", "_")){
  k=k+1
  if(i %in% dataseterrorE$lab){
    ID_bool <- which(dataseterrorE$lab == i)
    dataacpPLOT$Inference[k] <- ifelse(dataseterrorE$errorpercent[ID_bool]>0, "understimated", "overestimated")
    dataacpPLOT$Name[k] <- dataseterrorE$ID[ID_bool]
  }
  else{
    dataacpPLOT$Inference[k] <- ""
  }
}

# preparing the data for the Archetypes identifications



# veccol veccnames
colAA <- c()
labelAA <- c()
pchvec <- c()
iopp <- which.max(as.data.frame(AAtot$archetypes)$K) 
colAA[iopp] =  "#F8766D"
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
datatoadd <- cbind(sign(AXESTOREPRESENT[1])*matrixAAinPCA[, abs(AXESTOREPRESENT)[1]], 
                   sign(AXESTOREPRESENT[2])*matrixAAinPCA[, abs(AXESTOREPRESENT)[2]])
rownames(datatoadd)<- NULL
labels <- labelAA
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels, pchvec=pchvec)
datatoadd$Archetypes <- as.factor(labelAA)
datatoadd$colAA <- colAA
order_table <- order(datatoadd$Archetypes)
order_table <- c(which((factor(labelAA)) == "Fast"),which((factor(labelAA)) == "Slow"),which((factor(labelAA)) == "Intermediate"))
datatoadd <- datatoadd[order_table,]
order_table <- c(which(levels(factor(labelAA)) == "Fast"),which(levels(factor(labelAA)) == "Slow"),which(levels(factor(labelAA)) == "Intermediate"))
datatoadd$Archetypes <- factor(datatoadd$Archetypes, levels = levels(datatoadd$Archetypes)[order_table])


# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- traits

species_outliers <- stringr::str_replace(dataseterrorE$label, "_", " ")[which(dataseterrorE$Booleand_outliers==1)]
select.ind = list(name = dataphylo$Species[which(dataphylo$Genus %in% "Spratelloides")])
select.ind = list(name = dataphylo$Species[which(dataphylo$Species %in% species_outliers)])
custom_points <- data.frame(res.pca$x)
mycol = c("darkorchid4", "cyan3", "#4575b4", "#91bfdb", "#fee090", "#fc8d59","#d73027", "white", "midnightblue")#

p<- fviz_pca_biplot(res.pca,
                    label = c("var"),
                    pointshape=19,
                    arrowsize = 0.25, 
                    geom.ind = "none",
                    col.var = "darkblue", repel= T
)
p <- p+ 
  geom_point(
    data = datatoadd,
    aes(
      x = x, 
      y = y,
      shape = Archetypes
    ),
    fill="white",
    colour = "black",          
    size = 6,
    stroke = 1
  ) +
  scale_shape_manual(values = datatoadd$pchvec) 
p <- p +new_scale(new_aes = "shape")+
  scale_fill_manual(values = mycol, name="RMR0 cluster \ncentroid")+
  scale_colour_manual(values = mycol, name="RMR0 cluster \ncentroid")+
  scale_shape_manual(values=c(NA,25,23))+
  geom_point(data = custom_points,
             aes(x = PC1, y = PC2, colour= as.factor(round(dataacpPLOT$colACP, 2)),
                        fill = as.factor(round(dataacpPLOT$colACP, 2))), 
                    size = 3, pch=21, 
                    alpha=0.03)

p <- p +
  geom_point(data = custom_points,
                             aes(x = PC1, y = PC2, fill = as.factor(round(dataacpPLOT$colACP, 2)), 
                                 pch=(dataacpPLOT$Inference)),
                              size = 4, alpha = ifelse(dataacpPLOT$boolean, 1, 0), 
                              colour="black")

# p <- p +
#   geom_point(data = custom_points,
#              aes(x = PC1, y = PC2, fill = as.factor(round(dataacpPLOT$colACP, 2))),
#              size = 3, alpha = ifelse(dataacpPLOT$boolean, 1, 0),
#              pch=21, colour="black")
  # geom_point(data = custom_points[which(dataphylo$Species %in% select.ind$name),],
  #                   aes(x = PC1, y = PC2),
  #                   size = 3,
  #                   col=c("darkorchid4", "magenta", "#4575b4", "#91bfdb", "#fc8d59", "gold", "red",
  #                         "forestgreen", "pink", "black")[as.numeric(as.factor(dataphylo$Genus[which(dataphylo$Species %in% select.ind$name)]))])
                    #, col=c("darkorchid4", "#4575b4", "#fc8d59", "gold", "red", "forestgreen")) #, "forestgreen", "pink" "#91bfdb", , "pink", "black"


p <- p + geom_text_repel(data = custom_points[which(dataphylo$Species %in% select.ind$name),], 
              aes(x = PC1, y = PC2), 
              max.overlaps = 20,
              label=dataacpPLOT$Name[which(dataphylo$Species %in% select.ind$name)])+
    guides(fill = guide_legend(override.aes=list(shape = 21)), shape=guide_legend(title = "Inferred compared \nwith observed", reverse=TRUE))
  # new_scale_fill()+
p <- p +
  scale_fill_manual(values = mycol, name="RMR0 cluster \ncentroid") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black")))



# +
  # scale_fill_manual(values = datatoadd$colAA)


pdf(file = paste0(path_plots, "outliers_tot.pdf"), height = 8, width=6.8)
print(p)
dev.off()
  
