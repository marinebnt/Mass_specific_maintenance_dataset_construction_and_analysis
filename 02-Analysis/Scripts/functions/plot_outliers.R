# beneat marine 
# 14/10/24 
# PCA of the Teleostei species



dataseterrorE <- read.csv2( paste0(path_CV,  "/dataseterror", "SEM", "all", "c_m", ".csv"))

AXESTOREPRESENT=c(1,2)

# creating the PCA and the data frames associated
PCAtot <- runPCA(dataplot, traits)
plot(PCAtot$x[,1], PCAtot$x[,2])
plot(PCAtot)
res.pca <- PCAtot
# res.pca <- PCA(beta_iv1_tot)
dataacpPLOT <- dataacp_add_colorvector(dataphylo, kclusters=7, dataacp)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacpPLOT, AA=AAtot, PCA=PCAtot)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]

# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,1:2]
rownames(datatoadd)<- NULL
labels <- c("Periodic", "Equilibrium", "Opportunistic")
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- traits

species_outliers <- stringr::str_replace(dataseterrorE$label, "_", " ")[which(dataseterrorE$Booleand_outliers==1)]
select.ind = list(name = dataphylo$Species[which(dataphylo$Species %in% species_outliers)])
custom_points <- data.frame(res.pca$x)

p<- fviz_pca_biplot(res.pca,
                    label = c("var"),
                    pointshape=19,
                    arrowsize = 0.25, 
                    geom.ind = "none",
                    col.var = "darkblue", repel= T

)+
  scale_colour_manual(values = mycol, name="RMR0 cluster \ncentroid")
p <- p + geom_point(data = custom_points, 
                    aes(x = PC1, y = PC2), 
                    size = 3, 
                    col="lightgrey",
                    alpha=0.05)


p <- p + geom_point(data = custom_points[which(dataphylo$Species %in% select.ind$name),], 
                    aes(x = PC1, y = PC2), 
                    size = 3, 
                    col=c("darkorchid4", "magenta", "#4575b4", "#91bfdb", "#fc8d59", "gold", "red", 
                          "forestgreen", "pink", "black")[as.numeric(as.factor(dataphylo$Genus[which(dataphylo$Species %in% select.ind$name)]))])
                    #, col=c("darkorchid4", "#4575b4", "#fc8d59", "gold", "red", "forestgreen")) #, "forestgreen", "pink" "#91bfdb", , "pink", "black"


p <- p + geom_text_repel(data = custom_points[which(dataphylo$Species %in% select.ind$name),], 
              aes(x = PC1, y = PC2), 
              max.overlaps = 20,
              label=dataphylo$Species[which(dataphylo$Species %in% select.ind$name)])

p

pdf(file = paste0(path_plots, "/outliers_tot.pdf"), height = 8, width=6.8)
print(p)
dev.off()

