# beneat marine 
# 14/10/24 
# PCA of the Teleostei species


path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(getwd(), "/02-Analysis/Outputs/IMAGE_AA_FOR_ANALYSIS.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))


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
rownames(res.pca$rotation) <- c("Amat", "Amax", "M", "K", "Tlvl", "Hab")


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

p<- fviz_pca_biplot(res.pca,
                    label = c("var"),
                    pointshape=19,
                    arrowsize = 0.25, 
                    geom.ind = "none",
                    col.var = "darkblue", repel= T

)+
  scale_colour_manual(values = mycol, name="MSM cluster \ncentroid")
p <- p + geom_point(data = custom_points[which(dataphylo$Species %in% select.ind$name),], 
                    aes(x = PC1, y = PC2), 
                    size = 3, col=c("darkorchid4", "#4575b4", "#91bfdb", "#fc8d59", "gold", "red", "forestgreen", "pink", "black")) #, "forestgreen", "pink"

p <- p + geom_text_repel(data = custom_points[which(dataphylo$Species %in% select.ind$name),], 
              aes(x = PC1, y = PC2), 
              label=dataphylo$Species[which(dataphylo$Species %in% select.ind$name)])


pdf(file = paste0(path_plots, "/outliers_teleo.pdf"), height = 8, width=6.8)
print(p)
dev.off()
p
