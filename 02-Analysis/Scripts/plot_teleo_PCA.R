# beneat marine 
# 14/10/24 
# PCA of the Teleostei species


path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(getwd(), "/02-Analysis/Outputs/IMAGE_AA_FOR_ANALYSIS.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))



AXESTOREPRESENT = c(1,2)

# creating the PCA and the data frames associated
PCAteleo <- runPCA(dataplot_noelasmo, traits)
plot(PCAteleo$x[,1], PCAteleo$x[,2])
plot(PCAteleo)
res.pca <- PCAteleo
dataacp_noelasmoPLOT <- dataacp_add_colorvector(dataphylo_noelasmo, kclusters=6, dataacp_noelasmo)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacp_noelasmoPLOT, AA=AAteleo, PCA=PCAteleo)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]


# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,1:2]
rownames(datatoadd)<- NULL
labels <- c("Opportunistic", "Equilibrium", "Periodic")
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels)

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Amat", "Amax", "M", "K", "Tlvl", "Hab")

# plot
mycol =  c("darkorchid4", "cyan3", "#4575b4", "#91bfdb", "#fee090", "#fc8d59","#d73027") 
plot <- fviz_pca_biplot(res.pca, axes = AXESTOREPRESENT,
                        label = c("var"),
                        habillage = , 
                        col.ind = as.factor(round(dataacpPLOT$colACP, 2)), 
                        arrowsize = 1.5, pointshape=19, labelsize = 5, alpha = 0.5,
                        col.var = "darkblue", repel= T)+
  scale_colour_manual(values = mycol, name="MSM cluster \ncentroid") + #*  # Adjust the color scale
  geom_point(data = datatoadd, aes(x = x, y = y), pch=c(22,23,21),                 
             fill = c("tomato", "royalblue", "#00b159"), size = 6,
             stroke = 1) +
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
plot


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

