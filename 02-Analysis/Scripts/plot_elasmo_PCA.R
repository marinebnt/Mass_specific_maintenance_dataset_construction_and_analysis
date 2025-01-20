# beneat marine 
# 14/10/24 
# PCA for elasmobranchii species

path_phylosem_out <- paste0(getwd(), "/02-Analysis/Outputs")
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(path_phylosem_out, "/IMAGE_AA_FOR_ANALYSIS.RData")) # archetypal analysis and pca outputs
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))


AXESTOREPRESENT = c(1,2)

# creating the PCA and the data frames associated
PCAelasmo <- runPCA(dataplot_noteleo, traits)
plot(PCAelasmo$x[,AXESTOREPRESENT[1]], PCAelasmo$x[,AXESTOREPRESENT[2]])
plot(PCAelasmo)
res.pca <- PCAelasmo
dataacp_noteleoPLOT <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=10, dataacp_noteleo)
listforplot <- preparedataforplot(numbPCA1=1, numbPCA2=2, dataacp=dataacp_noteleoPLOT, AA=AAelasmo, PCA=PCAelasmo)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]


# preparing the data for the Archetypes identifications
datatoadd <- matrixAAinPCA[,1:2]
rownames(datatoadd)<- NULL
datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2])

# prepare the data for the arrows indentification
rownames(res.pca$rotation) <- c("Amat", "Amax", "M", "K", "Tlvl", "Hab")

# plot
mycol =  c("darkorchid4", "#4575b4", "#fee090", "#fc8d59") 
plot <- fviz_pca_biplot(res.pca, axes = AXESTOREPRESENT,
                        label = c("var"),  
                        col.ind= dataacpPLOT$c_m,
                        arrowsize = 1.5, pointshape=19,labelsize = 5,
                        col.var = "grey14", alpha = 0.5, repel= T)+
  scale_colour_gradient(
    low = "#4575b4",
    high = "red", 
    name="MSM color gradient"
  )+
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

plot


pdf(file = paste0(path_plots, "/PCA_elasmo_time.pdf"), width=3.14961, height=3.14961)
plot
dev.off()

ggsave(plot, file = paste0(path_plots, "/PCA_elasmo_time.pdf")) #, width = 80, height = 180, dpi = 600, units = "mm"


save(plot, file=paste0(path_plots, "/ELASMO_PCA_time.RData"))




