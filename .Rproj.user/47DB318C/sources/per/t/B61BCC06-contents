setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/tests/repository_for_zenodo")
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(getwd(), "/02-Analysis/Outputs/IMAGETOT_STD_Log_BodyDepth.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))



#***************************
#*
#*   TRAIT TIME RELATED
#*
pPCAelasmo_2 -> pPCAelasmo_tmp
pPCAteleo_2 -> pPCAteleo_tmp
pPCAtot_2 -> pPCAtot_tmp

#*
#*
#*



##############
## Life history traits 
##############


# treeteleo <- treeforpPCA(dataphylo_noelasmo)
rownames(dataplot_noelasmo) <- stringr::str_replace(dataphylo_noelasmo$Species, " ", "_")
# pPCAteleo_tmp <- runpPCA(dataplot_noelasmo, traits, tree=treeteleo)
eigen <- colSums(pPCAteleo_tmp$Evec)
eigenval <- c(abs(colSums(pPCAteleo_tmp$Evec))/sum(abs(eigen)))*100
n <- length(eigenval)
numbPCA1 = which.max(eigenval)
eigenval[numbPCA1] <- 0
numbPCA2 = which.max(eigenval)
plot(pPCAtot_tmp$L[,numbPCA1], -pPCAtot_tmp$L[,numbPCA2])
text(pPCAtot_tmp$L[,numbPCA1], -pPCAtot_tmp$L[,numbPCA2], rownames(pPCAtot_tmp$L))
plot(pPCAtot_tmp$S[,numbPCA1], -pPCAtot_tmp$S[,numbPCA2])
dataacp_noelasmo <- dataacp_add_colorvector(dataphylo_noelasmo, kclusters=7, dataacp_noelasmo)[[1]]
listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_noelasmo, AA=AAteleo_2, PCA=pPCAteleo_tmp)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]
pozarrow=3
a <- plotPCAclean(c("royalblue", "#00b159", "tomato"), c("EQU", "PER", "OPP"), pPCAteleo_tmp)
a
pdf(file = paste0(path_plots, "/pPCA_teleo.pdf"), height = 8, width=6.8)
a
dev.off()

# treeelasmo <- treeforpPCA(dataphylo_noteleo)
rownames(dataplot_noteleo) <- stringr::str_replace(dataphylo_noteleo$Species, " ", "_")
# pPCAelasmo_tmp <- runpPCA(dataplot_noteleo, traits, tree=treeelasmo)
eigen <- colSums(pPCAelasmo_tmp$Evec)
eigenval <- c(abs(colSums(pPCAelasmo_tmp$Evec))/sum(abs(eigen)))*100
n <- length(eigenval)
numbPCA1 = which.max(eigenval)
eigenval[numbPCA1] <- 0
numbPCA2 = which.max(eigenval)
dataacpPLOT <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=7, dataacp_noteleo)[[1]]
listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacpPLOT, AA=AAelasmo_2, PCA=pPCAelasmo_tmp)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]*100
pozarrow = 1
a <- plotPCAclean(c("royalblue", "tomato", "#00b159"), c("Demersal", "Pelagic", "Benthopelagic"), pPCAelasmo_tmp)
a
pdf(file = paste0(path_plots, "/pPCA_elasmo.pdf"), height = 8, width=6.8)
a
dev.off()

# treetot <- treeforpPCA(dataphylo)
rownames(dataplot) <- stringr::str_replace(dataphylo$Species, " ", "_")
# pPCAtot_tmp <- runpPCA(dataplot, traits, tree=treetot)
eigen <- colSums(pPCAtot_tmp$Evec)
eigenval <- c(abs(colSums(pPCAtot_tmp$Evec))/sum(abs(eigen)))*100
n <- length(eigenval)
numbPCA1 = which.max(eigenval)
eigenval[numbPCA1] <- 0
numbPCA2 = which.max(eigenval)
dataacp_notot <- dataacp_add_colorvector(dataphylo, kclusters=7, dataacp)[[1]]
listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_notot, AA=AAtot_2, PCA=pPCAtot_tmp)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]
arrows = 2
a <- plotPCAclean(c("tomato","royalblue",  "#00b159"), c("OPP","EQU","PER"), pPCAtot_tmp)
a
pdf(file = paste0(path_plots, "/pPCA_tot.pdf"), height = 8, width=6.8)
a
dev.off()
a
###########
# MORPHOMETRICS
###########
# 
# treeteleo <- treeforpPCA(dataphylo_noelasmo)
# rownames(dataplot_noelasmo) <- stringr::str_replace(dataphylo_noelasmo$Species, " ", "_")
# # pPCAteleo_tmp_morpho <- runpPCA(dataplot_noelasmo, traits_morpho, tree=treeteleo)
# eigen <- colSums(pPCAteleo_tmp_morpho$Evec)
# eigenval <- c(abs(colSums(pPCAteleo_tmp_morpho$Evec))/sum(abs(eigen)))
# n <- length(eigenval)
# numbPCA1 = which.max(eigenval)
# eigenval[numbPCA1] <- 0
# numbPCA2 = which.max(eigenval)
# dataacp_noelasmo_morpho <- dataacp_add_colorvector(dataphylo_noelasmo, kclusters=7, dataacp_noelasmo)
# listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_noelasmo_morpho, AA=AAtot_morpho, PCA=pPCAteleo_tmp_morpho)
# rotation = listforplot[[1]]
# matrixAAinpPCA = listforplot[[2]]
# dataacpPLOT = listforplot[[3]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
# eigenval = listforplot[[4]]
# plotPCAclean(veccol_tot_morpho, vecAA_tot_morpho, pPCAteleo_tmp_morpho)
# 
# treeelasmo <- treeforpPCA(dataphylo_noteleo)
# rownames(dataplot_noteleo) <- stringr::str_replace(dataphylo_noteleo$Species, " ", "_")
# #pPCAelasmo_tmp_morpho <- runpPCA(dataplot_noteleo, traits_morpho, tree=treeelasmo)
# eigen <- colSums(pPCAelasmo_tmp_morpho$Evec)
# eigenval <- c(abs(colSums(pPCAelasmo_tmp_morpho$Evec))/sum(abs(eigen)))
# n <- length(eigenval)
# numbPCA1 = which.max(eigenval)
# eigenval[numbPCA1] <- 0
# numbPCA2 = which.max(eigenval)
# dataacp_noteleo_morpho <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=7, dataacp_noteleo)
# listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_noteleo_morpho, AA=AAtot_morpho, PCA=pPCAelasmo_tmp_morpho)
# rotation = listforplot[[1]]
# matrixAAinpPCA = listforplot[[2]]
# dataacpPLOT = listforplot[[3]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
# eigenval = listforplot[[4]]
# plotPCAclean(veccol_tot_morpho, vecAA_tot_morpho, pPCAelasmo_tmp_morpho)
# 
# treetot <- treeforpPCA(dataphylo)
# rownames(dataplot) <- stringr::str_replace(dataphylo$Species, " ", "_")
# # pPCAtot_tmp_morpho <- runpPCA(dataplot, traits_morpho, tree=treetot)
# eigen <- colSums(pPCAtot_tmp$Evec)
# eigenval <- c(abs(colSums(pPCAtot_tmp$Evec))/sum(abs(eigen)))
# n <- length(eigenval)
# numbPCA1 = which.max(eigenval)
# eigenval[numbPCA1] <- 0
# numbPCA2 = which.max(eigenval)
# dataacp_notot_morpho <- dataacp_add_colorvector(dataphylo, kclusters=7, dataacp)
# listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_notot_morpho, AA=AAtot_morpho, PCA=pPCAtot_tmp_morpho)
# rotation = listforplot[[1]]
# matrixAAinpPCA = listforplot[[2]]
# dataacpPLOT = listforplot[[3]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
# eigenval = listforplot[[4]]
# # plotPCAclean(veccol_tot_morpho, vecAA_tot_morpho, pPCAtot_tmp_morpho)






#***************************
#*
#*   TRAIT NOT TIME RELATED
#*
pPCAelasmo -> pPCAelasmo_tmp
pPCAteleo -> pPCAteleo_tmp
pPCAtot -> pPCAtot_tmp

#*
#*
#*



##############
## Life history traits 
##############


# treeteleo <- treeforpPCA(dataphylo_noelasmo)
rownames(dataplot_noelasmo) <- stringr::str_replace(dataphylo_noelasmo$Species, " ", "_")
# pPCAteleo_tmp <- runpPCA(dataplot_noelasmo, traits, tree=treeteleo)
eigen <- colSums(pPCAteleo_tmp$Evec)
eigenval <- c(abs(colSums(pPCAteleo_tmp$Evec))/sum(abs(eigen)))
n <- length(eigenval)
numbPCA1 = which.max(eigenval)
eigenval[numbPCA1] <- 0
numbPCA2 = which.max(eigenval)
#plot(pPCAtot_tmp$x[,1], -pPCAtot_tmp$x[,2])
dataacp_noelasmo <- dataacp_add_colorvector(dataphylo_noelasmo, kclusters=7, dataacp_noelasmo)
listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_noelasmo, AA=AAtot, PCA=pPCAteleo_tmp)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]
pozarrow=3
plotPCAclean(c("royalblue", "#00b159", "tomato"), c("EQU", "PER", "OPP"), pPCAteleo_tmp)

# treeelasmo <- treeforpPCA(dataphylo_noteleo)
rownames(dataplot_noteleo) <- stringr::str_replace(dataphylo_noteleo$Species, " ", "_")
# pPCAelasmo_tmp <- runpPCA(dataplot_noteleo, traits, tree=treeelasmo)
eigen <- colSums(pPCAelasmo_tmp$Evec)
eigenval <- c(abs(colSums(pPCAelasmo_tmp$Evec))/sum(abs(eigen)))*100
n <- length(eigenval)
numbPCA1 = which.max(eigenval)
eigenval[numbPCA1] <- 0
numbPCA2 = which.max(eigenval)
dataacpPLOT <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=7, dataacp_noteleo)
listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacpPLOT, AA=AAtot, PCA=pPCAelasmo_tmp)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]*100
plotPCAclean(c("royalblue", "tomato", "#00b159"), c("Demersal", "Pelagic", "Benthopelagic"), pPCAelasmo_tmp)

# treetot <- treeforpPCA(dataphylo)
rownames(dataplot) <- stringr::str_replace(dataphylo$Species, " ", "_")
# pPCAtot_tmp <- runpPCA(dataplot, traits, tree=treetot)
eigen <- colSums(pPCAtot_tmp$Evec)
eigenval <- c(abs(colSums(pPCAtot_tmp$Evec))/sum(abs(eigen)))*100
n <- length(eigenval)
numbPCA1 = which.max(eigenval)
eigenval[numbPCA1] <- 0
numbPCA2 = which.max(eigenval)
dataacp_notot <- dataacp_add_colorvector(dataphylo, kclusters=7, dataacp)
listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_notot, AA=AAtot, PCA=pPCAtot_tmp)
rotation = listforplot[[1]]
matrixAAinPCA = listforplot[[2]]
dataacpPLOT = listforplot[[3]]
dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
eigenval = listforplot[[4]]
arrows = 2
plotPCAclean(c("tomato","royalblue",  "#00b159"), c("OPP","EQU","PER"), pPCAtot_tmp)


###########
# MORPHOMETRICS
###########



# treeteleo <- treeforpPCA(dataphylo_noelasmo)
# rownames(dataplot_noelasmo) <- stringr::str_replace(dataphylo_noelasmo$Species, " ", "_")
# # pPCAteleo_tmp_morpho <- runpPCA(dataplot_noelasmo, traits_morpho, tree=treeteleo)
# eigen <- colSums(pPCAteleo_tmp_morpho$Evec)
# eigenval <- c(abs(colSums(pPCAteleo_tmp_morpho$Evec))/sum(abs(eigen)))
# n <- length(eigenval)
# numbPCA1 = which.max(eigenval)
# eigenval[numbPCA1] <- 0
# numbPCA2 = which.max(eigenval)
# dataacp_noelasmo_morpho <- dataacp_add_colorvector(dataphylo_noelasmo, kclusters=7, dataacp_noelasmo)
# listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_noelasmo_morpho, AA=AAtot_morpho, PCA=pPCAteleo_tmp_morpho)
# rotation = listforplot[[1]]
# matrixAAinpPCA = listforplot[[2]]
# dataacpPLOT = listforplot[[3]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
# eigenval = listforplot[[4]]
# plotPCAclean(veccol_tot_morpho, vecAA_tot_morpho, pPCAteleo_tmp_morpho)
# 
# treeelasmo <- treeforpPCA(dataphylo_noteleo)
# rownames(dataplot_noteleo) <- stringr::str_replace(dataphylo_noteleo$Species, " ", "_")
# #pPCAelasmo_tmp_morpho <- runpPCA(dataplot_noteleo, traits_morpho, tree=treeelasmo)
# eigen <- colSums(pPCAelasmo_tmp_morpho$Evec)
# eigenval <- c(abs(colSums(pPCAelasmo_tmp_morpho$Evec))/sum(abs(eigen)))
# n <- length(eigenval)
# numbPCA1 = which.max(eigenval)
# eigenval[numbPCA1] <- 0
# numbPCA2 = which.max(eigenval)
# dataacp_noteleo_morpho <- dataacp_add_colorvector(dataphylo_noteleo, kclusters=7, dataacp_noteleo)
# listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_noteleo_morpho, AA=AAtot_morpho, PCA=pPCAelasmo_tmp_morpho)
# rotation = listforplot[[1]]
# matrixAAinpPCA = listforplot[[2]]
# dataacpPLOT = listforplot[[3]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
# eigenval = listforplot[[4]]
# plotPCAclean(veccol_tot_morpho, vecAA_tot_morpho, pPCAelasmo_tmp_morpho)
# 
# treetot <- treeforpPCA(dataphylo)
# rownames(dataplot) <- stringr::str_replace(dataphylo$Species, " ", "_")
# # pPCAtot_tmp_morpho <- runpPCA(dataplot, traits_morpho, tree=treetot)
# eigen <- colSums(pPCAtot_tmp$Evec)
# eigenval <- c(abs(colSums(pPCAtot_tmp$Evec))/sum(abs(eigen)))
# n <- length(eigenval)
# numbPCA1 = which.max(eigenval)
# eigenval[numbPCA1] <- 0
# numbPCA2 = which.max(eigenval)
# dataacp_notot_morpho <- dataacp_add_colorvector(dataphylo, kclusters=7, dataacp)
# listforplot <- preparedataforplotp(numbPCA1, numbPCA2, dataacp=dataacp_notot_morpho, AA=AAtot_morpho, PCA=pPCAtot_tmp_morpho)
# rotation = listforplot[[1]]
# matrixAAinpPCA = listforplot[[2]]
# dataacpPLOT = listforplot[[3]]
# dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
# eigenval = listforplot[[4]]
# # plotPCAclean(veccol_tot_morpho, vecAA_tot_morpho, pPCAtot_tmp_morpho)
# 
