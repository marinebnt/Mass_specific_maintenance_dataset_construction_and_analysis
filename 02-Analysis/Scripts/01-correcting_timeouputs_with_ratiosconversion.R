
setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/tests/repository_for_zenodo")
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")
load(paste0(getwd(), "/02-Analysis/Outputs/IMAGETOT_STD_Log_BodyDepth.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))








path_phylosem_out <- paste0(getwd(), "/01-Simulations/Outputs/phylosem_output")
path_output_genus <- paste0(getwd(), "/01-Simulations/Outputs/dataset_creation_output/dataset_for_phylosem/output_tot_stdmorpho")
path_analysis_out <- paste0(getwd(), "/02-Analysis/Outputs")


# load(paste0(path_analysis_out, "/IMAGETOT_STD_Log_BodyDepth.RData"))



# load(paste0(path_phylosem_out, "/IMAGETOT_STD_Log_BodyDepth.RData"))

dataphylo <- read.csv(paste0(path_phylosem_out, "/output_TLstdmecapsemFINALtot.csv"))
datagenus <- read.csv(paste0(path_output_genus, "/dataset_phylosemLOG.csv"))
datagenus$na <- is.na(datagenus$c_m)

# add a column with a number / habitat
dataphylo <- addhabitatcolumn(dataphylo)
dataphylo$Species <- str_replace(dataphylo[,"label"], "_", " ")
dataphylo <- left_join(dataphylo, datagenus[, c("Species", "na")])

# convert ratio of traits into trait value *************
dataphylo$Lm <- ratioTOtrait(dataphylo$Lm, dataphylo$Loo, IS_OPERATOR_Loo=TRUE)
dataphylo$tm <- ratioTOtrait(dataphylo$tm, dataphylo$M, IS_OPERATOR_Loo=FALSE) # 1/ M = longevity

#remove data that is unrealistic
# dataphylo <- dataphylo[-which(dataphylo$Lm>dataphylo$Loo),]
# dataphylo <- dataphylo[-which(dataphylo$tm>dataphylo$tmax),]

#Standardize
dataphylo[, -which(colnames(dataphylo) %in% c("Species", "X", "label", "hhabtot", "osmose", "na"))] <-
  decostand(dataphylo[, -which(colnames(dataphylo) %in% c("Species", "X", "label", "hhabtot", "osmose", "na"))], "standardize")

##### Create datasets to be used ######
dataphylo <- left_join(dataphylo[, - which(colnames(dataphylo) %in% c("X", "label"))],
                       datagenus[, which(colnames(datagenus) %in% c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode"))])
colnames(dataphylo)[which(colnames(dataphylo) == "hhabtot")]<- "Habitat"
colnames(dataphylo)[which(colnames(dataphylo) == "fecundity")]<- "Fecundity"
colnames(dataphylo)[which(colnames(dataphylo) == "TLDiet")]<- "Trophic.lvl"
colnames(dataphylo)[which(colnames(dataphylo) == "tmax")]<- "Age.max"
colnames(dataphylo)[which(colnames(dataphylo) == "tm")]<- "Age.mat"
colnames(dataphylo)[which(colnames(dataphylo) == "M")]<- "Mortality"
colnames(dataphylo)[which(colnames(dataphylo) == "Woo")]<- "Weight.Inf"


dataplot  <- data.frame(dataphylo$K, dataphylo$Mortality, dataphylo$Temperature,
                        dataphylo[, -which(colnames(dataphylo) %in%
                                             c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
dataacp   <- data.frame(dataphylo)

#********
# SELECTION CRITERIA = Elasmo
dataphylo_noelasmo <- dataphylo[-which(dataphylo$Class %in% c("Elasmobranchii")),]
dataplot_noelasmo  <- data.frame(dataphylo_noelasmo$K, dataphylo_noelasmo$Mortality, dataphylo_noelasmo$Temperature,
                                 dataphylo_noelasmo[, -which(colnames(dataphylo_noelasmo) %in%
                                                               c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
dataacp_noelasmo   <- data.frame(dataphylo_noelasmo)
#colnames(dataplot_noelasmo)[c(1,2,3)] <- c("K", "M", "Temperature")


#*****
# SELECTION CRITERIA = Teleo
dataphylo_noteleo <- dataphylo[-which(dataphylo$Class %in% c("Teleostei")),]
dataplot_noteleo  <- data.frame(dataphylo_noteleo$K, dataphylo_noteleo$Mortality, dataphylo_noteleo$Temperature,
                                dataphylo_noteleo[, -which(colnames(dataphylo_noteleo) %in%
                                                             c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
dataacp_noteleo   <- data.frame(dataphylo_noteleo)
#colnames(dataplot_noelasmo)[c(1,2,3)] <- c("K", "M", "Temperature")


#********
# SELECTION CRITERIA = Teleo Pelagic
dataphylo_noteleo_pela <- dataphylo[which(dataphylo$habitatpelagic > 0.5),]
dataplot_noteleo_pela  <- data.frame(dataphylo_noteleo_pela$K, dataphylo_noteleo_pela$Mortality, dataphylo_noteleo_pela$Temperature,
                                     dataphylo_noteleo_pela[, -which(colnames(dataphylo_noteleo_pela) %in%
                                                                       c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
dataacp_noteleo_pela   <- data.frame(dataphylo_noteleo_pela)
#colnames(dataplot_noelasmo)[c(1,2,3)] <- c("K", "M", "Temperature")


############################## TRAITS THV ##########################

traits = c("Weight.Inf", "Age.max", "Fecundity", "K", "Trophic.lvl", "Habitat")
traits_2 = c("Age.mat", "Age.max", "Mortality", "K", "Trophic.lvl", "Habitat")
veccol_teleo <- c("royalblue","darkgreen","tomato")
veccol_elasmo <- c("orange", "purple4", "pink")
veccol_tot <- c("royalblue", "tomato", "darkgreen")


vecAA_teleo <- c("Per", "Opp", "Equ")
vecAA_elasmo <- c("+ Age max", "Pelagic", "Demersal")
vecAA_elasmo_pela <- c("+K", "Pelagic", "-fec/+demer")

vecAA_eq <- c("a", "b", "c")

kmax=4

AAteleo <- runAA(dataplot_noelasmo, traits, kmax)
AAelasmo <- runAA(dataplot_noteleo, traits, kmax)
AAtot <- runAA(dataplot, traits, kmax)
# AAteleo2<- runAA(dataplot_noelasmo, traits, kmax=2)

dataphylo_noelasmo_2 <- dataphylo_noelasmo
dataphylo_noteleo_2  <- dataphylo_noteleo
dataphylo_2          <- dataphylo
dataplot_noelasmo_2 <- dataplot_noelasmo
dataplot_noteleo_2  <- dataplot_noteleo
dataplot_2          <- dataplot
dataacp_noelasmo_2 <- dataacp_noelasmo
dataacp_noteleo_2  <- dataacp_noteleo
dataacp_2          <- dataacp
AAteleo_2  <- runAA(dataplot_noelasmo_2, traits_2, kmax)
AAelasmo_2 <- runAA(dataplot_noteleo_2, traits_2, kmax)
AAtot_2    <- runAA(dataplot_2, traits_2, kmax)





treeteleo <- treeforpPCA(dataphylo_noelasmo)
rownames(dataplot_noelasmo) <- stringr::str_replace(dataphylo_noelasmo$Species, " ", "_")
pPCAteleo <- runpPCA(dataplot_noelasmo, traits, tree=treeteleo)

treeelasmo <- treeforpPCA(dataphylo_noteleo)
rownames(dataplot_noteleo) <- stringr::str_replace(dataphylo_noteleo$Species, " ", "_")
pPCAelasmo <- runpPCA(dataplot_noteleo, traits, tree=treeelasmo)

treetot <- treeforpPCA(dataphylo)
rownames(dataplot) <- stringr::str_replace(dataphylo$Species, " ", "_")
pPCAtot <- runpPCA(dataplot, traits, tree=treetot)


treeelasmo <- treeforpPCA(dataphylo_noteleo_2)
rownames(dataplot_noteleo_2) <- stringr::str_replace(dataphylo_noteleo_2$Species, " ", "_")
pPCAelasmo_2 <- runpPCA(dataplot_noteleo_2, traits_2, tree=treeelasmo)

treetot <- treeforpPCA(dataphylo_2)
rownames(dataplot_2) <- stringr::str_replace(dataphylo_2$Species, " ", "_")
pPCAtot_2 <- runpPCA(dataplot_2, traits_2, tree=treetot)

treeteleo <- treeforpPCA(dataphylo_noelasmo_2)
rownames(dataplot_noelasmo_2) <- stringr::str_replace(dataphylo_noelasmo_2$Species, " ", "_")
pPCAteleo_2 <- runpPCA(dataplot_noelasmo_2, traits_2, tree=treeteleo)









save.image(paste0(getwd(), "/02-Analysis/Outputs/IMAGETOT_STD_Log_BodyDepth_with_ratios_corrected.RData"))