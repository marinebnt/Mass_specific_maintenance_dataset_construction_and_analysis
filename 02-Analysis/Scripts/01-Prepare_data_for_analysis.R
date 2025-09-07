# beneat marine 
# 14/10/24 
############################
#* 
#* THIS SCRIPT GOAL IS TO RUN PCAs AND phylogenetic PCAs AND Archetypal Analysis
#* ITS OUTPUTS ARE USED AS INPUTS FOR THE PLOT CONSTRUCTIONS OF THE FOLLOWING SCRIPTS
#* 






#*****
####Preparing dataset####
#******
K_M_LOG = T  # **** ############### ARE K AND M in LOG ????
############################## Prepare dataset ##########################

path_phylosem_out <- paste0(getwd(), "/01-Dataset_construction/", OUTPUT_phylo)
path_output_genus <- paste0(getwd(), "/01-Dataset_construction/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV/output_tot_stdmorpho")
path_analysis_out <- paste0(getwd(), "/02-Analysis/", OUTPUT)
dir.create(path_analysis_out)


dataphylo <- read.csv(paste0(path_phylosem_out, "/output_SEMpsemFINALtot.csv"))
# if (K_M_LOG){
#   datagenus <- read.csv(paste0(path_output_genus, "/dataset_phylosemLOG.csv"))  
  IS_LOG_M = T
# }
# if (!K_M_LOG) {
  datagenus <- read.csv(paste0(path_output_genus, "/dataset_phylosem.csv")) 
#   IS_LOG_M = F
# }

datagenus$na <- is.na(datagenus$c_m)

# add a column with a number / habitat
dataphylo <- addhabitatcolumn(dataphylo)
dataphylo$Species <- str_replace(dataphylo[,"label"], "_", " ")
dataphylo <- left_join(dataphylo, datagenus[, c("Species", "na")])

# convert ratio of traits into trait value *************
dataphylo$Lm <- ratioTOtrait(dataphylo$Lm, dataphylo$Loo, IS_OPERATOR_Loo=TRUE, IS_LOG_M=IS_LOG_M)
dataphylo$tm <- ratioTOtrait(dataphylo$tm, dataphylo$M, IS_OPERATOR_Loo=FALSE, IS_LOG_M=IS_LOG_M) # 1/ M = longevity

#remove data that is unrealistic
if(length(which(dataphylo$Lm>dataphylo$Loo))>0) {dataphylo <- dataphylo[-which(dataphylo$Lm>dataphylo$Loo),]}
if(length(which(dataphylo$tm>dataphylo$tmax))>0) {dataphylo <- dataphylo[-which(dataphylo$tm>dataphylo$tmax),]}

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
# dataphylo_noelasmo <- dataphylo[-which(dataphylo$Class %in% c("Elasmobranchii")),]
# dataplot_noelasmo  <- data.frame(dataphylo_noelasmo$K, dataphylo_noelasmo$Mortality, dataphylo_noelasmo$Temperature,
#                                  dataphylo_noelasmo[, -which(colnames(dataphylo_noelasmo) %in%
#                                                                c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
# dataacp_noelasmo   <- data.frame(dataphylo_noelasmo)


#*****
# SELECTION CRITERIA = Teleo
# dataphylo_noteleo <- dataphylo[-which(dataphylo$Class %in% c("Teleostei")),]
# dataplot_noteleo  <- data.frame(dataphylo_noteleo$K, dataphylo_noteleo$Mortality, dataphylo_noteleo$Temperature,
#                                 dataphylo_noteleo[, -which(colnames(dataphylo_noteleo) %in%
#                                                              c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
# dataacp_noteleo   <- data.frame(dataphylo_noteleo)


#********
# SELECTION CRITERIA = Teleo Pelagic
# dataphylo_noteleo_pela <- dataphylo[which(dataphylo$habitatpelagic > 0.5),]
# dataplot_noteleo_pela  <- data.frame(dataphylo_noteleo_pela$K, dataphylo_noteleo_pela$Mortality, dataphylo_noteleo_pela$Temperature,
#                                      dataphylo_noteleo_pela[, -which(colnames(dataphylo_noteleo_pela) %in%
#                                                                        c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
# dataacp_noteleo_pela   <- data.frame(dataphylo_noteleo_pela)


############################## RUN Archetypal analysis ##########################


# veccol_teleo <- c("royalblue","darkgreen","tomato")
# veccol_elasmo <- c("orange", "purple4", "pink")
veccol_tot <- c("royalblue", "tomato", "darkgreen")


vecAA_teleo <- c("Per", "Opp", "Equ")
# vecAA_elasmo <- c("+ Age max", "Pelagic", "Demersal")
# vecAA_elasmo_pela <- c("+K", "Pelagic", "-fec/+demer")

vecAA_eq <- c("a", "b", "c")
kmax=3




# AAteleo <- runAA(dataplot_noelasmo, traits, kmax)
# AAelasmo <- runAA(dataplot_noteleo, traits, kmax)
AAtot <- runAA(dataplot, traits, kmax)




############################## SAVE outputs ##########################

save.image(paste0(path_analysis_out, "/IMAGE_AA_FOR_ANALYSIS.RData"))

save(AAtot, dataphylo, dataacp, dataplot, datagenus, 
          file=paste0(path_analysis_out, "/IMAGE_AA_CONSTRUCTED.RData"))

