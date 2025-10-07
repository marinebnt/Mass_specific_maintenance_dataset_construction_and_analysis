##*********************************************************************************
## TO USE PHYLOSEM I NEED A FEW THINGS :
## -> 1 - CREATE A DATASET WITH THE OXYGEN CONSUMPTION PARAMETERS
## - 2 - DEDUCE A RELATIONSHIP BETWEEN OXYGEN DATA AND WEIGHT AND TEMPERATURE
## - 3 - CREATE A SEPERATE DATASET WITH DATA EXTRACTED FROM FISHBASE
##        AND ADD C_M AND EPS_M DATA TO THE DATASET 
## - 4 - RUN PHYLOSEM
##*********************************************************************************

# Files needed : completerfishbase.csv

#### Load packages + path #####
library(rfishbase)
library(dplyr)
library(car)
library(stringr)
library(ape)
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(ggpubr)
library(worrms)


# setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/repository_for_zenodo")
path <- paste0(getwd(), "/01-Dataset_construction/Scripts")
pathoutput <- paste0(getwd(), "/01-Dataset_construction/Outputs/dataset_creation_output")
pathoutputplot <- paste0("01-Dataset_construction/Outputs/dataset_creation_output/plot_01")
dir.create(pathoutputplot, recursive = T)
##*************************
## TAXO DATA : organise + fill in ######
##*************************

##
ordersnames <- 
c("Clupeiformes" , "Carangiformes", "Scorpaeniformes", "Labriformes", "Perciformes", "Gobiiformes" , "Anguilliformes" , "Batrachoidiformes" , 
  "Pleuronectiformes", "Mugiliformes", "Lophiiformes", "Gadiformes", "Atheriniformes", "Scombriformes", "Aulopiformes", 
  "Beloniformes", "Tetraodontiformes", "Carcharhiniformes", "Rajiformes")


## extract taxo and species data
taxo_tot   <- load_taxa()
totspe_full<- species()
totspe_full     <- full_join(totspe_full, taxo_tot[,-which(colnames(taxo_tot)%in% c("Subfamily", "Genus"))], by=c("SpecCode"))
totspe_full     <- totspe_full[-which(totspe_full$SpecCode %in% c(0,1)),] # garbage
totspe_full$DemersPelag[totspe_full$DemersPelag %in% c("pelagic", "pelagic-oceanic", "pelagic-neritic", "bathypelagic")] <- c("pelagic")
totspe_full$DemersPelag[totspe_full$DemersPelag %in% c("reef-associated")] <- c("benthopelagic")
totspe_full$DemersPelag[totspe_full$DemersPelag %in% c("bathydemersal")] <- c("demersal")
# View(totspe_full[which(totspe_full$SpecCode %in% names(table(totspe_full$SpecCode))[table(totspe_full$SpecCode)==2]),])

# We are keeping only the Marine and Brack species, except if fresh+salwater species are known to be migratory with a longer time spent in marine ecosystem
sum_fresh  <- rowsum(as.data.frame(totspe_full[,c("Fresh", "Saltwater", "Brack")]), (totspe_full$Genus))
sum_fresh  <- sum_fresh[which(sum_fresh$Fresh<apply(sum_fresh[, c("Saltwater", "Brack")], 1, max)),]

fresh <- totspe_full$Genus[which(c(totspe_full$Fresh==1 & totspe_full$Saltwater==1 ))]
migratory <- (table(totspe_full$Genus, totspe_full$AnaCat))#as.data.frame(totspe[which(totspe$Fresh == 1), c("AnaCat")])
migratory <- migratory[which(rownames(migratory) %in% fresh),]
non_migr <- migratory[, "non-migratory"] + migratory[, 'potamodromous'] + migratory[, "catadromous"] # identify which species to exclude
migratory[,c("unknown")] <- 0
migratory[,c("non-migratory")] <- 0
migratory[,c("potamodromous")] <- ifelse(migratory[,c("potamodromous")]>0, c(-1000), 0) 
migratory[,c("catadromous")] <- ifelse(migratory[,c("potamodromous")]>0, c(-1000), 0) 
migratory[,c(" ")] <- 0
migr <- rowSums(migratory)
Fresh_genus_to_keep <- names(migr[(migr - non_migr )>0])
totspe_which_to_keep <- totspe_full[which(totspe_full$Genus %in% Fresh_genus_to_keep),]
# View(totspe_which_to_keep[which(totspe_which_to_keep$SpecCode %in% names(table(totspe_which_to_keep$SpecCode))[table(totspe_which_to_keep$SpecCode)==2]),])
# 1st : dataset with only brack and saltwater species
totspe_no_fresh <- totspe_full[which(c(totspe_full$Saltwater==1 & totspe_full$Fresh==0 & totspe_full$Brack==0) | c(totspe_full$Saltwater==1 & totspe_full$Fresh==0 & totspe_full$Brack==1)), ]
# View(totspe_no_fresh[which(totspe_no_fresh$SpecCode %in% names(table(totspe_no_fresh$SpecCode))[table(totspe_no_fresh$SpecCode)==2]),])
# 2nd : add the migratory freshwater species
totspe <- rbind(totspe_no_fresh, totspe_which_to_keep)
# verify I did not duplicate lines and remove them
# View(totspe[which(totspe$SpecCode %in% names(table(totspe$SpecCode))[table(totspe$SpecCode)==2]),])
totspe <- totspe[!duplicated(totspe$SpecCode),]
dim(totspe)

# keeping only useful variables
totspe     <- totspe[, names(totspe) %in% c("SpecCode", "DemersPelag", "Saltwater", "Brack", names(taxo_tot))]
totspe     <- totspe[-which(totspe$Genus %in% "Fundulus"),]

## fill the family name NAs
fillNAs <- read.csv2(paste0(getwd(), "/01-Dataset_construction/Inputs/completerfishbase.csv"))
famnames <- totspe$Family
for (i in seq_along(famnames)){
  if (is.na(famnames[i]) & !totspe$SpecCode[i]%in%c(1, 0)){
    ID    <- totspe$SpecCode[i]
    locID <- which(fillNAs$code==ID)
    famnames[i] <- fillNAs$fam[locID]
  }
}
totspe$Family <- famnames


## taxo_tot() data : many missing orders that are refereed to as "Eupercaria/misc" (only the ones from osmose)
## here we are replacing the missing data with the help of Worms package :
##******************************************************************************************
# 
# 1- get AphiaID
ordpercif        <- which(grepl("/", totspe$Order, fixed=TRUE))
data_missing_order <- totspe[ordpercif,]
frame_missing_order <- data.frame(Family = names(table(data_missing_order$Family)))
# 2- removing sources of errors for worrms
frame_missing_order[which(frame_missing_order$Family == "Scaridae"), "Order"] <- "Eupecaria"
frame_missing_order[which(frame_missing_order$Family == "Eulophiidae"), "Order"] <- "Perciformes"
frame_missing_order[which(frame_missing_order$Family == "Latilidae"), "Order"] <- "Eupecaria"
frame_missing_order[which(frame_missing_order$Family == "Neozoarcidae"), "Order"] <- "Perciformes"
frame_missing_order[which(frame_missing_order$Family == "Opisthocentridae"), "Order"] <- "Perciformes"
# 3- assign order to family
k=0
for (i in frame_missing_order$Family){ # loop to identify order for each family 
  k=k+1
  if(!is.na(frame_missing_order$Order[k])){next}
  data <- wm_children_(name=i)
  frame_missing_order$Order[k] = stringr::str_remove(data$order[1], pattern = " incertae sedis")
  if(length(table(data$order))>1){cat(i, "more than one order possible : ", names(table(data$order)))}
}
# 4 - assign order to full dataset
k=0
for (i in totspe$Family){
  k=k+1
  if (i %in% frame_missing_order$Family){
    totspe$Order[k] <- frame_missing_order$Order[which(frame_missing_order$Family == i)]
  }
}
##*********************************************************************





##****************************
##   OXYGEN DATA  :  merge 3 datasets   ######
##****************************

# * about the units : 
# * 
# * Wong 2021 :    resting metabolism # removed : included in Gravel 2024 already
# * Fishbase :     selected routine + standard + NA metabolism
# * Clark 2025 :   resting metabolism
# * Ikeda 2016 :   routine metabolism
# * Killen 2016b : resting metabolism
# * Gravel 2024 :  resting metabolism
# * 
# * Wong :     Temperature(K) : measurement temperature -> conversion in °C :  -273,15
# * Fishbase : Temperature (°C) : mean water temperature during the experiment
# * Clark :    Temperature (°C) : experimental temperature
# * Ikeda :    Temperature (°C) : in situ temperature
# * Killen :   Temperature (°C) : experimental temperature
# * Gravel :   Temperature (°C) : experimental temperature
# * 
# * Wong :     Weight (g) : measurement body mass -> conversion in kg : *10^-3
# * Fishbase : Weight of the test animal (g) or the mean weight if many
# individuals per experiment -> conversion in kg (*10^-3)
# * Clark :    Weight (g) : wet body mass -> conversion in kg *10^-3
# * Ikeda :    Weight (mg) : wet mass is being used here => needs to be
# *10^-6 to be in kg
# * Killen :   Weight (g) : wet body mass -> conversion in kg (*10^-3)
# * Gravel :   Weight (g) : wet body mass -> conversion in kg (*10^-3)
# * 
# * Wong : Watt =>  Oxygen consumption mg O2 . h^-1
# * Fishbase : Oxygen consumption (mg O2 . kg^-1 . h^-1) => needs to be
# * Weight(kg)-> mg O2 . h^-1
# * Clark :Watt =>  Oxygen consumption mg O2 . h^-1
# * Ikeda :    Routine respiration (microL O2 . h^-1)  => needs to be converted.
# Mass Volumic : 1.308*10^3 mg/L, so :
# (http://wiki.scienceamusante.net/index.php?title=dioxyg%C3%A8ne)
# *                                        1- microL-> *10^-6 -> L
# *                                        2- L     -> *1.308*10^3-> mg O2. h^-1
# * Killen :   mg O2/h ( # not the value adjusted to 1kg and 15°C)
# * Gravel :Watt =>  Oxygen consumption mg O2 . h^-1


##**********
#Deal with 2 datasets other than Fishbase : 
#extract data from Clarke and Johnston 1999 AND IOkeda et al 2016 and clean it 
oxygen25  <- read.csv(paste0(getwd(), "/01-Dataset_construction/Inputs/Clarke & Johnston 2025/metabolic_rates.csv"), skip = c(25))
oxygen25  <- oxygen25[which(oxygen25$Class == "Osteichthyes"),]
oxygen16  <- read.csv2(paste0(getwd(), "/01-Dataset_construction/Inputs/Ikeda16_ConsoO2.csv"))[1:102,]
oxygen16b <- read.csv(paste0(getwd(), "/01-Dataset_construction/Inputs/Killenetal16Ecological.csv"))
oxygen16b <- oxygen16b[-which(is.na(oxygen16b$RMRmass)),]
# oxygen21  <- read.csv('01-Dataset_construction/Inputs/Wong et al.,2021/SampleSizeDataset.csv') # already included in Gravel 2024
oxygen24  <- read.csv("01-Dataset_construction/Inputs/Gravel et al., 2024/data/FULL_Rmax_LHT_MR_ms.csv", header=TRUE, fileEncoding = 'UTF-8') %>% as.data.frame() 
oxygen24 <-  oxygen24 %>%
  filter(MRtype == 'RMR' & (DFtype == 'MeanT' | DFtype == 'Both')) %>%
  dplyr::select(-c("MRtype", "DFtype"))

# change units according to the units expected on fishbase
# clarke and johnston : first all in Watt
oxygen25$MR..d.[which(oxygen25$Unit..MR...e. == "mW")]      <-  
  oxygen25$MR..d.[which(oxygen25$Unit..MR...e. == "mW")]*10^(-3)
oxygen25$Unit..MR...e.[which(oxygen25$Unit..MR...e. == "mW")]      <-  "W"
oxygen25$MR..d.respi <- convert_watt_to_respi(oxygen25$MR..d.)
oxygen25$Mb..b.kg    <- as.numeric(oxygen25$Mb..b.)*10^-3
oxygen25$Genus       <- str_split_fixed(oxygen25$Species, " ", n=2)[,1]

oxygen16$kgWM        <- as.numeric(oxygen16$mgWM)*(10^-6)
oxygen16$R_L         <- as.numeric(oxygen16$R)*10^-6
oxygen16$OxygenCons  <- oxygen16$R_L*1.308*10^3

oxygen16b$RMRmasskg  <- oxygen16b$RMRmass*10^-3
oxygen16b$Species    <- stringr::str_replace(oxygen16b$species, "_", " ")

oxygen24$Species     <- oxygen24$ScientificName
oxygen24$MaxBodyMass_KG <- oxygen24$MRMassGrams * 10^-3
oxygen24$OxygenCons  <- convert_watt_to_respi(oxygen24$WholeOrganismMRWatts)
# oxygen21$Species <- oxygen21$ScientificName
# oxygen21$MaxBodyMass_KG <- oxygen21$MaxBodyMass_G * 10^-3
# oxygen21$MeasurementTemp_Celcius <- oxygen21$MeasurementTemp_Kelvin -273.15
# oxygen21$OxygenCons <- convert_watt_to_respi(oxygen21$WholeOrganismRMR_Watts)


# add species code to the data + complement the dataset that is not complete
oxygen16b    <- ox16b_correction(oxygen16b)
idsp16b      <- which(oxygen16b$Species %in%  taxo_tot$Species)
ox16bnospe   <- oxygen16b[-idsp16b,] #puting on the side the species without a complete species name
oxygen16b    <- oxygen16b[idsp16b,]
matchingspeccode16bnb     <- which(taxo_tot$Species %in% oxygen16b$Species)
matchingspeccode16b       <- taxo_tot[matchingspeccode16bnb, c("Species", "SpecCode")]
oxygen16b    <- full_join(oxygen16b, matchingspeccode16b)

oxygen16    <- ox16_correction(oxygen16)
idsp16      <- which(oxygen16$Species %in%  taxo_tot$Species)
ox16nospe   <- oxygen16[-idsp16,] #puting on the side the species without a complete species name
oxygen16    <- oxygen16[idsp16,]
matchingspeccode16nb     <- which(taxo_tot$Species %in% oxygen16$Species)
matchingspeccode16       <- taxo_tot[matchingspeccode16nb, c("Species", "SpecCode")]
oxygen16    <- full_join(oxygen16, matchingspeccode16)

oxygen25    <- ox99_correction(oxygen25) # this way are missing only misspecified only the species for which : no species specified
idsp99      <- which(oxygen25$Species %in%  taxo_tot$Species)
ox99nospe   <- oxygen25[-idsp99,] #puting on the side the species without a complete species name
oxygen25    <- oxygen25[idsp99,]
matchingspeccode99nb     <- which(taxo_tot$Species %in% oxygen25$Species)
matchingspeccode99       <- taxo_tot[matchingspeccode99nb, c("Species", "SpecCode")]
oxygen25    <- full_join(oxygen25, matchingspeccode99)

oxygen24    <- ox24_correction(oxygen24) # this way are missing only misspecified only the species for which : no species specified
idsp24      <- which(oxygen24$Species %in%  taxo_tot$Species)
ox24nospe   <- oxygen24[-idsp24,] #puting on the side the species without a complete species name
oxygen24    <- oxygen24[idsp24,]
matchingspeccode24nb     <- which(taxo_tot$Species %in% oxygen24$Species)
matchingspeccode24       <- taxo_tot[matchingspeccode24nb, c("Species", "SpecCode")]
oxygen24    <- full_join(oxygen24, matchingspeccode24)
# oxygen21    <- ox21_correction(oxygen21) # this way are missing only misspecified only the species for which : no species specified
# idsp21      <- which(oxygen21$ScientificName %in%  taxo_tot$Species)
# ox21nospe   <- oxygen21[-idsp21,] #puting on the side the species without a complete species name
# oxygen21    <- oxygen21[idsp21,]
# matchingspeccode21nb     <- which(taxo_tot$Species %in% oxygen21$Species)
# matchingspeccode21       <- taxo_tot[matchingspeccode21nb, c("Species", "SpecCode")]
# oxygen21    <- full_join(oxygen21, matchingspeccode21)
# complete genus and species and speccode columns
oxygen99$Genus           <- str_split_fixed(oxygen99$Species, " ", n=2)[,1]

void16                   <- which(oxygen16$Species == "")
oxygen16$Species[void16] <- paste0(oxygen16$Genus[void16], " sp")
oxygen16$Genus           <- str_split_fixed(oxygen16$Species, " ", n=2)[,1]




##***************
# extract oxygen data from fishbase and clean it
oxygen()      -> oxygenbase

oxygenbase$Weight     <- oxygenbase$Weight*10^-3    # Weight in g
oxygenbase$OxygenCons <- oxygenbase$OxygenCons*oxygenbase$Weight
oxygenbase <- oxygenbase[-which(oxygenbase$Weight == 0),]


##****************
#merge data from Fishbase, Clarke, Ikeda
Ref          <- c(rep("fishbase", dim(oxygenbase)[1]), rep("clarke", dim(oxygen25)[1]), rep("ikeda", dim(oxygen16)[1]), 
                  rep("killen", dim(oxygen16b)[1]), rep("gravel", dim(oxygen24)[1]))
OxygenCons   <- c(oxygenbase$OxygenCons, oxygen25$MR..d.respi, oxygen16$OxygenCons, oxygen16b$RMRsource, oxygen24$OxygenCons)
Weight       <- c(oxygenbase$Weight,     oxygen25$Mb..b.kg,   oxygen16$kgWM,   oxygen16b$RMRmasskg, oxygen24$MaxBodyMass_KG)
SpecCode     <- c(oxygenbase$SpecCode,      oxygen25$SpecCode,      oxygen16$SpecCode, oxygen16b$SpecCode, oxygen24$SpecCode)
Temperature  <- as.numeric(c(oxygenbase$Temp,       oxygen25$Tb..c., oxygen16$Temp, oxygen16b$RMRtemp, oxygen24$MRTempCelcius))
AppliedStress<- c(oxygenbase$AppliedStress, rep("NA", length(oxygen25$SpecCode)), rep("NA", length(oxygen16$SpecCode)), 
                  rep("NA", length(oxygen16b$SpecCode)), rep("NA", length(oxygen24$SpecCode)))
MetabolicLevel<-c(oxygenbase$MetabolicLevel, rep("standard", length(oxygen25$SpecCode)), rep("routine", length(oxygen16$SpecCode)), 
                  rep("standard", length(oxygen16b$SpecCode)), rep("standard", length(oxygen24$SpecCode)))
OxygenRefNo  <- as.numeric(c(oxygenbase$OxygenRefNo, rep("NA", length(oxygen25$SpecCode)), rep("NA", length(oxygen16$SpecCode)), 
                             rep("NA", length(oxygen16b$SpecCode)),  rep("NA", length(oxygen24$SpecCode))))
oxdata       <- data.frame(Ref, OxygenCons, Weight, Temperature, SpecCode, MetabolicLevel, AppliedStress, OxygenRefNo)


##****************
#remove the false replicates
specc_tab <- table(oxdata$SpecCode)
specc_vals <- names(specc_tab[specc_tab > 1])
nb_repetitions <- 0
ref_to_remove <- data.frame(
  SpecCode = character(),
  dataset_1 = character(),
  dataset_2 = character(),
  rowname_1 = character(),
  rowname_2 = character(),
  stringsAsFactors = FALSE
)
seen_keys <- c()
for (specc in specc_vals){
  
  data_specc <- oxdata[which(oxdata$SpecCode == specc),]
  if (length(names(table(data_specc$Ref)))==1){next}
  for (ref in names(table(data_specc$Ref))){
    data_specc_ref <- data_specc[which(data_specc$Ref == ref),]
    data_specc_not_ref <- data_specc[-which(data_specc$Ref == ref),]
    for (obs in 1:nrow(data_specc_ref)){
      ID_closest <- DescTools::Closest(x = data_specc_not_ref$OxygenCons,
                         a = data_specc_ref$OxygenCons[obs], which = T)
      ox_diff <- data_specc_not_ref[ID_closest,"OxygenCons"]-data_specc_ref[obs,"OxygenCons"]
      w_diff <- data_specc_not_ref[ID_closest,"Weight"]-data_specc_ref[obs,"Weight"]
      deg_diff <- data_specc_not_ref[ID_closest,"Temperature"]-data_specc_ref[obs,"Temperature"]
      
      if(abs(ox_diff[1])<1){
        if(abs(w_diff[1])<0.5){
          if(abs(deg_diff[1])<0.51){
            if(any(data_specc_not_ref[ID_closest, "MetabolicLevel"] %in% data_specc_ref[obs,"MetabolicLevel"])){
              nb_repetitions <- nb_repetitions+1
              # cat(unlist(data_specc_not_ref[ID_closest,]), "\nand\n",  unlist(data_specc_ref[obs,]), "\nvalues are really close\n")
              ref_to_remove <- rbind(ref_to_remove, data.frame("SpecCode"= specc, 
                                        "dataset_1"=data_specc_not_ref[ID_closest[1], "Ref"], 
                                        "dataset_2"=data_specc_ref[obs, "Ref"], 
                                        "rowname_1"=rownames(data_specc_not_ref[ID_closest[1],]), 
                                        "rowname_2"=rownames(data_specc_ref[obs,])))
              if(length(ID_closest)==2){
                cat(unlist(data_specc_not_ref[ID_closest,]), "\nand\n",  unlist(data_specc_ref[obs,]), "\nvalues are really close\n")
                ref_to_remove <- rbind(ref_to_remove, data.frame("SpecCode"= specc, 
                                                                 "dataset_1"=data_specc_not_ref[ID_closest[2], "Ref"], 
                                                                 "dataset_2"=data_specc_ref[obs, "Ref"], 
                                                                 "rowname_1"=rownames(data_specc_not_ref[ID_closest[2],]), 
                                                                 "rowname_2"=rownames(data_specc_ref[obs,])))
                rn1 <- rownames(data_specc_not_ref)[ID_closest[2]]
                rn2 <- rownames(data_specc_ref)[obs]
                key <- paste(sort(c(rn1, rn2)), collapse = "_")
                key2 <- paste(sort(c(rn2, rn1)), collapse = "_")
                if (!(key %in% seen_keys) || !(key2 %in% seen_keys)) {
                  seen_keys <- c(seen_keys, key)
                }
              }
              if(length(ID_closest)>2){cat("length ID_closest > 2\n")}
              
              rn1 <- rownames(data_specc_not_ref)[ID_closest[1]]
              rn2 <- rownames(data_specc_ref)[obs]
              key <- paste(sort(c(rn1, rn2)), collapse = "_")
              key2 <- paste(sort(c(rn2, rn1)), collapse = "_")
              if (!(key %in% seen_keys) || !(key2 %in% seen_keys)) {
                seen_keys <- c(seen_keys, key)
              }
            }
          }
        }
      }
    }
  }
}
ID_remove <- unlist(stringr::str_split(seen_keys, pattern = "_"))#[seq(1,length(seen_keys)*2, by=2)]
oxdata    <- oxdata[-which(rownames(oxdata) %in% ID_remove),]


specc_tab <- table(oxdata$SpecCode)
specc_vals <- names(specc_tab[specc_tab > 1])
nb_repetitions <- 0
ref_to_remove <- data.frame(
  SpecCode = character(),
  dataset_1 = character(),
  dataset_2 = character(),
  rowname_1 = character(),
  rowname_2 = character(),
  stringsAsFactors = FALSE
)
seen_keys <- c()
for (specc in specc_vals){
  
  data_specc <- oxdata[which(oxdata$SpecCode == specc),]
  if (length(names(table(data_specc$Ref)))==1){next}
  for (ref in names(table(data_specc$Ref))){
    data_specc_ref <- data_specc[which(data_specc$Ref == ref),]
    data_specc_not_ref <- data_specc[-which(data_specc$Ref == ref),]
    for (obs in 1:nrow(data_specc_ref)){
      ID_closest <- DescTools::Closest(x = data_specc_not_ref$OxygenCons,
                                       a = data_specc_ref$OxygenCons[obs], which = T)
      ox_diff <- (data_specc_not_ref[ID_closest,"OxygenCons"])*10^4-(data_specc_ref[obs,"OxygenCons"])*10^4
      w_diff <- data_specc_not_ref[ID_closest,"Weight"]-data_specc_ref[obs,"Weight"]
      deg_diff <- data_specc_not_ref[ID_closest,"Temperature"]-data_specc_ref[obs,"Temperature"]
      
      if(min(w_diff)==0){
        if(min(deg_diff)==0){
          if(any(data_specc_not_ref[ID_closest, "MetabolicLevel"] %in% data_specc_ref[obs,"MetabolicLevel"])){
            nb_repetitions <- nb_repetitions+1
            cat(unlist(data_specc_not_ref[ID_closest,]), "\nand\n",  unlist(data_specc_ref[obs,]), "\nvalues are really close\n")
            ref_to_remove <- rbind(ref_to_remove, data.frame("SpecCode"= specc, 
                                                             "dataset_1"=data_specc_not_ref[ID_closest[1], "Ref"], 
                                                             "dataset_2"=data_specc_ref[obs, "Ref"], 
                                                             "rowname_1"=rownames(data_specc_not_ref[ID_closest[1],]), 
                                                             "rowname_2"=rownames(data_specc_ref[obs,])))
            if(length(ID_closest)==2){
              cat(unlist(data_specc_not_ref[ID_closest,]), "\nand\n",  unlist(data_specc_ref[obs,]), "\nvalues are really close\n")
              ref_to_remove <- rbind(ref_to_remove, data.frame("SpecCode"= specc, 
                                                               "dataset_1"=data_specc_not_ref[ID_closest[2], "Ref"], 
                                                               "dataset_2"=data_specc_ref[obs, "Ref"], 
                                                               "rowname_1"=rownames(data_specc_not_ref[ID_closest[2],]), 
                                                               "rowname_2"=rownames(data_specc_ref[obs,])))
              # rn1 <- rownames(data_specc_not_ref)[ID_closest[2]]
              # rn2 <- rownames(data_specc_ref)[obs]
              # key <- paste(sort(c(rn1, rn2)), collapse = "_")
              # key2 <- paste(sort(c(rn2, rn1)), collapse = "_")
              # if (!(key %in% seen_keys) || !(key2 %in% seen_keys)) {
              #   seen_keys <- c(seen_keys, key)
              # }
            }
            if(length(ID_closest)>2){cat("length ID_closest > 2\n")}
            
            rn1 <- rownames(data_specc_not_ref)[ID_closest[1]]
            rn2 <- rownames(data_specc_ref)[obs]
            key <- paste(sort(c(rn1, rn2)), collapse = "_")
            key2 <- paste(sort(c(rn2, rn1)), collapse = "_")
            if (!(key %in% seen_keys) || !(key2 %in% seen_keys)) {
              seen_keys <- c(seen_keys, key)
            }
          }
        }
      }
    }
  }
}

ID_remove <- unlist(stringr::str_split(seen_keys, pattern = "_"))

oxdata -> oxdatabis
oxdatabis -> oxdata

for (rm in c(1:length(ID_remove[seq(1,length(seen_keys)*2, by=2)]))[1:7]){
  id_remove <- which(oxdata[which(rownames(oxdata) %in% ID_remove[c(rm, rm*2)]),"Ref"] %in% c("killen","ikeda", "gravel", "clarke"))[1]
  if (is.na(id_remove)){next}
  if(length(which(rownames(oxdata) %in% ID_remove[c(rm, rm*2)[id_remove]]))==0){next}
  oxdata    <- oxdata[-which(rownames(oxdata) %in% ID_remove[c(rm, rm*2)[id_remove]]),]
}




##****************
#organise the data set newly created
oxdata   <- left_join(oxdata, totspe_full[, which(colnames(totspe_full) %in% c("SpecCode", "Brack", "Saltwater", "DemersPelag", colnames(taxo_tot)))], by="SpecCode")
oxdata   <- oxdata[which(oxdata$AppliedStress %in% c("none specified", "NA")),]
oxdata   <- oxdata[which(oxdata$MetabolicLevel %in% c("standard", "routine", "NA") ),]
oxdata   <- left_join(oxdata, taxo_tot)


dataframe_boxplot <- data.frame(as.factor(Ref),Weight, 
                                Temperature, OxygenCons)
colnames(dataframe_boxplot) <- c("Ref", "Weight", "Temp", "OxygenCons")

dataplot <- data.frame("Ref" = as.factor(oxdata$Ref), "habitat" = oxdata$DemersPelag, "Weight^beta" = (oxdata$Weight*10^3)^(3/4),
                       "-1/Temp" = -1/(oxdata$Temperature), "log(OxygenCons)" = log(oxdata$OxygenCons), "Genus" = oxdata$Genus, 
                       "log(OxCons/Weight)" = log(oxdata$OxygenCons)-log((oxdata$Weight*10^3)^(3/4)))


# ##********** Plot scripts if needed
# ggplot(dataframe_boxplot, aes(x=Ref, y=log(Weight))) +
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), 
#                                  c("fishbase", "ikeda"), c("killen", "fishbase"), 
#                                  c("clarke", "killen"), c("killen", "ikeda"),
#                                  c("gravel", "killen"), c("gravel", "ikeda"),
#                                  c("gravel", "clarke"), c("gravel", "fishbase")),
#               map_signif_level=TRUE)
# ggplot(dataframe_boxplot, aes(x=Ref, y=Temp)) +
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), 
#                                  c("fishbase", "ikeda"), c("killen", "fishbase"), 
#                                  c("clarke", "killen"), c("killen", "ikeda"),
#                                  c("gravel", "killen"), c("gravel", "ikeda"),
#                                  c("gravel", "clarke"), c("gravel", "fishbase")),
#               map_signif_level=TRUE)
# ggplot(dataframe_boxplot, aes(x=Ref, y=log(OxygenCons))) +
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), 
#                                  c("fishbase", "ikeda"), c("killen", "fishbase"), 
#                                  c("clarke", "killen"), c("killen", "ikeda"),
#                                  c("gravel", "killen"), c("gravel", "ikeda"),
#                                  c("gravel", "clarke"), c("gravel", "fishbase")),
#               map_signif_level=TRUE)



##*******************
#add measures where species is not specified (none in the data from Clarke 1999)
vecnames<-c()
for (i in seq_len(dim(ox16nospe)[1])){
  fam <- ox16nospe$Genus[i]
  if(sum(oxdata$Family == fam)>1){
    vecnames<-c(vecnames,i)
  }
}
if (length(vecnames)>0){
  c("There are taxa without species names to add, and the code needs to be completed")
  demernospe <- c()
  for (i in seq_len(length(vecnames))){
    fam   <- ox16nospe$Genus[vecnames[i]]
    idfam <- which(oxdata$Family == fam)
    idmax <- max(table(oxdata$DemersPelag[idfam]))
    demernospe[i] <- names(table(oxdata$DemersPelag[idfam]))[idmax] 
  }
}

##*********************
#remove the individuals that make the slope of the plot oxygen~weight or oxygen~temperature being negative according to a sampling effort or habitat
tabsp <- table(dataplot$Genus)
r=0
sp=c()
genusNEG <- c()
for (i in c(1:length(tabsp))){
  ID_sp <- which(dataplot$Genus == names(tabsp)[i])
  data_tempo <- dataplot[ID_sp,]
  if (length(which(is.infinite(data_tempo$X.1.Temp)))>0){data_tempo <- data_tempo[-which(is.infinite(data_tempo$X.1.Temp)),]}
  
  if (length(data_tempo$habitat)<2){next}
  pa <- ggplot(data_tempo, aes(X.1.Temp, log.OxCons.Weight., colour = habitat)) +
    geom_point() +
    geom_smooth(span = 0.8, method='lm')
  pb <- ggplot(data_tempo, aes(Weight.beta, log.OxygenCons., colour = habitat)) +
    geom_point() +
    geom_smooth(span = 0.8, method = 'lm')
  r = r+1
  sp = c(sp, names(tabsp)[i])
  assign(x = paste0("plot", names(tabsp)[i]), value = ggarrange(pa, pb))
  
  if (length(table(data_tempo$habitat))==1){
    fit1 <- lm( log.OxCons.Weight.~X.1.Temp, data = data_tempo)
    fit2 <- lm( log.OxygenCons.~Weight.beta, data = data_tempo)
    ID1 <- 2
    ID2 <- 2
    sp_name <- names(tabsp)[i]
  } else {
    fit1 <- lm( log.OxCons.Weight.~X.1.Temp*habitat, data = data_tempo)
    fit2 <- lm( log.OxygenCons.~Weight.beta*habitat, data = data_tempo)
    ID1 <- c(2, grep(names(coef(fit1)), pattern = ":"))
    ID2 <- c(2, grep(names(coef(fit2)), pattern = ":"))
    # sp_name <- paste0(names(tabsp)[i], "_", 
           # unique(c(names(table(data_tempo$habitat))[coef(fit1)[ID1]<0], names(table(data_tempo$habitat))[coef(fit2)[ID2]<0])))
  }
  
  if (any(is.na(coef(fit1)[ID1]), is.na(coef(fit2)[ID2]))) {
    genusNEG <- c(genusNEG, names(tabsp)[i])
    next}
  if(any(coef(fit1)[ID1]<0) || any(coef(fit2)[ID2]<0)){
    genusNEG <- c(genusNEG, sp_name)}
}

sort_sp <- sort(sp)
pdf(file=paste0(pathoutputplot, "/genus_variables_check.pdf"))
for (k in sort_sp){
  plot <- get(paste0("plot", k))
  print(annotate_figure(plot, top = text_grob(k, face = "bold", size = 14)))
}
dev.off()



oxdata -> oxdatabis
oxdatabis -> oxdata

for (i in genusNEG){
  cat(dim(oxdata), "\n")
  if(length(which((oxdata$Genus == i)))==0){
    cat("no ", i, " in the dataset \n")
    next}
  if (!stringr::str_detect(i, pattern = "_")){
    oxdata <- oxdata[-which(oxdata$Genus == i),]
  }else {
    sp <- stringr::str_split(i, "_")[[1]][1]
    hab <- stringr::str_split(i, "_")[[1]][2]
    if( length(hab)==0){cat("issue with sp ", i)
      next}
    oxdata <- oxdata[-which(oxdata$Genus == sp & oxdata$DemersPelag == hab),]
  }
  if (dim(oxdata)[1] == 0){cat(i)}
}





##*******************
## Clean Oxygen dataset

#get rid of the species that do not have enough individuals
genustemp <- c()
genusnb   <- as.numeric(table(oxdata$Genus))
genusnames<- names(table(oxdata$Genus))
genustorm <- genusnames[which(genusnb<3)]  
if (length(genustorm)>0){oxdata    <- oxdata[-which(oxdata$Genus %in% genustorm),]}
#get rid of the genus that do not have a big enough range of temperature
genustemp <- c()
namestemp<- unlist(lapply(seq_along(table(oxdata$Genus)), function(i){lev <- oxdata$Temperature[which(oxdata$Genus == names(table(oxdata$Genus))[i])]
gap  <- max(lev)-min(lev)
range<- gap>1
return(range)}))
genustemp<- names(table(oxdata$Genus))[!namestemp]
if (length(genustemp)>0){oxdata    <- oxdata[-which(oxdata$Genus %in% genustemp),]}
#get rid of the genus that do not have a big enough range of weights
genustemp2 <- c()
namestemp2<- unlist(lapply(seq_along(table(oxdata$Genus)), function(i){lev <- oxdata$Weight[which(oxdata$Genus == names(table(oxdata$Genus))[i])]
  gap  <- max(lev)-min(lev)
  range<- gap>=(max(lev)*0.05)
  return(range)
  }))
genustemp2<- names(table(oxdata$Genus))[!namestemp2]
if (length(genustemp2)>0){oxdata    <- oxdata[-which(oxdata$Genus %in% genustemp2),]}
dim(oxdata)


##**************************
##    CREATE OUTPUTS   ######
##**************************
### create csv of oxygen data
write.csv(oxdata, paste0(pathoutput, "/dataset_oxygen.csv"))
write.csv(totspe, paste0(pathoutput, "/dataset_totspe.csv"))

