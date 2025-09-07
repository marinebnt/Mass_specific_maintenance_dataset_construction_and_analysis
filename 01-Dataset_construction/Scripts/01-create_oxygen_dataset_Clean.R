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
# * Fishbase : selected routine + standard + NA metabolism
# * Clark :    resting metabolism
# * Ikeda :    routine metabolism
# * Killen :   resting metabolism
# * 
# * Fishbase : Temperature (°C) : mean water temperature during the experiment
# * Clark :    Temperature (°C) : experimental temperature
# * Ikeda :    Temperature (°C) : in situ temperature
# * Killen :   Temperature (°C) : experimental temperature
# * 
# * Fishbase : Weight of the test animal (g) or the mean weight if many
# individuals per experiment -> conversion in kg (*10^-3)
# * Clark :    Weight (g) : wet body mass -> conversion in kg *10^-3
# * Ikeda :    Weight (mg) : wet mass is being used here => needs to be
# *10^-6 to be in kg
# * Killen :   Weight (g) : wet body mass -> conversion in kg (*10^-3)
# * 
# * Fishbase : Oxygen consumption (mg O2 . kg^-1 . h^-1) => needs to be
# * Weight(kg)-> mg O2 . h^-1
# * Clark :    Watt =>  Oxygen consumption mg O2 . h^-1
# * Ikeda :    Routine respiration (microL O2 . h^-1)  => needs to be converted.
# Mass Volumic : 1.308*10^3 mg/L, so :
# (http://wiki.scienceamusante.net/index.php?title=dioxyg%C3%A8ne)
# *                                        1- microL-> *10^-6 -> L
# *                                        2- L     -> *1.308*10^3-> mg O2. h^-1
# * Killen :   mg O2/h adjusted to 1kg and 15°C

source("01-Dataset_construction/Scripts/functions_data_collection.R")
##**********
#Deal with 2 datasets other than Fishbase : 
#extract data from Clarke and Johnston 1999 AND IOkeda et al 2016 and clean it 
oxygen99 <- read.csv(paste0(getwd(), "/01-Dataset_construction/Inputs/Clarke & Johnston 2025/metabolic_rates.csv"), skip = c(25))
oxygen99 <- oxygen99[which(oxygen99$Class == "Osteichthyes"),]
oxygen16 <- read.csv2(paste0(getwd(), "/01-Dataset_construction/Inputs/Ikeda16_ConsoO2.csv"))[1:102,]
oxygen16b <- read.csv(paste0(getwd(), "/01-Dataset_construction/Inputs/Killenetal16Ecological.csv"))
oxygen16b <- oxygen16b[-which(is.na(oxygen16b$RMRmass)),]
  
# change units according to the units expected on fishbase
# clarke and johnston : first all in Watt
oxygen99$MR..d.[which(oxygen99$Unit..MR...e. == "mW")]      <-  
  oxygen99$MR..d.[which(oxygen99$Unit..MR...e. == "mW")]*10^(-3)
oxygen99$Unit..MR...e.[which(oxygen99$Unit..MR...e. == "mW")]      <-  "W"
oxygen99$MR..d.respi <- convert_watt_to_respi(oxygen99$MR..d.)
oxygen99$Mb..b.kg    <- as.numeric(oxygen99$Mb..b.)*10^-3

oxygen16$kgWM        <- as.numeric(oxygen16$mgWM)*(10^-6)
oxygen16$R_L         <- as.numeric(oxygen16$R)*10^-6
oxygen16$OxygenCons  <- oxygen16$R_L*1.308*10^3

oxygen16b$RMRmasskg  <- oxygen16b$RMRmass*10^-3
oxygen16b$Species    <- stringr::str_replace(oxygen16b$species, "_", " ")

# add species code to the data + put aside the data that is not complete
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

oxygen99    <- ox99_correction(oxygen99) # this way are missing only misspecified only the species for which : no species specified
idsp99      <- which(oxygen99$Species %in%  taxo_tot$Species)
ox99nospe   <- oxygen99[-idsp99,] #puting on the side the species without a complete species name
oxygen99    <- oxygen99[idsp99,]
matchingspeccode99nb     <- which(taxo_tot$Species %in% oxygen99$Species)
matchingspeccode99       <- taxo_tot[matchingspeccode99nb, c("Species", "SpecCode")]
oxygen99    <- full_join(oxygen99, matchingspeccode99)

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
Ref          <- c(rep("fishbase", dim(oxygenbase)[1]), rep("clarke", dim(oxygen99)[1]), rep("ikeda", dim(oxygen16)[1]), rep("killen", dim(oxygen16b)[1]))
OxygenCons   <- c(oxygenbase$OxygenCons, oxygen99$MR..d.respi, oxygen16$OxygenCons, oxygen16b$RMRsource)
Weight       <- c(oxygenbase$Weight,     oxygen99$Mb..b.kg,   oxygen16$kgWM,   oxygen16b$RMRmasskg)
SpecCode     <- c(oxygenbase$SpecCode,      oxygen99$SpecCode,      oxygen16$SpecCode, oxygen16b$SpecCode)
Temperature  <- as.numeric(c(oxygenbase$Temp,       oxygen99$Tb..c., oxygen16$Temp, oxygen16b$RMRtemp))
AppliedStress<- c(oxygenbase$AppliedStress, rep("NA", length(oxygen99$SpecCode)), rep("NA", length(oxygen16$SpecCode)), rep("NA", length(oxygen16b$SpecCode)))
MetabolicLevel<-c(oxygenbase$MetabolicLevel, rep("standard", length(oxygen99$SpecCode)), rep("routine", length(oxygen16$SpecCode)), rep("standard", length(oxygen16b$SpecCode)))
OxygenRefNo  <- as.numeric(c(oxygenbase$OxygenRefNo, rep("NA", length(oxygen99$SpecCode)), rep("NA", length(oxygen16$SpecCode)), rep("NA", length(oxygen16b$SpecCode))))
oxdata       <- data.frame(Ref, OxygenCons, Weight, Temperature, SpecCode, MetabolicLevel, AppliedStress, OxygenRefNo)

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


##********** Plot scripts if needed
ggplot(dataframe_boxplot, aes(x=Ref, y=log(Weight))) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), c("fishbase", "ikeda"), c("killen", "fishbase"), c("clarke", "killen"), c("killen", "ikeda")),
              map_signif_level=TRUE)
ggplot(dataframe_boxplot, aes(x=Ref, y=Temp)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), c("fishbase", "ikeda"), c("killen", "fishbase"), c("clarke", "killen"), c("killen", "ikeda")),
              map_signif_level=TRUE)
ggplot(dataframe_boxplot, aes(x=Ref, y=log(OxygenCons))) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), c("fishbase", "ikeda"), c("killen", "fishbase"), c("clarke", "killen"), c("killen", "ikeda")),
              map_signif_level=TRUE)


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
# 1- plot it
tabsp <- table(dataplot$Genus)
r=0
sp=c()
for (i in 1:length(tabsp)){
  ID_sp <- which(dataplot$Genus == names(tabsp)[i])
  if (length(ID_sp)<2){next}
    pa <- ggplot(dataplot[ID_sp,], aes(X.1.Temp, log.OxCons.Weight., colour = habitat)) +
    geom_point() +
    geom_smooth(span = 0.8, method='lm')
    pb <- ggplot(dataplot[ID_sp,], aes(Weight.beta, log.OxygenCons., colour = habitat)) +
    geom_point() +
    geom_smooth(span = 0.8, method = 'lm')
    r = r+1
    sp = c(sp, names(tabsp)[i])
    assign(x = paste0("plot", names(tabsp)[i]), value = ggarrange(pa, pb))
}
sort_sp <- sort(sp)
pdf(file=paste0(pathoutputplot, "/genus_variables_check.pdf"))
for (k in sort_sp){
  
  plot <- get(paste0("plot", k))
  print(annotate_figure(plot, top = text_grob(k, face = "bold", size = 14)))
}
dev.off()


# 2- from the shape of the plots : manual selection of the species to exclude
# here : we haven't applied the species selection yet, taht's why you might find species that won't stay in the final dataset, such as Fundulus genus

genusNEG <- c("Anguilla", "Borostomias", "Conger", "Diplodus", "Echiichthys", 
              "Fundulus", "Gymnocephalus","Girella", "Hypanus",  "Labeo", "Lagodon", 
               "Lampanyctus", "Microstomus", "Myoxocephalus", "Orthodon", "Pimephales", "Planiliza", 
              "Scophthalmus",  "Scyliorhinus", "Stomias", "Syngnathus",
               "Tautogolabrus", "Torpedo", "Xiphophorus", "Zoarces")
oxdata -> oxdatabis
oxdatabis -> oxdata
dimox <-c()
for (i in seq_len(length(genusNEG))){   
  cat(i, "\t", dim(oxdata), "\n") # if 0 appear in the sequence, it is not normal
  idgenusNEG <- which(oxdata$Genus == genusNEG[i])
  if (length(idgenusNEG)>0){
    if (genusNEG[i] %in% c("Borostomias", "Conger", "Diplodus", "Echiichthys", 
                           "Fundulus","Gymnocephalus","Girella", "Hypanus", "Lagodon", 
                           "Lampanyctus", "Microstomus", "Myoxocephalus", "Orthodon",  "Pimephales", "Planiliza", 
                           "Scophthalmus",  "Stomias", "Syngnathus",
                           "Tautogolabrus", "Xiphophorus", "Zoarces")){
    oxdata <- oxdata[-idgenusNEG, ]
    }
    if (genusNEG[i] %in% c("Anguilla", "Scyliorhinus", 
                           "Torpedo", "Labeo")){
      if (genusNEG[i] %in% c("Scyliorhinus", "Anguilla")){
        idbenthopelag <- which(oxdata$DemersPelag == c("benthopelagic"))
        idRM          <- idgenusNEG[idgenusNEG %in% idbenthopelag]
        oxdata        <- oxdata[-idRM, ]
      }      
      else {
        iddemer <- which(oxdata$DemersPelag == c("demersal")) 
        idRM          <- idgenusNEG[idgenusNEG %in% iddemer]
        oxdata        <- oxdata[-idRM, ]
      }
    }
  dimox <- c(dimox, dim(oxdata)[1])
  }
  if(length(idgenusNEG)==0){
    oxdata <- oxdata
    dimox <- c(dimox, 0)}
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

