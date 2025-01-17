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


setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/repository_for_zenodo")
path <- paste0(getwd(), "/01-Simulations/Scripts")
pathoutput <- paste0(getwd(), "/01-Simulations/Outputs/dataset_creation_output")


##*************************
## TAXO DATA : organise + fill in ######
##*************************

##
ordersnames <- 
  ##
c("Clupeiformes" , "Carangiformes", "Scorpaeniformes", "Labriformes", "Perciformes", "Gobiiformes" , "Anguilliformes" , "Batrachoidiformes" , 
  "Pleuronectiformes", "Mugiliformes", "Lophiiformes", "Gadiformes", "Atheriniformes", "Scombriformes", "Aulopiformes", 
  "Beloniformes", "Tetraodontiformes", "Carcharhiniformes", "Rajiformes")


## extract taxo and species data
species()  -> totspe
totspe     <- totspe[-which(totspe$SpecCode %in% c(0,1)),]
totspe     <- totspe[which(totspe$Saltwater==1 | totspe$Brack==1),]
totspe     <- totspe[-which(totspe$Fresh==1 & totspe$Brack==1 & totspe$Saltwater==0),]
totspe     <- totspe[, names(totspe) %in% c("SpecCode", "DemersPelag", "Saltwater", "Brack")]
totspe$DemersPelag[totspe$DemersPelag %in% c("pelagic", "pelagic-oceanic", "pelagic-neritic", "bathypelagic")] <- c("pelagic")
totspe$DemersPelag[totspe$DemersPelag %in% c("reef-associated")] <- c("benthopelagic")
totspe$DemersPelag[totspe$DemersPelag %in% c("bathydemersal")] <- c("demersal")
taxo_tot   <- load_taxa()
taxo_tot   <- taxo_tot[-which(taxo_tot$SpecCode %in% c(0,1)),]
spetot     <- full_join(totspe, taxo_tot, by=c("SpecCode"))

IDNA       <- which(is.na(spetot$DemersPelag))
spetot     <- spetot[-IDNA, ]

## fill the family name NAs
fillNAs <- read.csv2(paste0(getwd(), "/01-Simulations/Inputs/completerfishbase.csv"))
famnames <- spetot$Family
for (i in seq_along(famnames)){
  if (is.na(famnames[i]) & !spetot$SpecCode[i]%in%c(1, 0)){
    ID    <- spetot$SpecCode[i]
    locID <- which(fillNAs$code==ID)
    famnames[i] <- fillNAs$fam[locID]
  }
}
spetot$Family <- famnames


## taxo_tot() data : many missing orders that are refereed to as "Eupercaria/misc" (only the ones from osmose)
## here we are replacing the missing data with what it should be according to a litterature research
##******************************************************************************************
famlabr          <- which(spetot$Family %in% c("Labridae", "Odacidae", "Scaridae"))
spetot$Order[famlabr]   <- rep("Labriformes", length(famlabr))
famloph<- which(spetot$Family %in% c("Antennariidae", "Brachionichthyidae", "Caulophrynidae", "Centrophrynidae", "Ceratiidae", 
                                     "Chaunacidae", "Diceratiidae", "Gigantactinidae", "Himantolophidae", "Linophrynidae", 
                                     "Lophichthyidae", "Lophiidae", "Melanocetidae", "Neoceratiidae", "Ogcocephalidae", "Oneirodidae",
                                     "Tetrabrachiidae", "Thaumatichthyidae"))
spetot$Order[famloph]   <- rep("Lophiiformes", length(famloph))
fambatra         <- which(spetot$Family %in% c("Batrachoididae"))
spetot$Order[fambatra]  <- rep("Batrachoidiformes", length(fambatra))
famaulo          <- which(spetot$Family %in% c(" Alepisauridae", "Anotopteridae", "Aulopidae", "Evermannellidae", "Giganturidae", "Ipnopidae",
                                               "Bathysauridae", "Bathysauroididae", "Bathysauropsidae", "Chlorophthalmidae",
                                               "Notosudidae", "Omosudidae", "Paralepididae", "Paraulopidae", "Pseudotrichonotidae", 
                                               "Scopelarchidae", "Synodontidae"))
spetot$Order[famaulo]   <- rep("Aulopiformes", length(famaulo)) 
famraj           <- which(spetot$Family %in% c("Anacanthobatidae", "Arhynchobatidae", "Rajidae", "Rhinobatidae"))
spetot$Order[famraj]    <- rep("Rajiformes", length(famraj)) 
fambelo          <- which(spetot$Family %in% c("Adrianichthyidae", "Belonidae", "Exocoetidae", "Hemiramphidae", 
                                               "Scomberesocidae", "Zenarchopteridae"))
spetot$Order[fambelo]   <- rep("Beloniformes", length(fambelo)) 
famscor          <- which(spetot$Family %in% c("Anoplopomatidae", " Abyssocottidae", " Apistidae",
                                               "Agonidae", "Aploactinidae", "Aploactinidae", "Bembridae", 
                                               "Caracanthidae", "Comephoridae", "Congiopodidae", "Cottidae",
                                               "Cottocomephoridae", "Cyclopteridae", "Dactylopteridae", 
                                               "Ereuniidae", "Eschmeyeridae", "Gnathanacanthidae", "Hemitripteridae",
                                               "Hexagrammidae", "Hoplichthyidae", "Liparidae", "Neosebastidae", 
                                               "Normanichthyidae", "Parabembridae ", "Pataecidae", "Peristediidae ", 
                                               "Platycephalidae", "Plectrogenidae", "Psychrolutidae", "Rhamphocottidae",
                                               "Scorpaenidae", "Sebastidae", "Setarchidae", "Synanceiidae", "Tetrarogidae",
                                               "Triglidae"))
spetot$Order[famscor]   <- rep("Scorpaeniformes", length(famscor)) 
famathe          <- which(spetot$Family %in% c("Atherinidae", "Bedotiidae", "Atherinopsidae", 
                                               "Dentatherinidae", "Isonidae", "Melanotaeniidae", 
                                               "Notocheiridae", "Phallostethidae", "Pseudomugilidae", "Telmatherinidae"))
spetot$Order[famathe]   <- rep("Atheriniformes", length(famathe)) 
ordpercif        <- which(grepl("Percifor", spetot$Order, fixed=TRUE))
spetot$Order[ordpercif] <- rep("Perciformes", length(ordpercif)) 
##*********************************************************************





##****************************
##   OXYGEN DATA  :  merge 3 datasets   ######
##****************************

# * about the units : 
# * 
# * Fishbase : selected routine + standard + NA metabolism
# * Clark :    standard metabolism
# * Ikeda :    routine metabolism
# * 
# * Fishbase : Temperature (°C) : mean water temperature during the experiment
# * Clark :    Temperature (°C) : experimental temperature
# * Ikeda :    Temperature (°C) : in situ temperature
# * 
# * Fishbase : Weight of the test animal (g) or the mean weight if many
# individuals per experiment -> conversion in kg (*10^-3)
# * Clark :    Weight (g) : wet body mass -> conversion in kg *10^-3
# * Ikeda :    Weight (mg) : wet mass is being used here => needs to be
# *10^-6 to be in kg
# * 
# * Fishbase : Oxygen consumption (mg O2 . kg^-1 . h^-1) => needs to be
# * Weight(kg)-> mg O2 . h^-1
# * Clark :    VO2 ? Resting metabolic rate (mmol O2 . h^-1)   => needs to be *
#  32        -> mg O2 . h^-1
# * Ikeda :    Routine respiration (microL O2 . h^-1)  => needs to be converted.
# Mass Volumic : 1.308*10^3 mg/L, so :
# (http://wiki.scienceamusante.net/index.php?title=dioxyg%C3%A8ne)
# *                                        1- microL-> *10^-6 -> L
# *                                        2- L     -> *1.308*10^3-> mg . h^-1


##**********
#Deal with 2 datasets other than Fishbase : 
#extract data from Clarke and Johnston 1999 AND IOkeda et al 2016 and clean it 
oxygen99 <- read.csv2(paste0(getwd(), "/01-Simulations/Inputs/Clarke99_Conso02.csv"))
oxygen16 <- read.csv2(paste0(getwd(), "/01-Simulations/Inputs/Ikeda16_ConsoO2.csv"))[1:102,]

# complete genus and species and speccode columns
oxygen99$Genus           <- str_split_fixed(oxygen99$Species, " ", n=2)[,1]

void16                   <- which(oxygen16$Species == "")
oxygen16$Species[void16] <- paste0(oxygen16$Genus[void16], " sp")
oxygen16$Genus           <- str_split_fixed(oxygen16$Species, " ", n=2)[,1]
taxo_tot                 <- load_taxa()

# add species code to the data + put aside the data that is not complete
idsp16      <- which(oxygen16$Species %in%  taxo_tot$Species)
ox16nospe   <- oxygen16[-idsp16,] #puting on the side the species without a complete species name
oxygen16    <- oxygen16[idsp16,]
matchingspeccode16nb     <- which(taxo_tot$Species %in% oxygen16$Species)
matchingspeccode16       <- taxo_tot[matchingspeccode16nb, c("Species", "SpecCode")]
oxygen16    <- full_join(oxygen16, matchingspeccode16)

idsp99      <- which(oxygen99$Species %in%  taxo_tot$Species)
ox99nospe   <- oxygen99[-idsp99,] #puting on the side the species without a complete species name
oxygen99    <- oxygen99[idsp99,]
matchingspeccode99nb     <- which(taxo_tot$Species %in% oxygen99$Species)
matchingspeccode99       <- taxo_tot[matchingspeccode99nb, c("Species", "SpecCode")]
oxygen99    <- full_join(oxygen99, matchingspeccode99)

# change units according to the units expected on fishbase
oxygen99$MedMb.kg.   <- as.numeric(oxygen99$MedMb.g.)*10^-3
oxygen99$RespiRate   <- as.numeric(oxygen99$VO2.MedMb.)
oxygen99$OxygenCons  <- as.numeric(oxygen99$RespiRate)*32

oxygen16$kgWM        <- as.numeric(oxygen16$mgWM)*(10^-6)
oxygen16$R_L         <- as.numeric(oxygen16$R)*10^-6
oxygen16$OxygenCons  <- oxygen16$R_L*1.308*10^3


##***************
# extract oxygen data from fishbase and clean it
oxygen() -> oxygenbase
oxygenbase$Weight     <- oxygenbase$Weight*10^-3    # Weight in g
oxygenbase$OxygenCons <- oxygenbase$OxygenCons*oxygenbase$Weight
oxygenbase <- oxygenbase[-which(oxygenbase$Weight == 0),]


##****************
#merge data from Fishbase, Clarke, Ikeda
Ref          <- c(rep("fishbase", dim(oxygenbase)[1]), rep("clarke", dim(oxygen99)[1]), rep("ikeda", dim(oxygen16)[1]))
OxygenCons   <- c(oxygenbase$OxygenCons, oxygen99$OxygenCons, oxygen16$OxygenCons)
Weight       <- c(oxygenbase$Weight,     oxygen99$MedMb.kg.,   oxygen16$kgWM)
SpecCode     <- c(oxygenbase$SpecCode,      oxygen99$SpecCode,      oxygen16$SpecCode)
Temperature  <- as.numeric(c(oxygenbase$Temp,       oxygen99$Temp.degC., oxygen16$Temp))
AppliedStress<- c(oxygenbase$AppliedStress, rep("NA", length(oxygen99$SpecCode)), rep("NA", length(oxygen16$SpecCode)))
MetabolicLevel<-c(oxygenbase$MetabolicLevel, rep("NA", length(oxygen99$SpecCode)), rep("NA", length(oxygen16$SpecCode)))
OxygenRefNo  <- as.numeric(c(oxygenbase$OxygenRefNo, rep("NA", length(oxygen99$SpecCode)), rep("NA", length(oxygen16$SpecCode))))
oxdata       <- data.frame(Ref, OxygenCons, Weight, Temperature, SpecCode, MetabolicLevel, AppliedStress, OxygenRefNo)

##****************
#organise the data set newly created
oxdata   <- full_join(oxdata, taxo_tot, by="SpecCode")
oxdata   <- full_join(oxdata, spetot[, which(colnames(spetot) %in% c("SpecCode", "Brack", "Saltwater", "DemersPelag"))], by="SpecCode")
oxdata   <- oxdata[which(oxdata$AppliedStress %in% c("none specified", "NA")),]
oxdata   <- oxdata[which(oxdata$MetabolicLevel %in% c("standard", "routine", "NA") ),]
oxdata   <- oxdata[which(oxdata$Saltwater == 1 | oxdata$Brack == 1), ]


dataframe_boxplot <- data.frame(as.factor(Ref),Weight, 
                                Temperature, OxygenCons)
colnames(dataframe_boxplot) <- c("Ref", "Weight", "Temp", "OxygenCons")


##********** Plot scripts if needed
# ggplot(dataframe_boxplot, aes(x=Ref, y=log(Weight))) + 
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), c("fishbase", "ikeda")), 
#               map_signif_level=TRUE)
# ggplot(dataframe_boxplot, aes(x=Ref, y=Temp)) + 
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), c("fishbase", "ikeda")), 
#               map_signif_level=TRUE)
# ggplot(dataframe_boxplot, aes(x=Ref, y=log(OxygenCons))) + 
#   geom_boxplot()+
#   geom_signif(comparisons = list(c("clarke", "fishbase"), c("clarke", "ikeda"), c("fishbase", "ikeda")), 
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
genusNEG <- c("Anguilla", "Anoplogaster", "Artediellus", "Borostomias", "Callionymus", "Chromis", "Conger", "Coryphaena", "Cyprinus",
              "Diaphus", "Fundulus", "Girella", "Hippocampus", "Hypanus", "Labeo", "Lagodon", "Lepomis", "Lutjanus", "Lycodes",
              "Microstomus", "Notothenia", "Ophiodon", "Orthodon", "Perca", "Rutilus", "Salmo", "Salvelinus", "Sander", "Scophthalmus", 
              "Scorpaena", "Scyliorhinus", "Solea", "Stomias", "Tarletonbeania", "Thunnus", "Torpedo", "Xiphophorus")

oxdata -> oxdatabis
oxdatabis -> oxdata
dimox <-c()
for (i in seq_len(length(genusNEG))){   
  idgenusNEG <- which(oxdata$Genus == genusNEG[i])
  if (length(idgenusNEG)>0){
    if (genusNEG[i] %in% c("Anoplogaster", "Artediellus", "Borostomias", "Callionymus", "Chromis", "Conger", "Coryphaena", "Cyprinus",
              "Diaphus", "Girella", "Hippocampus", "Hypanus", "Labeo", "Lagodon", "Lepomis", "Lutjanus", "Lycodes",
              "Microstomus", "Ophiodon", "Orthodon", "Perca", "Rutilus", "Salmo", "Sander", "Scophthalmus", 
              "Scorpaena", "Solea", "Stomias", "Tarletonbeania", "Thunnus", "Torpedo", "Xiphophorus")){
    oxdata <- oxdata[-idgenusNEG, ]
    }
    if (genusNEG[i] %in% c("Anguilla", "Notothenia", "Salmo", "Scyliorhinus")){
      if (genusNEG[i] %in% c("Anguilla", "Salmo", "Scyliorhinus")){
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
    if (genusNEG[i] %in% c("Fundulus", "Salvelinus", "Scyliorhinus")){
      if (genusNEG[i] %in% c("Salvelinus", "Fundulus")){
        idgenusNEG <- which(oxdata$Genus == genusNEG[i])
        id2120     <- which(oxdata$OxygenRefNo == c(2120))
        idRM          <- idgenusNEG[idgenusNEG %in% id2120]
        oxdata        <- oxdata[-idRM, ]
      }
      if (genusNEG[i] %in% c("Fundulus", "Scyliorhinus")){
        idgenusNEG <- which(oxdata$Genus == genusNEG[i])
        id2203 <- which(oxdata$OxygenRefNo == c(2203))
        idRM          <- idgenusNEG[idgenusNEG %in% id2203]
        oxdata        <- oxdata[-idRM, ]
      }
    }
  dimox <- c(dimox, dim(oxdata)[1])
  }
  else{
    oxdata <- oxdata
    dimox <- c(dimox, 0)}
}


##*******************
##Clean Oxygen dataset
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



##**************************
##    CREATE OUTPUTS   ######
##**************************
### create csv of oxygen data
write.csv(oxdata, paste0(pathoutput, "/dataset_oxygen.csv"))
write.csv(spetot, paste0(pathoutput, "/dataset_spetot.csv"))

