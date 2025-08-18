######### prepare dataset for phylosem usage : with all data from fishbase ################
# except : fish from not salted water
# needed  : output of the create_dataset_fishbase_c_m_eps_m.R file, called :
# spetot_fishabse_c_m_eps_m.csv
####################################################################


#  INITIALIZATION
##############################################################################################
library(rfishbase)
library(stringr)
library(dplyr)
library(ggplot2)
library(grDevices)
library(ggpubr)

#############
path <- paste0("01-Dataset_construction/Scripts")
pathoutput <- paste0("01-Dataset_construction/Outputs/dataset_creation_output")
pathinput <- paste0("01-Dataset_construction/Inputs")
#############
spetot <- read.csv(paste0(pathoutput, "/spetot_fishabse_c_m_eps_mLMnoHABTREF.csv"))
spetot <- spetot[,-which(colnames(spetot) %in% c("X.1", "X"))]


#my species 
osmosespnames <- c("Alosa alosa", "Alosa fallax", "Anguilla anguilla",
"Argyrosomus regius", "Aristaeomorpha foliacea", "Aristeus antennatus",
"Atherina boyeri", "Auxis rochei", "Belone belone", "Boops boops",
"Caranx crysos", "Chelidonichthys lucerna",
"Coris julis", "Coryphaena hippurus", "Crangon crangon",
"Crystallogobius linearis", "Dentex dentex", "Dentex gibbosus",
"Dentex maroccanus", "Dicentrarchus labrax", "Diplodus annularis",
"Diplodus cervinus", "Diplodus puntazzo", "Diplodus sargus",
"Diplodus vulgaris",
"Eledone cirrhosa", "Engraulis encrasicolus", "Epinephelus aeneus",
"Epinephelus marginatus", "Etrumeus sadina", "Eutrigla gurnardus",
"Galeus melastomus", "Gobius niger", "Halobatrachus didactylus",
"Illex coindetii", "Lepidorhombus whiffiagonis", "Chelon auratus",
"Chelon ramada", "Chelon saliens", "Loligo vulgaris", "Lophius budegassa",
"Lophius piscatorius", "Merlangius merlangus", "Merluccius merluccius",
"Micromesistius poutassou", "Mugil cephalus", "Mullus barbatus",
"Mullus surmuletus", "Mustelus mustelus", "Nephrops norvegicus",
"Octopus vulgaris", "Pagellus acarne", "Pagellus erythrinus", "Pagrus pagrus",
"Palaemon serratus", "Palinurus elephas", "Parapenaeus longirostris",
"Penaeus kerathurus", "Phycis phycis", "Platichthys flesus",
"Pleuronectes platessa", "Pomatomus saltatrix", "Pomatoschistus marmoratus",
"Pomatoschistus minutus", "Rhinobatos rhinobatos",
"Sarda sarda", "Sardina pilchardus", "Sardinella aurita",
"Saurida undosquamis", "Sciaena umbra",
"Scomber japonicus", "Scomber scombrus", "Scophthalmus maximus",
"Scorpaena notata", "Scyliorhinus canicula",
"Sepia officinalis", "Seriola dumerili", "Serranus atricauda",
"Solea solea", "Sparus aurata",
"Sphyraena sphyraena", "Sphyraena viridensis", "Spicara maena",
"Spicara smaris", "Spondyliosoma cantharus",
"Sprattus sprattus", "Squilla mantis", "Stephanolepis diaspros",
"Thunnus alalunga", "Thunnus thynnus",
"Trachurus mediterraneus", "Trachurus picturatus", "Trachurus trachurus",
"Trachyrincus scabrus", "Trigla lyra",
"Trisopterus luscus", "Trisopterus minutus", "Upeneus moluccensis",
"Xiphias gladius", "Gobius ophiocephalus", "Ost euphausiids")
##############################################################################################

# Now we have the oxygen consumption parameters for some genus, we can infer the missing parameters for the remaining ones
########
#* complementary data and dataset creation for phylosem
#* https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1095-8649.2000.tb00870.x to choose the length indicator
#* https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.1890/0012-9658(2001)082[0784%3AAASAFR]2.0.CO%3B2 This article says :
#* Reminder :  Tr is  the  age  when  50%  of  individuals attain  first  reproduction, K is the von  Bertalanffy growth coefficient,
#* Lr is the  body  length  at  which  50%  of  individuals  attain  first  reproduction,  and Linf is  the  asymptotic  body  length
#* when  age  is infinit. 
#* Also : Linf correlated at 0.94 to Lr
#* https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.1002/eap.1606 p.8 => Linf is correlated at 0.5 to Lmat
#* Check the following writer : pecuchet matthew mcclean ester beukhof 

#* lifespan data
#* according to Pecuchet 2017 : https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12587
#* "Life span is defined as the theoretical maximum expected age for a species and is estimated
#* within FishBase using the growth (K) and length at infinity (Linf) parameters from the 
#* Von Bertalanffy growth equation"
#* On fishbase : https://www.fishbase.se/manual/key%20facts.htm
#* "Life span is the approximate maximum age (tmax) that fish of a given population would reach. 
#* Following Taylor (1958), it is calculated as the age at 95% of Linf, using the parameters of 
#* the von Bertalanffy growth function as estimated above, viz.: tmax = t0 + 3 / K. "

#* What we already have : 
#* spetot contains oxygen data + water column position data + phylogeny +c_m + eps_m

# data from fishlife datasets : the optimal temperature
# Thorson 2023 & Phylosem. To access its data : 
#########







#     EXTRACT DATA
##############################################################################################

#extract data from Price et al 2022
############
# morphometric data
data_morpho <- read.csv2(paste0(pathinput, "/Morphometrics_Price.csv"))

# Cleaning Price 2022 data
##############
# Fit a species code to the morphometrics species
taxo_tot    <- load_taxa()
data_morpho$Species <- str_replace(data_morpho$Tree_name, "_", " ")
idsp        <- which(data_morpho$Species %in%  taxo_tot$Species)
datanospe   <- data_morpho[-idsp,] 
#puting on the side the species without a name in fishbase
dataspe             <- data_morpho[idsp,]
matchingspeccodenb  <- which(taxo_tot$Species %in% dataspe$Species)
matchingspeccode    <- taxo_tot[matchingspeccodenb, c("Species", "SpecCode")]
datawithcode        <- full_join(dataspe, matchingspeccode, by="Species")
morpho_price        <- datawithcode[, c("Standard_length", "Min_caudalpeduncle_depth", "SpecCode", 
                                 "Lower_jaw_length", "Max_body_depth", "Max_fish_width")]
colnames(morpho_price)<- c("Standard_length", "Min_caudalpeduncle_depth", "SpecCode", 
                            "Lower_jaw_length", "Max_body_depth", "Max_body_width")
#Scale morpho data to the geometric mean, following Mosimann 1970 and Price 2019
# geometric mean
data_for_geometric_mean    <- as.data.frame(sapply(morpho_price[, c("Standard_length", "Max_body_depth", "Max_body_width")], as.numeric)) 
morpho_price$geomtric_mean <- (data_for_geometric_mean$Standard_length*data_for_geometric_mean$Max_body_depth*data_for_geometric_mean$Max_body_width)^{1/3}
# Morphometric data scaled by the geometric data
morpho_price[, c("Standard_length", "Min_caudalpeduncle_depth", 
                 "Lower_jaw_length", "Max_body_depth", "Max_body_width")] <- 
  as.data.frame(sapply(morpho_price[, c("Standard_length", "Min_caudalpeduncle_depth", 
                                        "Lower_jaw_length", "Max_body_depth", "Max_body_width")], as.numeric))  / morpho_price$geomtric_mean

# extract data from fishbase
############
# estimates from Fishbase
estimate    <- estimate()
estimateTL  <- estimate[, c("SpecCode", "Troph")]
estimateTem <- estimate[, c("SpecCode", "TempPrefMean")] 
  # estimate of the temperature : https://aquamaps.org/main/
  #FB_Book_MarineAquaMaps_062023.pdf#page=1
estimateLoo <- estimate[, c("SpecCode", "Linf")]
# ecology data
ecology     <- ecology()
ecology     <- ecology[, c("SpecCode", "DietTroph", "FoodTroph")]
# growth dta
growth      <- popgrowth()
growthLm    <- growth[, c("SpecCode", "Lm", "TypeLm")] 
  # growthLm    <- growthLm[-which(growth$Sex %in% c("unsexed", "Unsexed")),]
NA          -> growthLm$Lm[-which(growthLm$TypeLm %in% "TL")]
growthLoo   <- growth[, c("SpecCode", "Loo", "Type")]  
  # growthLoo   <- growthLoo[-which(growth$Sex %in% c("unsexed", "Unsexed")),]
NA          -> growthLoo$Loo[-which(growthLoo$Type %in% "TL")]
growthother <- growth[, c("SpecCode", "K", "tmax", "Winfinity", "M")]   
  # growthother <- growthother[-which(growth$Sex %in% c("unsexed", "Unsexed")),]
  # repro      <- reproduction()
  # repro      <- repro[, c("SpecCode", "RepGuild1")]
#fecundity data
fecun      <- fecundity()
fecun      <- fecun[, c("SpecCode", "FecundityMean", "FecundityMin", "FecundityMax")]
#offspring data (not used in the end)
offsp      <- spawning()
offsp      <- offsp[, c("SpecCode", "LengthOffspringMin")]
#swim data (not used in the end)
swim       <- swimming()
swim       <- swim[, c("SpecCode", "AdultMode", "AspectRatio")]
# morphology data
morpho     <- morphology()
morpho     <- morpho[, names(morpho) %in% c("SpecCode", "PeduncleDepth")]
# morphometric data
# mmtrcs     <- morphometrics()
# mmtrcs     <- mmtrcs[, names(mmtrcs) %in% c("SpecCode", "BD")]
# maturity data
matur      <- maturity()
matur      <- matur[, c("SpecCode", "tm", "AgeMatMin", "AgeMatMin2")]
################


# Cleaning Fishbase data
################
c("balistiform")     -> swim$AdultMode[which(swim$AdultMode %in% c("Balistiform"))]
c("anguilliform")    -> swim$AdultMode[which(swim$AdultMode %in% c("Anguilliform"))]
c("subcarangiform")  -> swim$AdultMode[which(swim$AdultMode %in% c("Subcarangiform"))]
c("tetraodontiform") -> swim$AdultMode[which(swim$AdultMode %in% c("Tetraodontiform"))]
spetot[, "Genus"]    -> gen
spetot$SpecCode[which(gen %in% "Bathypterois")] -> code
c("subcarangiform")  -> swim$AdultMode[which(swim$SpecCode %in% code)]
##############

##############################################################################################





#     PREPARE DATA FOR DATASET
##############################################################################################

# Because there are replicates, data should be simplified :
# creating new data frames with the TL, morpho, climate and ecology data
##############
#_____function_______
removenumericreplicates <- function(speccode, data, minmaxmean){
  min = (minmaxmean == "min")
  max = (minmaxmean == "max")
  mean =(minmaxmean == "mean")
  j=0
  vec    <- c()
  for (i in as.numeric(levels(as.factor(speccode)))) { 
    j=j+1
    ID     <- which(speccode==i)
    datana <- na.omit(data[ID])
    if (length(datana)>0) {
      if (min) {vec[j] <- min(datana)}
      if (max) {vec[j] <- max(datana)}
      if (mean){vec[j] <- mean(datana)}
    }
    else {vec[j] <- NA}
  }
  frame       <- data.frame(as.numeric(levels(as.factor(speccode))), vec)
  return(frame)
}
############## Fishbase data
#____ growth data
lengthdata <- removenumericreplicates(growthLoo$SpecCode, growthLoo$Loo, "mean")
names(lengthdata) <- c("SpecCode", "Loo")
looestim   <- removenumericreplicates(estimateLoo$SpecCode, estimateLoo$Linf, "mean")
names(looestim)   <- c("SpecCode", "Looesti")
lengthdata <- full_join(lengthdata, looestim)
lengthdata$Loo[which(is.na(lengthdata$Loo))] <- lengthdata$Looesti[which(is.na(lengthdata$Loo))]
lengthdata <- lengthdata[,  c("SpecCode", "Loo")]

Lmdata     <- removenumericreplicates(growthLm$SpecCode, growthLm$Lm, "mean")
colnames(Lmdata) <- c("SpecCode", "Lm")
tmaxdata   <- removenumericreplicates(growthother$SpecCode, growthother$tmax, "mean")
colnames(tmaxdata) <- c("SpecCode", "tmax")
Kdata      <- removenumericreplicates(growthother$SpecCode, growthother$K, "mean")
colnames(Kdata) <- c("SpecCode", "K")
Woodata    <- removenumericreplicates(growthother$SpecCode, growthother$Winfinity, "mean")
colnames(Woodata) <- c("SpecCode", "Woo")
Mdata      <- removenumericreplicates(growthother$SpecCode, growthother$M, "mean")
colnames(Mdata) <- c("SpecCode", "M")

framelength        <- full_join(lengthdata, Lmdata)
framelength        <- full_join(framelength, tmaxdata) 
framelength        <- full_join(framelength, Kdata) 
framelength        <- full_join(framelength, Woodata)
framelength        <- full_join(framelength, Mdata)
#____ TL data
frameTLdiet <- removenumericreplicates(ecology$SpecCode, ecology$DietTroph, "mean")
frameTLfood <- removenumericreplicates(ecology$SpecCode, ecology$FoodTroph, "mean")
frameTLesti <- removenumericreplicates(estimate$SpecCode, estimate$Troph, "mean")
names(frameTLesti) <- c("SpecCode", "TLesti")
names(frameTLfood) <- c("SpecCode", "TLFood")
names(frameTLdiet) <- c("SpecCode", "TLDiet")
SpecCode <- as.data.frame(unique(c(frameTLdiet$SpecCode, frameTLesti$SpecCode, frameTLfood$SpecCode)))
colnames(SpecCode) <- c("SpecCode")
frameTLdiet <- full_join(SpecCode, frameTLdiet, by="SpecCode")
frameTL     <- full_join(frameTLdiet, frameTLfood, by="SpecCode")
frameTL     <- full_join(frameTL, frameTLesti, by="SpecCode")
frameTL$TLDiet[which(is.na(frameTL$TLDiet))] <- frameTL$TLFood[which(is.na(frameTL$TLDiet))]
frameTL$TLDiet[which(is.na(frameTL$TLDiet))] <- frameTL$TLesti[which(is.na(frameTL$TLDiet))]
#______Fecundity data
framefecmeannew <- removenumericreplicates(c(fecun$SpecCode, fecun$SpecCode), c(fecun$FecundityMin, fecun$FecundityMax), "mean")
framefecmean    <- removenumericreplicates(fecun$SpecCode, fecun$FecundityMean, "mean")
names(framefecmeannew)  <- c("SpecCode", "fecunditmeannew") 
names(framefecmean) <- c("SpecCode", "fecundity")
SpecCode <- as.data.frame(unique(c(framefecmeannew$SpecCode, framefecmeannew$SpecCode)))
colnames(SpecCode) <- c("SpecCode")
framefecmean <- full_join(SpecCode, framefecmean, by="SpecCode")
framefec     <- full_join(framefecmean, framefecmeannew, by="SpecCode")
framefec$fecundity[which(is.na(framefec$fecundity))] <- framefec$fecunditmeannew[which(is.na(framefec$fecundity))]
#______Offspring data
frameoff <- removenumericreplicates(offsp$SpecCode, offsp$LengthOffspringMin, "mean")
names(frameoff)     <- c("SpecCode", "lengthOffspring") 
#______Maturity data
framematur     <- removenumericreplicates(matur$SpecCode, matur$tm, "mean")
names(framematur)     <- c("SpecCode", "tm") 
frameminmatur  <- removenumericreplicates(matur$SpecCode, matur$AgeMatMin, "mean")
names(frameminmatur)     <- c("SpecCode", "tmmin") 
frameminmatur2 <- removenumericreplicates(matur$SpecCode, matur$AgeMatMin2, "mean")
names(frameminmatur2)     <- c("SpecCode", "tmmin2") 
SpecCode <- as.data.frame(unique(c(framematur$SpecCode, frameminmatur$SpecCode, frameminmatur2$SpecCode)))
colnames(SpecCode) <- c("SpecCode")
framematur <- full_join(SpecCode, framematur, by="SpecCode")
framemat   <- full_join(framematur, frameminmatur, by="SpecCode")
framemat   <- full_join(framemat, frameminmatur2, by="SpecCode")
framemat$tm[which(is.na(framemat$tm))] <- framemat$tmmin[which(is.na(framemat$tm))]
framemat$tm[which(is.na(framemat$tm))] <- framemat$tmmin2[which(is.na(framemat$tm))]
#______Aspect Ratio data
framear <- removenumericreplicates(swim$SpecCode, swim$AspectRatio, "mean")
names(framear)     <- c("SpecCode", "AspectRatio")
#_______Temperature data
frameTemp <- removenumericreplicates(estimateTem$SpecCode, estimateTem$TempPrefMean, "mean")
names(frameTemp)     <- c("SpecCode", "Temperature") 
frameTemp$Temperature[which(frameTemp$Temperature==0)] <- 0.1

############## Price et al data
# #_______Morpho from Price publication : remove replicates
peduncle_depth <- removenumericreplicates(morpho_price$SpecCode, as.numeric(morpho_price$Min_caudalpeduncle_depth), "mean")
colnames(peduncle_depth) <- c("SpecCode", "Min_caudalpeduncle_depth")
Lower_jaw_length         <- removenumericreplicates(morpho_price$SpecCode, as.numeric(morpho_price$Lower_jaw_length), "mean")
colnames(Lower_jaw_length) <- c("SpecCode", "Lower_jaw_length")
Max_body_depth <- removenumericreplicates(morpho_price$SpecCode, as.numeric(morpho_price$Max_body_depth), "mean")
colnames(Max_body_depth) <- c("SpecCode", "Max_body_depth")
Max_body_width           <- removenumericreplicates(morpho_price$SpecCode, as.numeric(morpho_price$Max_body_width), "mean")
colnames(Max_body_width) <- c("SpecCode", "Max_body_width")
merge_morpho <- full_join(peduncle_depth, Lower_jaw_length)
merge_morpho <- full_join(merge_morpho, Max_body_depth)
morpho_tot   <- full_join(merge_morpho, Max_body_width)

# length(which(morpho_tot$SpecCode %in% taxo_tot$SpecCode[which(taxo_tot$Class == "Elasmobranchii")]))
# length(which(morpho$SpecCode %in% taxo_tot$SpecCode[which(taxo_tot$Class == "Elasmobranchii")]))
# length(which(peduncle_depth_price$SpecCode %in% taxo_tot$SpecCode[which(taxo_tot$Class == "Elasmobranchii")]))
# #############


# Rule : Linf/Lmax in [0.66777, 1.3333], otherwise remove Linf
##############
framelength[,c("Loo", "SpecCode")] -> Loo 
length_length() -> Lmax
sum(!is.na(Loo$Loo))
Lmax[, c("LengthMax", "SpecCode")] -> Lmax
Lmax <- removenumericreplicates(Lmax$SpecCode, Lmax$LengthMax, "mean")
names(Lmax)<-c("SpecCode", "Lmax")
lengthratio <- full_join(Loo, Lmax)
lengthratio[, "ratio"] <- lengthratio$Loo/lengthratio$Lmax
l <- length(lengthratio$ratio)
for (i in 1:l){
  if (!is.na(lengthratio$ratio[i])){
    if (lengthratio$ratio[i]>1.34 | lengthratio$ratio[i]<0.68){
      lengthratio$Loo[i] <- NA
    }
  }
}
newLooframe <- lengthratio[, c("Loo", "SpecCode")]
##################


#      MERGE DATA 
##############################################################################################
#merge the data
spemerged  <- full_join(framemat[, c(1,2)], framelength[, -which(colnames(framelength)=="Loo")], by="SpecCode")
dim(spemerged)
sum(duplicated(spemerged$SpecCode))
spemerged  <- full_join(spemerged, frameTL[, c("SpecCode", "TLDiet")], by="SpecCode")
dim(spemerged)
sum(duplicated(spemerged$SpecCode))
spemerged  <- full_join(spemerged, newLooframe, by="SpecCode")
dim(spemerged)
#sum(duplicated(spemerged$SpecCode))
#spemerged  <- full_join(spemerged, frameoff, by="SpecCode")   #not needed in the end
#dim(spemerged)
sum(duplicated(spemerged$SpecCode))
spemerged  <- full_join(spemerged, framefec[, c(1,2)], by="SpecCode")
dim(spemerged)
sum(duplicated(spemerged$SpecCode))
#spemerged  <- full_join(spemerged, swim[, c("SpecCode", "AdultMode")], by="SpecCode") #, "AspectRatio" #not needed
#dim(spemerged)
sum(duplicated(spemerged$SpecCode))
spemerged  <- full_join(spemerged, morpho_tot[, c("SpecCode", "Min_caudalpeduncle_depth", "Lower_jaw_length", "Max_body_depth", "Max_body_width")], by="SpecCode") 
dim(spemerged)
sum(duplicated(spemerged$SpecCode))
# spemerged  <é- full_join(spemerged, framear, by="SpecCode")
spemerged  <- full_join(spemerged, frameTemp, by="SpecCode")
dim(spemerged)
sum(duplicated(spemerged$SpecCode))
spemerged  <- left_join(spetot[,-which(colnames(spetot) == c("eps_m"))], spemerged, by="SpecCode")
dataset <- spemerged
dim(spemerged)
sum(duplicated(spemerged$SpecCode))
#View(dataset[which(duplicated(dataset$Species)),])
##############################################################################################




#      CLEAN DATA
##############################################################################################
#remove Loo and tmax measurement when Loo<Lm or tmax<tm
dataset$Loo[which(dataset[, c("Loo")]<dataset[, c("Lm")], arr.ind = T)]           <- NA
dataset$tmax[which(dataset[, c("tmax")]<dataset[, c("tm")],)] <- NA
#remove species that are not in the class Teleostean or Chondrychtian
dataset <- dataset[which(dataset$Class %in% c("Elasmobranchii", "Teleostei")),]
##############################################################################################



#       PREPARE DATA TO BE USED 
##############################################################################################
#### Transform the data set (dataandox) to fit the package Phylosem requirements : binomial
#######
#_____function_______
convertchardataintobinomial <- function(prefix, colname, datatot){
  data    <- datatot [, colname]
  nbc     <- length(table(data))
  nbl     <- length(data)
  dataNEW <- matrix(nrow=nbl, ncol=nbc)
  dataNEW <- lapply(1:nbc, function(i){
    name        <- names(table(data))[i]
    dataNEW[,i] <- data==name
  })
  names(dataNEW)<- paste0(prefix, names(table(data)))
  datatot <- cbind(dataNEW, datatot)
  datatot <- datatot[,-which(names(datatot) == colname)]
  return(datatot)
}
#Demerpela
dataset <- convertchardataintobinomial("habitat", "DemersPelag", dataset)
#######

#### Control dataset and remove data that is not realistic (<0 or =0 or NAN)
colnames(dataset)
boxplot(na.omit(dataset$Min_caudalpeduncle_depth))
boxplot(na.omit(dataset$Max_body_depth))
boxplot(na.omit(dataset$Max_body_width))
boxplot(na.omit(dataset$Lower_jaw_length))
dataset[which(dataset$Lower_jaw_length > 200),]
dataset[, c("Loo", "Lm", "K", "tmax", "Woo", "M", "TLDiet", "fecundity", "tm"  #"PeduncleDepth", "lengthOffspring", "AspectRatio"
            )] [dataset[, c("Loo", "Lm", "K", "tmax", "Woo", "M", "TLDiet", "fecundity", "tm" #"PeduncleDepth",  "lengthOffspring", "AspectRatio"
                                                          )]<=0] <- NA # constraint due to log transformation 
#### Last criterias to select data 
row.names(dataset)        <- chartr(" ", "_", dataset$Species)
dataset[sapply(dataset, is.nan)] <- NA
dataset$fecundity <- as.numeric(dataset$fecundity)
dataset$Brack     <- as.logical(dataset$Brack)
colnames(dataset)[which(colnames(dataset) == "Brack")] <- c("waterbrack")
dataset$Saltwater <- as.logical(dataset$Saltwater)
colnames(dataset)[which(colnames(dataset) == "Saltwater")] <- c("watersaltwater")

#######
#### Exploration : see correlation between traits to see which trait to log or not 
reg <- function(x, y, ...) {
  points(x,y, ...)
  abline(lm(y~x)) 
  abline(a=0,b=1,lty=1,col="red")
}# made to draw lowess line instead of regression line (lm)
panel.hist <- function(x, ...){
  #from help of pairs
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y)
}
datasetlog <- dataset
for (i in 1:dim(datasetlog)[2]){
  if (!(colnames(dataset)[i] %in% c("SpecCode", "Temperature", "TLDiet", "K", "M"))){
    if(is.numeric(datasetlog[,i])){
    datasetlog[,i] <- log(datasetlog[,i])
    }
  }
}

dir.create(paste0(pathoutput, "/plot_03"), showWarnings = F)
jpeg(file=paste0(pathoutput, "/plot_03/logornotlogplotHALREVISITED.jpeg"), width=1500, height=1000, res=100)
pairs(datasetlog[, c("Loo", "Lm", "K", "tmax", "Woo", "M", "TLDiet", "fecundity", "tm", "c_m",
                     "Min_caudalpeduncle_depth", "Max_body_width", "Max_body_depth", "Lower_jaw_length", "Temperature")], 
      col =  c(adjustcolor(4, .4), adjustcolor(3, .3))[as.factor(datasetlog$Class)], lower.panel = panel.smooth, upper.panel =reg, 
      diag.panel = panel.hist)
dev.off()
jpeg(file=paste0(pathoutput, "/plot_03/logornotlogplotnolog.jpeg"), width=1500, height=1000, res=100)
pairs(dataset[, c("Loo", "Lm", "K", "tmax", "Woo", "M", "TLDiet", "fecundity", "tm", "c_m",
                  "Min_caudalpeduncle_depth", "Max_body_width", "Max_body_depth", "Lower_jaw_length", "Temperature")], 
      col = c(adjustcolor(4, .4), adjustcolor(3, .3))[as.factor(dataset$Class)], lower.panel = panel.smooth, upper.panel =reg, 
      diag.panel = panel.hist)
dev.off()
datasetlog <- dataset
for (i in 1:dim(datasetlog)[2]){
    if(is.numeric(datasetlog[,i])){
      datasetlog[,i] <- log(datasetlog[,i])
    }
}
jpeg(file=paste0(pathoutput, "/plot_03/logornotlogplotTotlog.jpeg"), width=1500, height=1000, res=100)
pairs(datasetlog[, c("Loo", "Lm", "K", "tmax", "Woo", "M", "TLDiet", "fecundity", "tm", "c_m",
                     "Min_caudalpeduncle_depth", "Max_body_width", "Max_body_depth", "Lower_jaw_length", "Temperature")], 
      col = c(adjustcolor(4, .4), adjustcolor(3, .3))[as.factor(datasetlog$Class)], lower.panel = panel.smooth, upper.panel =reg, 
      diag.panel = panel.hist)
dev.off()
datalmfec <- na.omit(cbind(datasetlog$Lm, datasetlog$fecundity, datasetlog$Species))
plot(datalmfec, label=datalmfec[,3])
text(datalmfec,labels=datalmfec[,3])
#######

#*Transform the data set (dataandox) to fit the package requirements : log neparian numeric data
#########
for (i in 1:dim(dataset)[2]){
  if (!(colnames(dataset)[i] %in% c("SpecCode", "Temperature", "TLDiet"))){
    if(is.numeric(dataset[,i])){
      dataset[,i] <- log(dataset[,i])
    }
  }
    if(is.logical(dataset[,i])){
    dataset[,i] <- as.numeric(dataset[,i])
  }
}
summary(dataset)
############
##############################################################################################



#    2 DATASETS : DATASET TOT AND THE ONES ONLY WITH OSMOSE DATA + DATA WITH C_M
##############################################################################################
#################
#________________
#   DATASET TOT
dataset_TOT <- dataset
#create dataset_trait dataset with only the data needed for phylosem, without taxonomy
dataset_traits_TOT    <- dataset_TOT[, sapply(dataset_TOT, is.numeric)]
dataset_traits_TOT    <- dataset_traits_TOT[,-which(names(dataset_traits_TOT)==c("SpecCode"))]

#################
#________________
#  DATASET OSMOSE
#include only species that have data for c_m and eps_m or species from osmose
dataset_osm      <- dataset
rownames(dataset_osm) <- dataset$Species
dataset_ID       <- which(rownames(dataset_osm) %in% osmosespnames)
dataset_ID_gen   <- which(dataset_osm$Genus %in% dataset_osm$Genus[dataset_ID])

dataosmose       <- dataset_osm[dataset_ID_gen, ]
datatoadd        <- dataset_osm[-dataset_ID_gen,]
datatoadd        <- datatoadd[which(!is.na(datatoadd$c_m)), ]
dataset_GE_ge    <- rbind(datatoadd, dataosmose)
#create dataset_trait dataset with only the data needed for phylosem, without taxonomy
#osmose data : all genus matching with osmose data
dataset_traits_GE_ge <- dataset_GE_ge[, sapply(dataset_GE_ge, is.numeric)]
dataset_traits_GE_ge <- dataset_traits_GE_ge[,-which(names(dataset_traits_GE_ge)==c("SpecCode"))]

dataosmose       <- dataset_osm[dataset_ID, ]
datatoadd        <- dataset_osm[-dataset_ID,]
datatoadd        <- datatoadd[which(!is.na(datatoadd$c_m)), ]
dataset_GE       <- rbind(datatoadd, dataosmose)
#create dataset_trait dataset with only the data needed for phylosem, without taxonomy
#osmose data : only species from osmose
dataset_traits_GE <- dataset_GE[, sapply(dataset_GE, is.numeric)]
dataset_traits_GE <- dataset_traits_GE[,-which(names(dataset_traits_GE)==c("SpecCode"))]
##############################################################################################



#    VISUALIZE DATA
##############################################################################################
# plot the number of NAs
namena2 <- c()
na2 <- c()
j2=0
for (i in 1:dim(dataset_traits_TOT)[2]){
  j2=j2+1
  namena2[j2] <- colnames(dataset_traits_TOT[i]) 
  na2[j2] <- sum(is.na(dataset_traits_TOT[i]))
}
dataframe_na <- data.frame(namena2, na2)
ordna <- order(dataframe_na$na2)
a<-ggplot(dataframe_na, aes(x = namena2, y = na2)) +
  geom_col() +
  geom_text(aes(label = na2), vjust = -0.5, size = 4)+
  ylab("Number of NAs in dataset complete")

namena2 <- c()
na2 <- c()
j2=0
for (i in 1:dim(dataset_traits_GE)[2]){
  j2=j2+1
  namena2[j2] <- colnames(dataset_traits_GE[i]) 
  na2[j2] <- sum(is.na(dataset_traits_GE[i]))
}
dataframe_na <- data.frame(namena2, na2)
ordna <- order(dataframe_na$na2)
b<-ggplot(dataframe_na, aes(x = namena2, y = na2)) +
  geom_col() +
  geom_text(aes(label = na2), vjust = -0.5, size = 4)+
  ylab("Number of NAs in dataset for OSMOSE")

d<-ggarrange(a, b, ncol = 2, nrow = 1, labels = "AUTO",  common.legend = TRUE, legend="right")
ggsave(d, filename=paste0(pathoutput, "/plot_03/", "boxplot_NAsLOG.png"), width=c(30), height =c(15))


# count NAs per family
namena <- c()
na <- c()
na_cm <- c()
nb_fam <- c()
j=0
for (i in seq_along(table(dataset_GE$Family))){
  j=j+1
  namena[j] <- names(table(dataset_GE$Family))[i]
  nb_fam <- length(which(dataset_GE$Family == namena[j]))
  na[j] <- sum(is.na(dataset_GE[which(dataset_GE$Family == namena[j]), ]))/(nb_fam*dim(dataset_GE)[2])
  na_cm[j] <- sum(is.na(dataset_GE$c_m[which(dataset_GE$Family == namena[j])]))/nb_fam
}
dataframe_na <- data.frame(namena, na, na_cm)
ordna <- order(dataframe_na$na)
b<-ggplot(dataframe_na, aes(x = namena, y = na)) +
  geom_col() +
  geom_text(aes(label = na), vjust = -0.5, size = 4)+
  ylab("% of NAs per family")
c<-ggplot(dataframe_na, aes(x = namena, y = na_cm)) +
  geom_col() +
  geom_text(aes(label = na_cm), vjust = -0.5, size = 4)+
  ylab("% of NAs in c_m per family")
ggsave(b, filename=paste0(pathoutput, "/plot_03/", "boxplot_NAs_perfamily_GELOG.png"), width=c(30), height =c(15))
ggsave(c, filename=paste0(pathoutput, "/plot_03/", "boxplot_c_m_NAs_perfamily_GELOG.png"), width=c(30), height =c(15))


# count NAs per family
namena <- c()
na <- c()
na_cm <- c()
nb_fam <- c()
j=0
for (i in seq_along(table(dataset_TOT$Family))){
  j=j+1
  namena[j] <- names(table(dataset_TOT$Family))[i]
  nb_fam <- length(which(dataset_TOT$Family == namena[j]))
  na[j] <- sum(is.na(dataset_TOT[which(dataset_TOT$Family == namena[j]), ]))/(nb_fam*dim(dataset_TOT)[2])
  na_cm[j] <- sum(is.na(dataset_TOT$c_m[which(dataset_TOT$Family == namena[j])]))/nb_fam
}
dataframe_na <- data.frame(namena, na, na_cm)
ordna <- order(dataframe_na$na)
b<-ggplot(dataframe_na, aes(x = namena, y = na)) +
  geom_col() +
  geom_text(aes(label = na), vjust = -0.5, size = 4)+
  ylab("% of NAs per family")
c<-ggplot(dataframe_na, aes(x = namena, y = na_cm)) +
  geom_col() +
  geom_text(aes(label = na_cm), vjust = -0.5, size = 4)+
  ylab("% of NAs in c_m per family")
ggsave(b, filename=paste0(pathoutput, "/plot_03/", "boxplot_NAs_perfamily_TOTLOG.png"), width=c(30), height =c(15))
ggsave(c, filename=paste0(pathoutput, "/plot_03/", "boxplot_c_m_NAs_perfamily_TOTLOG.png"), width=c(30), height =c(15))


#osmose species witha value for c_m and eps_m
dataset_GE$Species[which(dataset_GE$Species[!is.na(dataset_GE$c_m)] %in% osmosespnames)]

sum(is.na(dataset_GE$c_m))
sum(is.na(dataset_GE$eps_m))
##############################################################################################




#    OUTPUT
##############################################################################################
##############
# write a csv file with this output
##############
options("scipen")

dir.create(paste0(pathoutput, "/dataset_for_phylosem_NOUNITCV/output_genus_stdmorpho"), showWarnings = F, recursive = T)
dir.create(paste0(pathoutput, "/dataset_for_phylosem_NOUNITCV/output_tot_stdmorpho"), showWarnings = F, recursive = T)

write.csv(dataset_GE, file = paste0(pathoutput, "/dataset_for_phylosem_NOUNITCV/output_genus_stdmorpho/dataset_phylosem.csv"))
write.csv(dataset_traits_GE, file = paste0(pathoutput, "/dataset_for_phylosem_NOUNITCV/output_genus_stdmorpho/dataset_traits_phylosem.csv"))

write.csv(dataset_TOT, file = paste0(pathoutput, "/dataset_for_phylosem_NOUNITCV/output_tot_stdmorpho/dataset_phylosem.csv"))
write.csv(dataset_traits_TOT, file = paste0(pathoutput, "/dataset_for_phylosem_NOUNITCV/output_tot_stdmorpho/dataset_traits_phylosem.csv"))

save.image( paste0(pathoutput, "/dataset_for_phylosem_NOUNITCV/IMAGELOG.RData"))
dev.off()
##############################################################################################
