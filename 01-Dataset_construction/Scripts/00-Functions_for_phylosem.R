
# packages
library(phylosem)
library(dplyr)
library(ape)
library(ggplot2)
library(ggtree)
library(phylopath)
library(stringr)
library(scales)
library(ggrepel)
require(grid)
library(ggpubr)
library(scales)
library(plotROC)
library(pROC)
library(fishtree)
library(rotl)
library(rsample)




#############################
############################
# Cross-validation ########

psem <- function(semID, trait, nameCVi){
  
  if (nameCVi == "all"){
    data_for_CV <- dataset_traits
    k_CV <- rsample::vfold_cv(as.data.frame(dataset), v=nbCV, repeats = 1)
    trait_name <- ""
  }
  if (nameCVi %in% c("c_m", "spe")){
    trait_name <- trait
    data_only_trait <- dataset_traits[-which(is.na(dataset_traits[,trait])), ]
    data_for_CV     <- dataset_traits[-which(is.na(dataset_traits[,trait])), ]
    totvar          <- var(data_only_trait[,trait])
    if (nameCVi  == "c_m"){
      data_Genus       <- dataset[-which(is.na(dataset[,trait])), ]
      genus_counts     <- table(data_Genus$Genus)
    }
    repeat {
      k_CV  <- rsample::vfold_cv(as.data.frame(data_only_trait), v=nbCV, repeats = 1)
      if (nameCVi  == "c_m"){k_CV  <- rsample::vfold_cv(as.data.frame(genus_counts), v=nbCV, repeats = 1)}
      CVvar  <- c()
      CVmed  <- c()
      for (i in 1:nbCV){
        if(nameCVi  == "c_m"){
          genus    <- rsample::analysis(k_CV$splits[[i]])$Var1
          sampleID <- which(data_Genus$Genus %in% genus)
        } else {
          sampleID <- k_CV$splits[[i]]$in_id 
        }
        CVvar[i] <- var(data_for_CV[-sampleID, trait])
        CVmed[i] <- median(data_for_CV[-sampleID, trait])
      }
      if (abs(100*max(abs(CVvar))/totvar - 100)<10 ){ break } #& abs(100*max(abs(CVmed))/totmed - 100)<10
    }
  }

  for (semIDi in semID){
    for (i in 1:nbCV){
      if(nameCVi  == "c_m"){
        genus      <- rsample::analysis(k_CV$splits[[i]])$Var1
        whichnotNA <- stringr::str_replace(data_Genus$Species[which(data_Genus$Genus %in% genus)], " ", "_")
      } else {
          whichnotNA <- rownames(rsample::analysis(k_CV$splits[[i]]))
      }
      

      data_CV <- dataset_traits
      data_CV[-which(rownames(data_CV) %in% whichnotNA),trait] <- NA
      psem = phylosem(sem = as.character(modellist[[semIDi]]),
                      data = data_CV,
                      tree = P,
                      family = c(rep("binomial", 3),  rep("fixed", 15) ),
                      estimate_ou = FALSE,
                      estimate_lambda = FALSE,
                      estimate_kappa = FALSE,
                      covs = colnames(data_CV)

      )
      matrixCSV[semIDi, i] <-  psem$opt$convergence
      assign(paste0("psemFINAL", "_rep", i, "_sem", semIDi), psem) # "_rep", i,
      df.1 = as_phylo4d(psem)
      df.2 = as(df.1, "data.frame")
      df.3 = df.2[order(df.2$label),]
      df.4 = df.3[df.3$node.type == "tip", ]
      df = df.4[,!names(df.4) %in% c("node", "ancestor", "edge.length", "node.type")]
      
      if (any(file.exists(paste0(pathoutput, "/CrossValidation/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCVi, trait_name, ".csv")))){
        file.remove(paste0(pathoutput, "/CrossValidation/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCVi, trait_name, ".csv"))
      }
      write.csv(df, paste0(pathoutput, "/CrossValidation/", "output", i, "_", modelname[semIDi], "psemFINAL", nameCVi, trait_name, ".csv"))


    }
    
  }
  
  write.csv(df, paste0(pathoutput, "/CrossValidation/", "Convergence", i, "_", "CV", nameCVi, trait_name, ".csv"))
  
  return(list(sample=k_CV, dataset_for_CV=data_for_CV))
}








# function to create a dataset with the cross validation outputs : data expected/data infered by phylosem/exected-infered
checkphylosemdata <- function(fileID, semID, trait, name, sample, maxCV){
  
  if (name != "all"){
    traitname = trait
    traittotest = trait
  }
  if (name == "all") {    
    traittotest = trait
    traitname = c("TOT")	
  }
  
  whichNA   <- sample[floor((1+floor(maxCV)*(fileID-1))):(floor(maxCV)*fileID)]
  if (maxCV!=floor(maxCV) & fileID==nbCV){whichNA<-c(whichNA, sample[length(sample)])}
  
  rest      <- dataset_traits[-whichnotNA,]
  data      <- read.csv(paste0(pathoutput, "/CrossValidation/output", fileID, "_", modelname[semID], "psemFINAL", name, traitname, ".csv"))
  
  # correct 'rest' considering the infered data is a ratio, and not the observed data (Lm and tm traits)
  # rest$Lm <- rest$Lm+rest$Loo
  # rest$tm <- log(exp(rest$tm)/rest$M)
  if(traittotest == "tm"){
    rest$tm <- rest$tm - rest$M
    data$tm <- data$tm - data$M
  }
  if(traittotest == "Lm"){
    rest$Lm <- rest$Lm + rest$Loo
    data$Lm <- data$Lm + data$Loo
  }
  
  #identify data to compare
  IDnonarest <- which(!is.na(rest[, traittotest]))
  species    <- rownames(rest[IDnonarest,])
  IDphylosem <- which(data$label %in% species)
  
  #extract data
  dataphylosem <- data[IDphylosem, c("label", traittotest)]
  names(dataphylosem) <- c("label", paste0(traittotest, "phylosem"))
  datafishbase <- rest[IDnonarest, c(traittotest) ]
  datafishbase <- cbind(datafishbase, rownames(rest[IDnonarest,]))
  colnames(datafishbase) <- c(traittotest, "label")
  
  dataphylosemextract <- dataphylosem[which(dataphylosem$label %in% datafishbase[,2]),]
  
  datacompare  <- full_join(dataphylosemextract, as.data.frame(datafishbase), by="label")
  datacompare[,traittotest] <- as.numeric(datacompare[,traittotest])
  
  #create comparaison column
  datacompare <- cbind(datacompare, ((datacompare[, traittotest]) - datacompare[, paste0(traittotest , "phylosem")] )*100/(datacompare[, traittotest]))
  colnames(datacompare) <- c(colnames(datacompare)[1:3], paste0(traittotest, "errorpercent"))
  
  return(datacompare)
}







#############################
############################
# CHANGE C_M AND EPS_M UNIT ###########
#*
#******explanation : ****************
#* 
#* We have data in Watt
#* we want it converted in mg O2 . kg^-1 . h^-1
#* 
#*


convert_watt_to_respi <- function(watt){
  # from watt (J * s^-1) to J * h^-1
  joules = watt*3600
  
  # J * h^-1 -> mmol of O2 * h^-1  ==> / 434 
  #* (Clarke & Johnston 1999)
  mmol = joules/434
  
  #* mmol of O2 * h^-1 -> mg of O2 * h^-1 ==> * 32
  #* (molar mass) 
  mg = mmol*32
  
  return(mg)
}



# minor species names corrections 
# clarke and Johnston 2025
ox99_correction <- function(oxygen99){
  oxygen99$Species[which(oxygen99$Species == "Anguilla anguillai")] <- c("Anguilla anguilla") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Barbus aeneus")] <- c("Labeobarbus aeneus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Crenimugil labrosus")] <- c("Chelon labrosus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Cyprinis carpio")] <- c("Cyprinus carpio") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Gymnelis viridis")] <- c("Gymnelus viridis") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Gymnoscopelus oplsthopterus")] <- c("Gymnoscopelus opisthopterus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Ictalurus nebucosus")] <- c("Ameiurus nebulosus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Liparis keofoedi")] <- c("Liparis fabricii") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Liparius atlanticus")] <- c("Liparis atlanticus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Liza dumerili")] <- c("Chelon dumerili") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Liza richardsonii")] <- c("Chelon richardsonii") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Mugil cephalas")] <- c("Mugil cephalus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Mystus cavensius")] <- c("Mystus cavasius") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Myxocephalus scorpius")] <- c("Myoxocephalus scorpius") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Notothenia angustifrons")] <- c("Gobionotothen angustifrons") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Notothenia nudifrons")] <- c("Nototheniops nudifrons") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Oreochromis alcalicus grahami")] <- c("Oreochromis grahami") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Oreochromis niltoicus")] <- c("Oreochromis niloticus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Oronchynchus nerka")] <- c("Oncorhynchus nerka") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Pseudocrenilabrius  multicolor")] <- c("Pseudocrenilabrus multicolor") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Salmo gairdneri")] <- c("Oncorhynchus mykiss") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Saratherodon mossambicus")] <- c("Oreochromis mossambicus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Tilapia mossambica")] <- c("Oreochromis mossambicus") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Tilapia nilotica")] <- c("Oreochromis niloticus") # correction wrong sp names
  # some species have no specified species name (just genus)
  # we are giving them false ones for the method to work. 
  # it will not impact the results, except for the position in the water column. We are hoping to have a match. 
  oxygen99$Species[which(oxygen99$Species == "Brachirus sp")] <- c("Brachirus orientalis") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Cynoglossus sp")] <- c("Cynoglossus gilchristi") # correction wrong sp names
  oxygen99$Species[which(oxygen99$Species == "Synaptura sp")] <- c("Synaptura commersonnii") # correction wrong sp names
  return(oxygen99)
}


ox16_correction <- function(oxygen16){
  # all the species names are right. 
  # some species have no specified species name (just genus)
  # we are giving them false ones for the method to work. 
  # it will not impact the results, except for the position in the water column. We are hoping to have a match. 
  oxygen16$Species[which(oxygen16$Species == "Ambassis sp.")] <- c("Ambassis gymnocephalus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Amblyeleotris sp.")] <- c("Amblyeleotris japonica") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Amblygobius sp.")] <- c("Amblygobius phalaena") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Ammodytes sp.")] <- c("Ammodytes marinus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Apogon sp.")] <- c("Apogon imberbis") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Callionymus sp.")] <- c("Callionymus lyra") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Caracanthus sp.")] <- c("Caracanthus typicus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Caulophrynidae sp.")] <- c("Caulophryne pelagica") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Cypselurus sp.")] <- c("Cypselurus naresii") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Galeoides sp.")] <- c("Galeoides decadactylus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Gerres sp.")] <- c("Gerres oyena") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Herklotsichthys sp.")] <- c("Herklotsichthys quadrimaculatus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Hypoatherina sp.")] <- c("Hypoatherina golanii") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Leiognathus sp.")] <- c("Leiognathus equula") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Lethrinus sp.")] <- c("Lethrinus atlanticus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Monacanthidae sp.")] <- c("Monacanthus chinensis") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Mullidae sp.")] <- c("Mulloidichthys ayliffe") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Neomyxus sp.")] <- c("Neomyxus leuciscus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Omobranchus sp.")] <- c("Omobranchus zebra") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Oneirodes sp.")] <- c("Oneirodes haplonema") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Pomacentrus sp.")] <- c("Pomacentrus caeruleus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Scomberesocidae sp.")] <- c("Cololabis adoceta") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Siganus sp.")] <- c("Siganus corallinus") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Sphyraena sp.")] <- c("Sphyraena japonica") # correction wrong sp names
  oxygen16$Species[which(oxygen16$Species == "Terapon sp.")] <- c("Terapon theraps") # correction wrong sp names
  return(oxygen16)
}



ox16b_correction <- function(oxygen16b){
  oxygen16b$Species[which(oxygen16b$Species == "Aphanius dispar")] <- c("Aphaniops dispar") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Aristichthys nobilis")] <- c("Hypophthalmichthys nobilis") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Crenimugil labrosus")] <- c("Chelon labrosus") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Epalzeorhynchos frenatus")] <- c("Epalzeorhynchos frenatum") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Gadus ogac")] <- c("Gadus macrocephalus") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Ictalurus nebulosus")] <- c("Ameiurus nebulosus") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Leiostomus xanthrus")] <- c("Leiostomus xanthurus") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Chelon aurata")] <- c("Liza aurata") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Macrozoraces americanus")] <- c("Zoarces americanus") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Notothenia nudifrons")] <- c("Nototheniops nudifrons") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Paracirrhites acratus")] <- c("Paracirrhites arcatus") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Stizostedion vitreum")] <- c("Sander vitreus") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Trematomus centronotus")] <- c("Trematomus pennellii") # correction wrong sp names
  # all the species names are right. 
  # some species have no specified species name (just genus)
  # we are giving them false ones for the method to work. 
  # it will not impact the results, except for the position in the water column. We are hoping to have a match. 
  oxygen16b$Species[which(oxygen16b$Species == "Hypophthalmichthys sp")] <- c("Hypophthalmichthys molitrix") # correction wrong sp names
  oxygen16b$Species[which(oxygen16b$Species == "Onychostoma sp")] <- c("Onychostoma elongatum") # correction wrong sp names
  
  return(oxygen16b)
}


ox21_correction <- function(oxygen21){
  oxygen21$ScientificName[which(oxygen21$ScientificName == "Gadus ogac")] <- c("Gadus macrocephalus") # correction wrong sp names
  oxygen21$ScientificName[which(oxygen21$ScientificName == "Leucoraja erinacea")] <- c("Leucoraja erinaceus") # correction wrong sp names
  return(oxygen21)
}


ox24_correction <- function(oxygen24){
  oxygen24$ScientificName[which(oxygen24$ScientificName == "Cephaloscyllium isabellum")] <- c("Cephaloscyllium isabella") # correction wrong sp names
  oxygen24$ScientificName[which(oxygen24$ScientificName == "Leucoraja erinacea")] <- c("Leucoraja erinaceus") # correction wrong sp names
  return(oxygen24)
}




##################
##################
#  ELSE    #######

osmosespnames <- c("Alosa alosa",                "Alosa fallax",                "Anguilla anguilla",           "Argyrosomus regius",          "Aristaeomorpha foliacea",
                   "Aristeus antennatus",         "Atherina boyeri",             "Auxis rochei",          "Belone belone",               "Boops boops",
                   "Caranx crysos",               "Chelidonichthys lucerna",     "Coris julis",                 "Coryphaena hippurus",         "Crangon crangon",
                   "Crystallogobius linearis",    "Dentex dentex",               "Dentex gibbosus",             "Dentex maroccanus",           "Dicentrarchus labrax",
                   "Diplodus annularis",          "Diplodus cervinus",           "Diplodus puntazzo",           "Diplodus sargus",       "Diplodus vulgaris",
                   "Eledone cirrhosa",            "Engraulis encrasicolus",      "Epinephelus aeneus",          "Epinephelus marginatus",      "Etrumeus sadina",
                   "Eutrigla gurnardus",          "Galeus melastomus",           "Gobius niger",                "Halobatrachus didactylus",    "Illex coindetii",
                   "Lepidorhombus whiffiagonis",  "Chelon auratus",                 "Chelon ramada",                 "Chelon saliens",                "Loligo vulgaris",
                   "Lophius budegassa",           "Lophius piscatorius",         "Merlangius merlangus",        "Merluccius merluccius",       "Micromesistius poutassou",
                   "Mugil cephalus",              "Mullus barbatus",             "Mullus surmuletus",           "Mustelus mustelus",           "Nephrops norvegicus",
                   "Octopus vulgaris",            "Pagellus acarne",             "Pagellus erythrinus",         "Pagrus pagrus",               "Palaemon serratus",
                   "Palinurus elephas",           "Parapenaeus longirostris",    "Penaeus kerathurus",          "Phycis phycis",               "Platichthys flesus",
                   "Pleuronectes platessa",       "Pomatomus saltatrix",         "Pomatoschistus marmoratus",   "Pomatoschistus minutus",      "Rhinobatos rhinobatos",
                   "Sarda sarda",                 "Sardina pilchardus",          "Sardinella aurita",           "Saurida undosquamis",         "Sciaena umbra",
                   "Scomber japonicus",          "Scomber scombrus",            "Scophthalmus maximus",        "Scorpaena notata",            "Scyliorhinus canicula",
                   "Sepia officinalis",           "Seriola dumerili",            "Serranus atricauda",          "Solea solea",                 "Sparus aurata",
                   "Sphyraena sphyraena",         "Sphyraena viridensis",        "Spicara maena",               "Spicara smaris",              "Spondyliosoma cantharus",
                   "Sprattus sprattus",           "Squilla mantis",              "Stephanolepis diaspros",      "Thunnus alalunga",            "Thunnus thynnus",
                   "Trachurus mediterraneus",     "Trachurus picturatus",        "Trachurus trachurus",         "Trachyrincus scabrus",        "Trigla lyra",
                   "Trisopterus luscus",          "Trisopterus minutus",         "Upeneus moluccensis",         "Xiphias gladius",             "Gobius ophiocephalus",
                   "Ost euphausiids")
