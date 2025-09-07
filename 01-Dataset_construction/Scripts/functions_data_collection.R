#* CHANGE C_M AND EPS_M UNIT
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