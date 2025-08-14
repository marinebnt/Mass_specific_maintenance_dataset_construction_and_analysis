###################### use dataext with phylosem : with all data from fishbase #######################
# needed  : output of the create_dataext_forphylosem_genus.R file, called : output_tot/dataext_traits_phylosem.csv
################################################################################################

#COMPARE SEM with the AIC indicator


#Load packages and files
# path <- c("Z:/phylosem_newtry_CHANGEDSEM")
# path <- c("/home1/datahome/mbeneat/phylosem_newtry_CHANGEDSEM")
# pathoutput <-  c("Z:/phylosem_newtry_CHANGEDSEM")
# pathoutput <- c("/home1/scratch/mbeneat/phylosem_newtry_CHANGEDSEM")
pathoutput <- paste0("02-Analysis/Outputs/compareSEM")
dir.create(pathoutput)

source(paste0("01-Dataset_construction/Scripts/00-Functions_for_phylosem.R"))
# library(TMBhelper)
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


dataset        <- read.csv("01-Dataset_construction/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV/output_tot_stdmorpho/dataset_phylosem.csv")
rownames(dataset) <- dataset$X
dataset_traits <-  read.csv("01-Dataset_construction/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV/output_tot_stdmorpho/dataset_traits_phylosem.csv")
rownames(dataset_traits) <- dataset_traits$X


dataset        <- dataset[, -which(grepl("X", colnames(dataset)))]
dataset_traits <- dataset_traits[, -which(grepl("X", colnames(dataset_traits)))]

# cormat <- matrix(-99, dim(dataset_traits)[2],dim(dataset_traits)[2])
# colnames(cormat) <- colnames(dataset_traits)
# rownames(cormat) <- colnames(dataset_traits)
# for (i in 1:dim(dataset_traits)[2]){
#   a <- which(is.na(dataset_traits[,i]))
#   for (j in 1:dim(dataset_traits)[2]){
#     b <- which(is.na(dataset_traits[,j]))
#     c <- as.numeric(names(table(c(a,b))))
#     cormat[i,j]  <- cor(dataset_traits[-c,i], dataset_traits[-c, j])
#   }
# }
# write.csv(cormat, paste0(pathoutput, "/cormatrix.csv"))

dataset_traits <- dataset_traits[,-which(colnames(dataset_traits) %in% c("eps_m", "AspectRatio",
                                                                         "waterbrack", "watersaltwater", "lengthOffspring"))]
if(sum(grepl(pattern = "swim", colnames(dataset_traits))>0)){
  dataset_traits <- dataset_traits[,-grep(pattern = "swim", colnames(dataset_traits))]
}
# WITH THE RATIO
boxplot(exp(dataset_traits$Lm)/exp(dataset_traits$Loo))
dataset_traits$Lm <- dataset_traits$Lm-dataset_traits$Loo
boxplot(exp(dataset_traits$tm)*(dataset_traits$M))
dataset_traits$tm <- dataset_traits$tm+dataset_traits$M


dataset -> dataext
dataset_traits -> dataext_traits



#####
stdevol0 = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Min_caudalpeduncle_depth, p11
Min_caudalpeduncle_depth -> TLDiet, p12
Max_body_depth-> TLDiet, p13
Max_body_width-> TLDiet,  p14
Lower_jaw_length-> TLDiet, p15
Loo -> Woo, p16
K -> Lm, p17
K -> tm, p18
M -> Lm, p19
M -> tmax, p20
M -> tm, p22
M -> c_m, p23
Max_body_depth-> c_m, p24
Max_body_width-> c_m,  p25
Lower_jaw_length-> c_m, p26
Min_caudalpeduncle_depth  -> c_m, p27
Woo -> c_m, p28
Woo -> fecundity, p29
Woo -> TLDiet, p30
"
#####

#####
stdevol00 = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Min_caudalpeduncle_depth, p11
Min_caudalpeduncle_depth -> TLDiet, p12
Loo -> Woo, p13
K -> Lm, p14
K -> tm, p15
M -> Lm, p16
M -> tmax, p17
M -> tm, p18
M -> c_m, p19
Max_body_depth-> c_m, p20
Max_body_width-> c_m,  p21
Lower_jaw_length-> c_m, p22
Min_caudalpeduncle_depth  -> c_m, p23
Woo -> c_m, p24
Woo -> fecundity, p25
Woo -> TLDiet, p26
"
#####

#####
stdevol1 = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Woo, p11
K -> Lm, p12
K -> tm, p13
M -> Lm, p14
M -> tmax, p15
M -> tm, p16
M -> c_m, p17
Max_body_depth-> c_m, p18
Max_body_width-> c_m,  p19
Lower_jaw_length-> c_m, p20
Min_caudalpeduncle_depth  -> c_m, p21
Min_caudalpeduncle_depth  -> TLDiet, p22
Woo -> c_m, p23
Woo -> fecundity, p24
Woo -> TLDiet, p25
"
#####

#####
stdevol2 = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Woo, p11
K -> Lm, p12
K -> tm, p13
M -> Lm, p14
M -> tmax, p15
M -> tm, p16
M -> c_m, p17
Max_body_depth-> c_m, p18
Max_body_width-> c_m,  p19
Lower_jaw_length-> c_m, p20
Min_caudalpeduncle_depth  -> c_m, p21
Woo -> c_m, p23
Woo -> fecundity, p24
Woo -> TLDiet, p25
"
#####

stdmeca = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Min_caudalpeduncle_depth, p11
Loo -> Woo, p12
Max_body_depth-> c_m, p13
Max_body_width-> c_m,  p14
Lower_jaw_length-> c_m, p15
Min_caudalpeduncle_depth  -> c_m, p16
c_m -> K, p18
c_m -> M, p19
c_m -> Woo, p20
K -> Lm, p21
K -> tm, p22
M -> tmax, p23
M -> tm, p24
M -> Lm, p25
Woo -> TLDiet, p26
Woo -> fecundity, p27
"
####

TLstdevol = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Min_caudalpeduncle_depth, p11
Loo -> Woo, p12
K -> Lm, p13
K -> tm, p14
M -> Lm, p15
M -> tmax, p16
M -> tm, p17
M -> c_m, p18
Max_body_depth-> c_m, p19
Max_body_width-> c_m,  p20
Lower_jaw_length-> c_m, p21
Min_caudalpeduncle_depth  -> c_m, p22
Max_body_depth-> TLDiet, p23
Max_body_width-> TLDiet, p24
Lower_jaw_length-> TLDiet, p25
Min_caudalpeduncle_depth  -> TLDiet, p26
Woo -> c_m, p27
Woo -> fecundity, p28
Woo -> TLDiet, p29
"
#####

#MODEL 2 -> MECANISTIC
TLstdmeca = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Min_caudalpeduncle_depth, p11
Loo -> Woo, p12
Max_body_depth-> c_m, p13
Max_body_width-> c_m,  p14
Lower_jaw_length-> c_m, p15
Min_caudalpeduncle_depth  -> c_m, p16
Max_body_depth-> TLDiet, p17
Max_body_width-> TLDiet,  p18
Lower_jaw_length-> TLDiet, p19
Min_caudalpeduncle_depth  -> TLDiet, p20
c_m -> K, p21
c_m -> M, p22
c_m -> Woo, p23
K -> Lm, p24
K -> tm, p25
M -> tmax, p26
M -> tm, p27
M -> Lm, p28
Woo -> TLDiet, p29
Woo -> fecundity, p30
"
####

#####
TLWoo = "
habitatbenthopelagic -> c_m, p1
habitatdemersal -> c_m, p2
habitatpelagic -> c_m, p3
Temperature-> K, p4
Temperature-> M, p5
Max_body_depth-> TLDiet, p6
Max_body_width-> TLDiet, p7
Lower_jaw_length-> TLDiet, p8
Min_caudalpeduncle_depth  -> TLDiet, p9
Max_body_depth-> c_m, p10
Max_body_width-> c_m,  p11
Lower_jaw_length-> c_m, p12
Min_caudalpeduncle_depth -> c_m, p13
c_m -> K, p14
c_m -> M, p15
c_m -> Loo, p16
K -> Lm, p17
K -> tm, p18
M -> tmax, p19
M -> tm, p20
M -> Lm, p21
Loo -> TLDiet, p22
Loo -> fecundity, p23
"
####

#####
TLWooWoo = "
habitatbenthopelagic -> c_m, p1
habitatdemersal -> c_m, p2
habitatpelagic -> c_m, p3
Temperature-> K, p4
Temperature-> M, p5
Max_body_depth-> TLDiet, p6
Max_body_width-> TLDiet, p7
Lower_jaw_length-> TLDiet, p8
Min_caudalpeduncle_depth  -> TLDiet, p9
Max_body_depth-> c_m, p10
Max_body_width-> c_m,  p11
Lower_jaw_length-> c_m, p12
Min_caudalpeduncle_depth -> c_m, p13
c_m -> K, p14
c_m -> M, p15
c_m -> Loo, p16
K -> Lm, p17
K -> tm, p18
M -> tmax, p19
M -> tm, p20
M -> Lm, p21
Loo -> TLDiet, p22
Loo -> fecundity, p23
Loo -> Woo, p24
"
#####

Woo = "
habitatbenthopelagic -> c_m, p1
habitatdemersal -> c_m, p2
habitatpelagic -> c_m, p3
Temperature-> K, p4
Temperature-> M, p5
Max_body_depth-> c_m, p6
Max_body_width-> c_m,  p7
Lower_jaw_length-> c_m, p8
Min_caudalpeduncle_depth-> c_m, p9
c_m -> K, p10
c_m -> M, p11
c_m -> Loo, p12
K -> Lm, p13
K -> tm, p14
M -> tmax, p15
M -> tm, p16
M -> Lm, p17
Loo -> TLDiet, p18
Loo -> fecundity, p19
"
#####
WooWoo = "
habitatbenthopelagic -> c_m, p1
habitatdemersal -> c_m, p2
habitatpelagic -> c_m, p3
Temperature-> K, p4
Temperature-> M, p5
Max_body_depth-> c_m, p6
Max_body_width-> c_m,  p7
Lower_jaw_length-> c_m, p8
Min_caudalpeduncle_depth-> c_m, p9
c_m -> K, p10
c_m -> M, p11
c_m -> Loo, p12
K -> Lm, p13
K -> tm, p14
M -> tmax, p15
M -> tm, p16
M -> Lm, p17
Loo -> TLDiet, p18
Loo -> fecundity, p19
Loo -> Woo, p20
"
####
stdLoo = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Min_caudalpeduncle_depth, p11
Loo -> Woo, p12
Loo -> c_m, p13
Max_body_depth -> c_m, p14
Max_body_width -> c_m,  p15
Lower_jaw_length -> c_m, p16
Min_caudalpeduncle_depth -> c_m, p17
Min_caudalpeduncle_depth -> TLDiet, p18
c_m -> K, p19
c_m -> M, p20
c_m -> Woo, p21
K -> Lm, p22
K -> tm, p23
M -> tmax, p24
M -> tm, p25
M -> Lm, p26
Woo -> TLDiet, p27
Woo -> fecundity, p28
"
####


#####
TLstdLoo = "
habitatbenthopelagic -> c_m, p0
habitatdemersal -> c_m, p1
habitatpelagic -> c_m, p2
Temperature-> K, p3
Temperature-> M, p4
Temperature-> Loo, p5
Loo -> K, p6
Loo -> M, p7
Loo -> Max_body_width, p8
Loo -> Max_body_depth, p9
Loo -> Lower_jaw_length, p10
Loo -> Min_caudalpeduncle_depth, p11
Loo -> Woo, p12
Loo -> c_m, p13
Max_body_depth -> c_m, p14
Max_body_width -> c_m,  p15
Lower_jaw_length -> c_m, p16
Min_caudalpeduncle_depth -> c_m, p17
Max_body_depth -> TLDiet, p18
Max_body_width -> TLDiet, p19
Lower_jaw_length -> TLDiet, p20
Min_caudalpeduncle_depth -> TLDiet, p21
c_m -> K, p22
c_m -> M, p23
c_m -> Woo, p24
K -> Lm, p25
K -> tm, p26
M -> tmax, p27
M -> tm, p28
M -> Lm, p29
Woo -> TLDiet, p30
Woo -> fecundity, p31
"

#####

Woo2 = "
habitatbenthopelagic -> c_m, p1
habitatdemersal -> c_m, p2
habitatpelagic -> c_m, p3
Temperature-> K, p4
Temperature-> M, p5
Max_body_depth-> c_m, p6
Max_body_width-> c_m,  p7
Lower_jaw_length-> c_m, p8
Min_caudalpeduncle_depth-> c_m, p9
Min_caudalpeduncle_depth-> TLDiet, p10
c_m -> K, p11
c_m -> M, p12
c_m -> Loo, p13
K -> Lm, p14
K -> tm, p15
M -> tmax, p16
M -> tm, p17
M -> Lm, p18
Loo -> TLDiet, p19
Loo -> fecundity, p20
"
#####
WooWoo2 = "
habitatbenthopelagic -> c_m, p1
habitatdemersal -> c_m, p2
habitatpelagic -> c_m, p3
Temperature-> K, p4
Temperature-> M, p5
Max_body_depth-> c_m, p6
Max_body_width-> c_m,  p7
Lower_jaw_length-> c_m, p8
Min_caudalpeduncle_depth-> c_m, p9
Min_caudalpeduncle_depth-> TLDiet, p10
c_m -> K, p11
c_m -> M, p12
c_m -> Loo, p13
K -> Lm, p14
K -> tm, p15
M -> tmax, p16
M -> tm, p17
M -> Lm, p18
Loo -> TLDiet, p19
Loo -> fecundity, p20
Loo -> Woo, p21
"


####

# dataext <- dataext[1:200,]
# dataext_traits <- dataext_traits[1:200,]

# phylogenetic tree
######
nb         <- dim(dataext)[1]
frm        = ~SuperClass/Class/Order/Family/Genus/Species
phylo      <- c()
phylo      = as.data.frame(dataext[which(names(dataext) %in% c("SuperClass", "Class", "Order", "Family", "Genus", "Species"))], stringsAsFactors = TRUE)
phylo$Species        = stringr::str_replace(phylo$Species, " ", "_")
phylo      <- mutate_if(phylo, is.character, as.factor)
phylo_tree_med             = as.phylo(x = frm, data = phylo, collapse = FALSE, use.labels = TRUE)
phylo_tree_med$edge.length = rep(1,length(phylo_tree_med$edge))

par(mfrow=c(1,1))
write.tree(phy = phylo_tree_med, file = paste0(pathoutput, "/P.tree"), append = FALSE, digits = 10, tree.names = FALSE) # format understood by the *ape* package
P = read.tree(paste0(pathoutput, "/P.tree"))
if (any(is.na(P))) {
  print("P has NA values - not allowed for")
} else {
  print("P looks OK") }


#*
#* RUN MODEL
#*

###########*************PARAMETERS
###########*
#trait = c("c_m") # choose the trait you want to cross validate ****
nbCV <- 10
semID=c(2)

listsem <- list(stdevol1, stdevol2, stdmeca, TLstdevol, TLstdmeca, TLWoo, TLWooWoo, Woo, WooWoo, TLstdLoo) #stdLoo
listsem <- list(stdmeca, TLstdmeca, stdevol00, stdevol0, Woo2, TLWoo, WooWoo2, TLWooWoo, stdLoo, TLstdLoo)#, stdLoo, TLstdLoo)
# listsem <- list(stdevol0)
max     <- length(rownames(dataext))
IDNA_c_m<- which(!is.na(dataext_traits$c_m))
IDtot   <- length(dataext_traits$c_m)
max25   <- length(IDNA_c_m)/4
max33   <- length(IDNA_c_m)/nbCV
max25tot   <- IDtot/4
max33tot   <- IDtot/nbCV
modelname <- c("stdevol1", "stdevol2", "stdmeca", "TLstdevol", "TLstdmeca", "TLWoo", "TLWooWoo", "Woo", "WooWoo", "TLstdLoo")
modelname <- c("stdmeca", "TLstdmeca", "stdevol00", "stdevol0", "Woo2", "TLWoo", "WooWoo2", "TLWooWoo", "stdLoo", "TLstdLoo")#, "stdLoo", "TLstdLoo")
# modelname <- c("stdevol0")
#modellist <-list(modelmeca, modelmeca)
matrixCSV <- matrix("", 2,10)
matrixCSVtot <- matrix("", 2,1)
nameCV <- c("c_m", "all", "spe")
ncol = 3
nrow = 4 
###########*
###########*

comparison <- compare_phylosem(
  sem_set=listsem, #[5:length(listsem)]
  tree=P,
  data=dataext_traits,
  family =c(rep("binomial", 3),  rep("fixed", 18) ),
  covs=colnames(dataext_traits),
  estimate_ou = FALSE,
  estimate_lambda = FALSE,
  estimate_kappa = FALSE,
  control = phylosem_control()
)

# phylosem::phylosem(
#   sem=listsem[[1]], #[5:length(listsem)]
#   tree=P,
#   data=dataext_traits,
#   family =c(rep("binomial", 3),  rep("fixed", 18) ),
#   covs=colnames(dataext_traits),
#   estimate_ou = FALSE,
#   estimate_lambda = FALSE,
#   estimate_kappa = FALSE,
#   control = phylosem_control(),
# )

cat("\n---------->\n")
best.compare_phylosem <-
  function( x ) {
    AICs <- vapply(x, AIC, 1)
    x[[which.min(AICs)]]
    return(AICs)
  }
best.compare_phylosem(comparison)

save.image(paste0(pathoutput, "/comparesemLAST.RData"))





#06/11
"stdmeca",
114329.3 
"TLstdmeca", ****
  114276.0 
"stdevol00", 
114391.1 
"stdevol0", 
114369.0 
"Woo2", 
117379.3 
"TLWoo", 
117347.6 
"WooWoo2",
114574.5 
"TLWooWoo", 
114543.1 
"stdLoo", 
114316.3 
"TLstdLoo"
114294.7


# 29/04
   "stdmeca", "TLstdmeca", "stdevol00", "stdevol0", "Woo2", "TLWoo", "WooWoo2", "TLWooWoo", "stdLoo", "TLstdLoo"
    115913.6  115893.3     115924.2     115902.1    119500.0   119463.3   116774.2 116738.0 115917.3   115895.1
