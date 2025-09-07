
###################### use dataset with phylosem : with all data from fishbase #######################
# needed  : output of the create_dataset_forphylosem_genus.R file, called : output_tot/dataset_traits_phylosem.csv
################################################################################################

#Load packages and files
path <- c("X:/PHYLOSEM")
path <- c("/home1/datahome/mbeneat/PHYLOSEM")
pathoutput <-  c("Z:/PHYLOSEM")
pathoutput <- c("/home1/scratch/mbeneat/PHYLOSEM/TLstdLoo")
source(paste0(path, "/functionsphylosem.R"))

# prepare output folders 
dir.create(pathoutput, recursive = TRUE)
dir.create(paste0(pathoutput, "/else"))
dir.create(paste0(pathoutput, "/CrossValidation"))
dir.create(paste0(pathoutput, "/phylosem"))
dir.create(paste0(pathoutput, "/Phylogeny"))
dir.create(paste0(pathoutput, "/CrossValidation/plot"))
dir.create(paste0(pathoutput, "/CrossValidation/boxplot"))
source(paste0(path, "/functionsphylosem.R"))


# read dataset
dataset                  <- read.csv(paste0(path, "/output_tot/dataset_phylosem.csv"), row.names=c(1))
dataset_traits           <- read.csv(paste0(path, "/output_tot/dataset_traits_phylosem.csv"), row.names=c(1))

# clean dataset from unused traits
dataset_traits <- dataset_traits[,-which(colnames(dataset_traits) %in% 
                                           c("eps_m", "AspectRatio", "waterbrack", "watersaltwater", "lengthOffspring"))]

# WITH THE RATIO : convert Lm and tm into the corresponding and less variable 
# Beverton-Holt invariant trait 
boxplot(exp(dataset_traits$Lm)/exp(dataset_traits$Loo))
dataset_traits$Lm <- dataset_traits$Lm-dataset_traits$Loo
boxplot(exp(dataset_traits$tm)*(dataset_traits$M))
dataset_traits$tm <- dataset_traits$tm*dataset_traits$M

#*Define the SEM model
#*
#* # TLstdLoo
SEM =  "
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
####




# load tree
######
P = read.tree(paste0(path, "/Taxo_17003.tree"))


#*
#* RUN MODEL
#*

###########*************PARAMETERS
###########*
#trait = c("c_m") # choose the trait you want to cross validate ****
nbCV <- 10
rep  <- 1
semID=c(1)

modelname    <- c("SEM")
modellist    <- list(SEM)
matrixCSV    <- matrix("", 2,10)
matrixCSVtot <- matrix("", 2,1)
nameCV       <- c("c_m", "all", "spe")

###########*
###########*

#*
#* RUN PHYLOSEM 
#*

err <- matrix()
sem <- matrix()

for (semID in semID){
  data_CV <- dataext_traits
  psem = phylosem(sem = as.character(modellist[[semID]]),
                  data = data_CV,
                  tree = P,
                  # family = c(rep("fixed", 18)),
                  family = c(rep("binomial", 3),  rep("fixed", 15) ), #rep("binomial", 11),
                  estimate_ou = FALSE,
                  estimate_lambda = FALSE,
                  estimate_kappa = FALSE,
                  covs = colnames(data_CV)
                  # method="BFGS",newtonsteps = 1
  )
  matrixCSVtot[semID] <-  psem$opt$message
  assign(paste0("psemFINALtot", "_sem", semID), psem)
  
  # output 1 : dataset completed
  df.1 = as_phylo4d(psem)
  df.2 = as(df.1, "data.frame")
  df.3 = df.2[order(df.2$label),]
  df.4 = df.3[df.3$node.type == "tip", ]
  df = df.4[,!names(df.4) %in% c("node", "ancestor", "edge.length", "node.type")]
  labels  = stringr::str_replace_all(rownames(df), "_", " ")
  IDosmose= as.numeric(labels %in% osmosespnames)
  df$osmose <- IDosmose
  
  # output 2 : SEM and path coeffients
  semmodel = summary(psem)$coefficients
  #
  # output 3 : std and not std path coefficients
  coefnotstd <- coef(psem, standardized =FALSE)
  coefstd <- coef(psem, standardized =TRUE)
  coef = cbind(coefnotstd, coefstd[,3])
  colnames(coef) <- c("Path", "Parameter", "EstimateNotStd", "EstimateStd")
  
  
  if (file.exists(paste0(pathoutput, "/phylosem/", "output", "_", modelname[semID], "psemFINAL", ".csv"))){
    file.remove(paste0(pathoutput, "/phylosem/", "output", "_", modelname[semID], "psemFINAL", ".csv"))
  }
  if (file.exists(paste0(pathoutput, "/phylosem/", "c_mosmose", "_", modelname[semID], "psemFINAL", ".csv"))){
    file.remove(paste0(pathoutput, "/phylosem/", "c_mosmose", "_", modelname[semID], "psemFINAL", ".csv"))
  }
  if (file.exists(paste0(pathoutput, "/phylosem/", "semmodel", "_", ".csv"))){
    file.remove(paste0(pathoutput, "/phylosem/", "semmodel", "_", ".csv"))
  }
  if (file.exists(paste0(pathoutput, "/phylosem/", "MLE", "_", ".csv"))){
    file.remove(paste0(pathoutput, "/phylosem/", "MLE", "_", ".csv"))
  }
  
  write.csv(semmodel, paste0(pathoutput, "/phylosem/", "semmodel", "_", modelname[semID],".csv"))
  write.csv(coef, paste0(pathoutput, "/phylosem/", "coefnostd_std", "_", modelname[semID],".csv"))
  write.csv(df, paste0(pathoutput, "/phylosem/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
  # output4 : only OSMOSE values
  write.csv(df[IDosmose, c("c_m")], paste0(pathoutput, "/phylosem/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
  
  #* PLOT SEM WITH CORRELATION
  #*
  my_fitted_DAG = as_fitted_DAG(psem)
  coef_plot( my_fitted_DAG, error_bar = "se")
  pdf(file=paste0(pathoutput, "/phylosem/", modelname[semID], "SEM_DEG.pdf"), paper='a4r', width=200)
  plt <- plot(my_fitted_DAG, las = 2 )
  print(plt)
  dev.off()
}


CVconverged <- c(matrixCSV[1,], matrixCSV[2,], matrixCSVtot[1,], matrixCSVtot[2,])
CVconv <- data.frame(CVconverged)
write.csv(CVconv, paste0(pathoutput, "/phylosem/CSVconverged", ".csv"))

save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))



#*
#* CROSS-VALIDATION RUN
#*


##### RUN cross validation

# load(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))
# source(paste0(path, "/functionsphylosem"))

names_var_tot <- c("Benthopelagic", "Demersal", "Pelagic ", "MSM", "Maturation age", "Maturation size", 
                   "Max age", "K", "Infinite weight", "Mortality", "Trophic level", "Infinite length", "Fecundity",  "Peduncle depth", 
                   "Jaw length", "Body depth", "Body width", "Temperature")

# cross-validation on all traits at the same time
sampletot           <- psem(semID, trait=colnames(dataset_traits), nameCV[2]) #CV on all 
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on c_m only : genus by genus
samplec_mspe  <- psem(semID, trait="c_m", nameCV[1])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on Loo only
sampleLoospe <- psem(semID, trait="Loo", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on M only
sampleMspe <- psem(semID, trait="M", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatbenthopelagic only
samplehabitatbenthopelagicspe <- psem(semID, trait="habitatbenthopelagic", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatdemersal only
samplehabitatdemersalspe <- psem(semID, trait="habitatdemersal", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatpelagic only
samplehabitatpelagicspe <- psem(semID, trait="habitatpelagic", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on tm only
sampletmspe <- psem(semID, trait="tm", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatbenthopelagic only
sampleLmspe <- psem(semID, trait="Lm", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatdemersal only
sampleLmaxspe <- psem(semID, trait="Lmax", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatpelagic only
sampleKspe <- psem(semID, trait="K", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on M only
sampleWoospe <- psem(semID, trait="Woo", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatbenthopelagic only
sampleTLspe <- psem(semID, trait="TL", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatdemersal only
samplefecspe <- psem(semID, trait="fecundity", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatpelagic only
sampleMin_caudalpeduncle_depthspe <- psem(semID, trait="Min_caudalpeduncle_depth", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on M only
sampleLower_jaw_lengthspe <- psem(semID, trait="Lower_jaw_length", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatbenthopelagic only
sampleMax_body_widthspe <- psem(semID, trait="Max_body_width", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatdemersal only
sampleMax_body_depthspe <- psem(semID, trait="Max_body_depth", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))

# cross-validation on habitatdemersal only
sampleTemperaturespe <- psem(semID, trait="Temperature", nameCV[3])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))




# # #*
# # #* CROSS-VALIDATION CHECK & PLOTS
# # #*





# ncol=5
# nrow=4
# plot_checkphylosemdata(semID, trait=c("c_m", "K", "M", "Loo", "habitatpelagic", "habitatbenthopelagic", "habitatdemersal", "Woo"), name=nameCV[3],
#                        sample=list(samplec_mc_m$sample, sampleKspe$sample, sampleMspe$sample, sampleLoospe$sample,
#                                    samplehabitatpelagicspe$sample, samplehabitatbenthopelagicspe$sample, samplehabitatdemersalspe$sample,
#                                    sampleWoospe$sample),
#                        maxCV=list(samplec_mc_m$maxCV, sampleMspe$maxCV, sampleKspe$maxCV, sampleLoospe$maxCV,
#                                   samplehabitatpelagicspe$maxCV, samplehabitatbenthopelagicspe$maxCV, samplehabitatdemersalspe$maxCV,
#                                   sampleWoospe$maxCV), names_var)

# plot_checkphylosemdata(semID, trait=c("c_m", "K", "M", "Loo", "habitatpelagic", "habitatbenthopelagic", "habitatdemersal", "Woo"), name=nameCV[2], 
#                        sample=list(sampletot$sample), 
#                        maxCV=list(sampletot$maxCV), names_var)
# ncol = 4
# nrow = 5
# plot_checkphylosemdata(semID, trait=colnames(dataset_traits), name=nameCV[2], 
#                        sample=list(sampletot$sample), 
#                        maxCV=list(sampletot$maxCV), names_var)



# 
# # 
# 
# 
# #*
# #* RUN PHYLOSEM WITH THE BEST SEM: evol (semID = 2)
# #*
# 
# err <- matrix()
# sem <- matrix()
# 
# for (semID in semID){
#   data_CV <- dataset_traits
#   psem = phylosem(sem = as.character(modellist[[semID]]),
#                   data = data_CV,
#                   tree = P,
#                   family = c(rep("binomial", 3),  rep("fixed", 18) ), #rep("binomial", 11),
#                   estimate_ou = FALSE,
#                   estimate_lambda = FALSE,
#                   estimate_kappa = FALSE,
#                   covs = colnames(data_CV)
#                   # method="BFGS",newtonsteps = 1
#   )
#   matrixCSVtot[semID] <-  psem$opt$convergence
#   assign(paste0("psemFINALtot", "_sem", semID), psem)
#   
#   # output 1 : dataset completed
#   df.1 = as_phylo4d(psem)
#   df.2 = as(df.1, "data.frame")
#   df.3 = df.2[order(df.2$label),]
#   df.4 = df.3[df.3$node.type == "tip", ]
#   df = df.4[,!names(df.4) %in% c("node", "ancestor", "edge.length", "node.type")]
#   labels  = stringr::str_replace_all(rownames(df), "_", " ")
#   IDosmose= as.numeric(labels %in% osmosespnames)
#   df$osmose <- IDosmose
#   
#   # output 2 : SEM and path coeffients
#   semmodel = summary(psem)$coefficients
#   
#   # output 3 : std and not std path coefficients
#   coefnotstd <- coef(psem, standardized =FALSE)
#   coefstd <- coef(psem, standardized =TRUE)
#   coef = cbind(coefnotstd, coefstd[,3])
#   colnames(coef) <- c("Path", "Parameter", "EstimateNotStd", "EstimateStd")
#   
#   
#   if (file.exists(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
#   }
#   if (file.exists(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
#   }
#   if (file.exists(paste0(pathoutput, "/", "semmodel", "_", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "semmodel", "_", ".csv"))
#   }
#   if (file.exists(paste0(pathoutput, "/", "MLE", "_", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "MLE", "_", ".csv"))
#   }
#   
#   write.csv(semmodel, paste0(pathoutput, "/", "semmodel", "_", modelname[semID],".csv"))
#   write.csv(coef, paste0(pathoutput, "/", "coefnostd_std", "_", modelname[semID],".csv"))
#   write.csv(df, paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
#   # output4 : only OSMOSE values
#   write.csv(df[IDosmose, c("c_m")], paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
#   
#   #* PLOT SEM WITH CORRELATION
#   #*
#   my_fitted_DAG = as_fitted_DAG(psem)
#   coef_plot( my_fitted_DAG, error_bar = "se")
#   pdf(file=paste0(pathoutput, "/", modelname[semID], "SEM_DEG.pdf"), paper='a4r', width=200)
#   plt <- plot(my_fitted_DAG, las = 2 )
#   print(plt)
#   dev.off()
# }
# 
# 
# CVconverged <- c(matrixCSV[1,], matrixCSV[2,], matrixCSVtot[1,], matrixCSVtot[2,])
# CVconv <- data.frame(CVconverged)
# write.csv(CVconv, paste0(pathoutput, "/CSVconverged", ".csv"))
# 




#* 
#* COMPARE PHYLOGENETIC METHODS
#* 
# 

Grid = expand.grid( "OU" = c(FALSE,TRUE),
                    "lambda" = c(FALSE,TRUE),
                    "kappa" = c(FALSE,TRUE) )
psem_grid = NULL
for( i in 1:4){ #nrow(Grid)
  data_CV <- dataset_traits
  psem_grid[[i]] = phylosem(sem = as.character(modellist[[semID]]),
                            data = data_CV,
                            tree = P,
                            family = c(rep("binomial", 3),  rep("fixed", 18) ),
                            estimate_ou = Grid[i,'OU'],
                            estimate_lambda = Grid[i,'lambda'],
                            estimate_kappa = Grid[i,'kappa'],
                            quiet = TRUE,
                            covs = colnames(data_CV),
                            newtonsteps = 1)
}
# Extract AIC for each model and rank-order by parsimony
Grid$AIC = sapply( psem_grid, \(m) m$opt$AIC )
Grid = Grid[order(Grid$AIC,decreasing=FALSE),]

write.csv2(as.data.frame(Grid), paste0(pathoutput, "/Phylogeny/phylogeneticmethods.csv"))



# ###################### use dataext with phylosem : with all data from fishbase #######################
# # needed  : output of the create_dataext_forphylosem_genus.R file, called : output_tot/dataext_traits_phylosem.csv
# ################################################################################################
# 
# # #Load packages and files
# # path <- c("X:/phylosem_newtry_CHANGEDSEM")
# # path <- c("/home1/datahome/mbeneat/phylosem_newtry_CHANGEDSEM")
# # pathoutput <-  c("Z:/phylosem_newtry_CHANGEDSEM/tot_stdmorph_WooLoo_ratio")
# # pathoutput <- c("/home1/scratch/mbeneat/phylosem_newtry_CHANGEDSEM/tot_stdmorph_WooLoo_ratio")
# # 
# # 
# # load(file=paste0(pathoutput, "/", "imageworkspace.RData"))
# 
# #Load packages and files
# path <- c("X:/PHYLOSEM")
# path <- c("/home1/datahome/mbeneat/PHYLOSEM")
# pathoutput <-  c("Z:/PHYLOSEM")
# pathoutput <- c("/home1/scratch/mbeneat/PHYLOSEM")
# source(paste0(path, "/functionsphylosem.R"))
# 
# 
# osmosespnames <- c("Alosa alosa",                "Alosa fallax",                "Anguilla anguilla",           "Argyrosomus regius",          "Aristaeomorpha foliacea",
#                    "Aristeus antennatus",         "Atherina boyeri",             "Auxis rochei",          "Belone belone",               "Boops boops",
#                    "Caranx crysos",               "Chelidonichthys lucerna",     "Coris julis",                 "Coryphaena hippurus",         "Crangon crangon",
#                    "Crystallogobius linearis",    "Dentex dentex",               "Dentex gibbosus",             "Dentex maroccanus",           "Dicentrarchus labrax",
#                    "Diplodus annularis",          "Diplodus cervinus",           "Diplodus puntazzo",           "Diplodus sargus",       "Diplodus vulgaris",
#                    "Eledone cirrhosa",            "Engraulis encrasicolus",      "Epinephelus aeneus",          "Epinephelus marginatus",      "Etrumeus sadina",
#                    "Eutrigla gurnardus",          "Galeus melastomus",           "Gobius niger",                "Halobatrachus didactylus",    "Illex coindetii",
#                    "Lepidorhombus whiffiagonis",  "Chelon auratus",                 "Chelon ramada",                 "Chelon saliens",                "Loligo vulgaris",
#                    "Lophius budegassa",           "Lophius piscatorius",         "Merlangius merlangus",        "Merluccius merluccius",       "Micromesistius poutassou",
#                    "Mugil cephalus",              "Mullus barbatus",             "Mullus surmuletus",           "Mustelus mustelus",           "Nephrops norvegicus",
#                    "Octopus vulgaris",            "Pagellus acarne",             "Pagellus erythrinus",         "Pagrus pagrus",               "Palaemon serratus",
#                    "Palinurus elephas",           "Parapenaeus longirostris",    "Penaeus kerathurus",          "Phycis phycis",               "Platichthys flesus",
#                    "Pleuronectes platessa",       "Pomatomus saltatrix",         "Pomatoschistus marmoratus",   "Pomatoschistus minutus",      "Rhinobatos rhinobatos",
#                    "Sarda sarda",                 "Sardina pilchardus",          "Sardinella aurita",           "Saurida undosquamis",         "Sciaena umbra",
#                    "Scomber japonicus",          "Scomber scombrus",            "Scophthalmus maximus",        "Scorpaena notata",            "Scyliorhinus canicula",
#                    "Sepia officinalis",           "Seriola dumerili",            "Serranus atricauda",          "Solea solea",                 "Sparus aurata",
#                    "Sphyraena sphyraena",         "Sphyraena viridensis",        "Spicara maena",               "Spicara smaris",              "Spondyliosoma cantharus",
#                    "Sprattus sprattus",           "Squilla mantis",              "Stephanolepis diaspros",      "Thunnus alalunga",            "Thunnus thynnus",
#                    "Trachurus mediterraneus",     "Trachurus picturatus",        "Trachurus trachurus",         "Trachyrincus scabrus",        "Trigla lyra",
#                    "Trisopterus luscus",          "Trisopterus minutus",         "Upeneus moluccensis",         "Xiphias gladius",             "Gobius ophiocephalus",
#                    "Ost euphausiids")
# 
# 
# dataset                  <- read.csv(paste0(path, "/output_tot/dataset_phylosem.csv"))
# rownames(dataset)        <- dataset$X
# dataset_traits           <- read.csv(paste0(path, "/output_tot/dataset_traits_phylosem.csv"))
# rownames(dataset_traits) <- dataset_traits$X
# dataset                  <- dataset[, -which(grepl("X", colnames(dataset)))]
# dataset_traits           <- dataset_traits[, -which(grepl("X", colnames(dataset_traits)))]
# dataset_traits <- dataset_traits[,-which(colnames(dataset_traits) %in% 
#                                            c("eps_m", "AspectRatio", "waterbrack", "watersaltwater", "lengthOffspring"))]
# 
# 
# range(dataset$c_m, na.rm = T)
# range(dataset_traits$c_m, na.rm = T)
# colnames(dataset)
# colnames(dataset_traits)
# 
# 
# # cormat <- matrix(-99, dim(dataset_traits)[2],dim(dataset_traits)[2])
# # colnames(cormat) <- colnames(dataset_traits)
# # rownames(cormat) <- colnames(dataset_traits)
# # for (i in 1:dim(dataset_traits)[2]){
# #   a <- which(is.na(dataset_traits[,i]))
# #   for (j in 1:dim(dataset_traits)[2]){
# #     b <- which(is.na(dataset_traits[,j]))
# #     c <- as.numeric(names(table(c(a,b))))
# #     cormat[i,j]  <- cor(dataset_traits[-c,i], dataset_traits[-c, j])
# #   }
# # }
# # write.csv(cormat, paste0(pathoutput, "/cormatrix.csv"))
# 
# 
# # WITH THE RATIO
# boxplot(exp(dataset_traits$Lm)/exp(dataset_traits$Loo))
# dataset_traits$Lm <- dataset_traits$Lm-dataset_traits$Loo
# boxplot(exp(dataset_traits$tm)*(dataset_traits$M))
# dataset_traits$tm <- dataset_traits$tm*dataset_traits$M
# 
# 
# dataset -> dataext
# dataset_traits -> dataext_traits
# 
# 
# #*Define the SEM model
# #*
# #* # TLstdLoo
# SEM =  "
# habitatbenthopelagic -> c_m, p0
# habitatdemersal -> c_m, p1
# habitatpelagic -> c_m, p2
# Temperature-> K, p3
# Temperature-> M, p4
# Temperature-> Loo, p5
# Loo -> K, p6
# Loo -> M, p7
# Loo -> Max_body_width, p8
# Loo -> Max_body_depth, p9
# Loo -> Lower_jaw_length, p10
# Loo -> Min_caudalpeduncle_depth, p11
# Loo -> Woo, p12
# Loo -> c_m, p13
# Max_body_depth -> c_m, p14
# Max_body_width -> c_m,  p15
# Lower_jaw_length -> c_m, p16
# Min_caudalpeduncle_depth -> c_m, p17
# Max_body_depth -> TLDiet, p18
# Max_body_width -> TLDiet, p19
# Lower_jaw_length -> TLDiet, p20
# Min_caudalpeduncle_depth -> TLDiet, p21
# c_m -> K, p22
# c_m -> M, p23
# c_m -> Woo, p24
# K -> Lm, p25
# K -> tm, p26
# M -> tmax, p27
# M -> tm, p28
# M -> Lm, p29
# Woo -> TLDiet, p30
# Woo -> fecundity, p31
# "
# ####
# 
# range(dataext$c_m[!is.na(dataext$c_m)])
# range(dataext_traits$c_m[!is.na(dataext_traits$c_m)])
# 
# 
# 
# # phylogenetic tree
# ######
# nb         <- dim(dataext)[1]
# frm        = ~SuperClass/Class/Order/Family/Genus/Species
# phylo      <- c()
# phylo      = as.data.frame(dataext[which(names(dataext) %in% c("SuperClass", "Class", "Order", "Family", "Genus", "Species"))], stringsAsFactors = TRUE)
# phylo$Species        = stringr::str_replace(phylo$Species, " ", "_")
# phylo      <- mutate_if(phylo, is.character, as.factor)
# phylo_tree_med             = as.phylo(x = frm, data = phylo, collapse = FALSE, use.labels = TRUE)
# phylo_tree_med$edge.length = rep(1,length(phylo_tree_med$edge))
# 
# par(mfrow=c(1,1))
# write.tree(phy = phylo_tree_med, file = paste0(path, "/P.tree"), append = FALSE, digits = 10, tree.names = FALSE) # format understood by the *ape* package
# P = read.tree(paste0(path, "/P.tree"))
# if (any(is.na(P))) {
#   print("P has NA values - not allowed for")
# } else {
#   print("P looks OK") }
# 
# #*
# #* RUN MODEL
# #*
# 
# ###########*************PARAMETERS
# ###########*
# #trait = c("c_m") # choose the trait you want to cross validate ****
# nbCV <- 10
# semID=c(1)
# 
# listsem      <- list(TLstdmeca)
# max          <- length(rownames(dataext))
# IDNA_c_m     <- which(!is.na(dataext_traits$c_m))
# IDtot        <- length(dataext_traits$c_m)
# max          <- length(IDNA_c_m)/nbCV
# maxtot       <- IDtot/nbCV
# modelname    <- c("SEM")
# modellist    <-list(TLstdmeca)
# matrixCSV    <- matrix("", 2,10)
# matrixCSVtot <- matrix("", 2,1)
# nameCV       <- c("c_m", "all", "spe")
# 
# ###########*
# ###########*
# 
# 
# #*
# #* RUN PHYLOSEM 
# #*
# 
# err <- matrix()
# sem <- matrix()
# 
# for (semID in semID){
#   data_CV <- dataext_traits
#   psem = phylosem(sem = as.character(modellist[[semID]]),
#                   data = data_CV,
#                   tree = P,
#                   # family = c(rep("fixed", 18)),
#                   family = c(rep("binomial", 3),  rep("fixed", 18) ), #rep("binomial", 11),
#                   estimate_ou = FALSE,
#                   estimate_lambda = FALSE,
#                   estimate_kappa = FALSE,
#                   covs = colnames(data_CV)
#                   # method="BFGS",newtonsteps = 1
#   )
#   matrixCSVtot[semID] <-  psem$opt$message
#   assign(paste0("psemFINALtot", "_sem", semID), psem)
#   
#   # output 1 : dataset completed
#   df.1 = as_phylo4d(psem)
#   df.2 = as(df.1, "data.frame")
#   df.3 = df.2[order(df.2$label),]
#   df.4 = df.3[df.3$node.type == "tip", ]
#   df = df.4[,!names(df.4) %in% c("node", "ancestor", "edge.length", "node.type")]
#   labels  = stringr::str_replace_all(rownames(df), "_", " ")
#   IDosmose= as.numeric(labels %in% osmosespnames)
#   df$osmose <- IDosmose
#   
#   # output 2 : SEM and path coeffients
#   semmodel = summary(psem)$coefficients
#   #
#   # output 3 : std and not std path coefficients
#   coefnotstd <- coef(psem, standardized =FALSE)
#   coefstd <- coef(psem, standardized =TRUE)
#   coef = cbind(coefnotstd, coefstd[,3])
#   colnames(coef) <- c("Path", "Parameter", "EstimateNotStd", "EstimateStd")
#   
#   
#   if (file.exists(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
#   }
#   if (file.exists(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
#   }
#   if (file.exists(paste0(pathoutput, "/", "semmodel", "_", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "semmodel", "_", ".csv"))
#   }
#   if (file.exists(paste0(pathoutput, "/", "MLE", "_", ".csv"))){
#     file.remove(paste0(pathoutput, "/", "MLE", "_", ".csv"))
#   }
#   
#   write.csv(semmodel, paste0(pathoutput, "/", "semmodel", "_", modelname[semID],".csv"))
#   write.csv(coef, paste0(pathoutput, "/", "coefnostd_std", "_", modelname[semID],".csv"))
#   write.csv(df, paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
#   # output4 : only OSMOSE values
#   write.csv(df[IDosmose, c("c_m")], paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
#   
#   #* PLOT SEM WITH CORRELATION
#   #*
#   my_fitted_DAG = as_fitted_DAG(psem)
#   coef_plot( my_fitted_DAG, error_bar = "se")
#   pdf(file=paste0(pathoutput, "/", modelname[semID], "SEM_DEG.pdf"), paper='a4r', width=200)
#   plt <- plot(my_fitted_DAG, las = 2 )
#   print(plt)
#   dev.off()
# }
# 
# 
# CVconverged <- c(matrixCSV[1,], matrixCSV[2,], matrixCSVtot[1,], matrixCSVtot[2,])
# CVconv <- data.frame(CVconverged)
# write.csv(CVconv, paste0(pathoutput, "/CSVconverged", ".csv"))
# 
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# 
# 
# 
# #*
# #* CROSS-VALIDATION RUN
# #*
# 
# 
# #RUN cross validation
# 
# load(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# source(paste0(path, "/functionsphylosem.R"))
# 
# names_var_tot <- c("Benthopelagic", "Demersal", "Pelagic ", "MSM", "Maturation age", "Maturation size", 
#                    "Max age", "K", "Infinite weight", "Mortality", "Trophic level", "Infinite length", "Fecundity",  "Peduncle depth", 
#                    "Jaw length", "Body depth", "Body width", "Temperature")
# 
# sampletot <- psem(semID, trait=colnames(dataext_traits), nameCV[2])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# ncol=5
# nrow=4
# plot_checkphylosemdata(semID, trait=colnames(dataext_traits), name=nameCV[2], sample=list(sampletot$sample), maxCV=list(sampletot$maxCV), names_var=names_var_tot)
# 
# samplec_mspe  <- psem(semID, trait="c_m", nameCV[3])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("c_m"), name=nameCV[3], sample=list(samplec_mspe$sample), maxCV=list(samplec_mspe$maxCV), names_var=c("MSM"))
# 
# samplec_mc_m  <- psem(semID, trait="c_m", nameCV[1])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("c_m"), name=nameCV[1], sample=list(samplec_mc_m$sample), maxCV=list(samplec_mc_m$maxCV), names_var=c("MSM"))
# 
# sampleMspe <- psem(semID, trait="M", nameCV[3])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("M"), name=nameCV[3], sample=list(sampleMspe$sample), maxCV=list(sampleMspe$maxCV), names_var=c("Mortality"))
# 
# sampleLoospe <- psem(semID, trait="Loo", nameCV[3])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("Loo"), name=nameCV[3], sample=list(sampleLoospe$sample), maxCV=list(sampleLoospe$maxCV), names_var=c("Infinite length"))
# 
# sampleWoospe <- psem(semID, trait="Woo", nameCV[3])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("Woo"), name=nameCV[3], sample=list(sampleWoospe$sample), maxCV=list(sampleWoospe$maxCV), names_var=c("Infinite weight"))
# 
# samplehabitatbenthopelagicspe <- psem(semID, trait="habitatbenthopelagic", nameCV[3])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("habitatbenthopelagic"), name=nameCV[3], sample=list(samplehabitatbenthopelagicspe$sample), maxCV=list(samplehabitatbenthopelagicspe$maxCV), names_var=c("Benthopelagic"))
# 
# samplehabitatdemersalspe <- psem(semID, trait="habitatdemersal", nameCV[3])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("habitatdemersal"), name=nameCV[3], sample=list(samplehabitatdemersalspe$sample), maxCV=list(samplehabitatdemersalspe$maxCV), names_var=c("Demersal"))
# 
# samplehabitatpelagicspe <- psem(semID, trait="habitatpelagic", nameCV[3])
# save.image(file=paste0(pathoutput, "/", "imageworkspaceEND.RData"))
# plot_checkphylosemdata(semID, trait=c("habitatpelagic"), name=nameCV[3], sample=list(samplehabitatpelagicspe$sample), maxCV=list(samplehabitatpelagicspe$maxCV), names_var=c("Pelagic"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # # #*
# # # #* CROSS-VALIDATION CHECK & PLOTS
# # # #*
# 
# 
# 
# 
# 
# # ncol=5
# # nrow=4
# # plot_checkphylosemdata(semID, trait=c("c_m", "K", "M", "Loo", "habitatpelagic", "habitatbenthopelagic", "habitatdemersal", "Woo"), name=nameCV[3],
# #                        sample=list(samplec_mc_m$sample, sampleKspe$sample, sampleMspe$sample, sampleLoospe$sample,
# #                                    samplehabitatpelagicspe$sample, samplehabitatbenthopelagicspe$sample, samplehabitatdemersalspe$sample,
# #                                    sampleWoospe$sample),
# #                        maxCV=list(samplec_mc_m$maxCV, sampleMspe$maxCV, sampleKspe$maxCV, sampleLoospe$maxCV,
# #                                   samplehabitatpelagicspe$maxCV, samplehabitatbenthopelagicspe$maxCV, samplehabitatdemersalspe$maxCV,
# #                                   sampleWoospe$maxCV), names_var)
# 
# # plot_checkphylosemdata(semID, trait=c("c_m", "K", "M", "Loo", "habitatpelagic", "habitatbenthopelagic", "habitatdemersal", "Woo"), name=nameCV[2], 
# #                        sample=list(sampletot$sample), 
# #                        maxCV=list(sampletot$maxCV), names_var)
# # ncol = 4
# # nrow = 5
# # plot_checkphylosemdata(semID, trait=colnames(dataset_traits), name=nameCV[2], 
# #                        sample=list(sampletot$sample), 
# #                        maxCV=list(sampletot$maxCV), names_var)
# 
# 
# 
# # 
# # # 
# # 
# # 
# # #*
# # #* RUN PHYLOSEM WITH THE BEST SEM: evol (semID = 2)
# # #*
# # 
# # err <- matrix()
# # sem <- matrix()
# # 
# # for (semID in semID){
# #   data_CV <- dataext_traits
# #   psem = phylosem(sem = as.character(modellist[[semID]]),
# #                   data = data_CV,
# #                   tree = P,
# #                   family = c(rep("binomial", 3),  rep("fixed", 18) ), #rep("binomial", 11),
# #                   estimate_ou = FALSE,
# #                   estimate_lambda = FALSE,
# #                   estimate_kappa = FALSE,
# #                   covs = colnames(data_CV)
# #                   # method="BFGS",newtonsteps = 1
# #   )
# #   matrixCSVtot[semID] <-  psem$opt$convergence
# #   assign(paste0("psemFINALtot", "_sem", semID), psem)
# #   
# #   # output 1 : dataset completed
# #   df.1 = as_phylo4d(psem)
# #   df.2 = as(df.1, "data.frame")
# #   df.3 = df.2[order(df.2$label),]
# #   df.4 = df.3[df.3$node.type == "tip", ]
# #   df = df.4[,!names(df.4) %in% c("node", "ancestor", "edge.length", "node.type")]
# #   labels  = stringr::str_replace_all(rownames(df), "_", " ")
# #   IDosmose= as.numeric(labels %in% osmosespnames)
# #   df$osmose <- IDosmose
# #   
# #   # output 2 : SEM and path coeffients
# #   semmodel = summary(psem)$coefficients
# #   
# #   # output 3 : std and not std path coefficients
# #   coefnotstd <- coef(psem, standardized =FALSE)
# #   coefstd <- coef(psem, standardized =TRUE)
# #   coef = cbind(coefnotstd, coefstd[,3])
# #   colnames(coef) <- c("Path", "Parameter", "EstimateNotStd", "EstimateStd")
# #   
# #   
# #   if (file.exists(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))){
# #     file.remove(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
# #   }
# #   if (file.exists(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))){
# #     file.remove(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
# #   }
# #   if (file.exists(paste0(pathoutput, "/", "semmodel", "_", ".csv"))){
# #     file.remove(paste0(pathoutput, "/", "semmodel", "_", ".csv"))
# #   }
# #   if (file.exists(paste0(pathoutput, "/", "MLE", "_", ".csv"))){
# #     file.remove(paste0(pathoutput, "/", "MLE", "_", ".csv"))
# #   }
# #   
# #   write.csv(semmodel, paste0(pathoutput, "/", "semmodel", "_", modelname[semID],".csv"))
# #   write.csv(coef, paste0(pathoutput, "/", "coefnostd_std", "_", modelname[semID],".csv"))
# #   write.csv(df, paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
# #   # output4 : only OSMOSE values
# #   write.csv(df[IDosmose, c("c_m")], paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
# #   
# #   #* PLOT SEM WITH CORRELATION
# #   #*
# #   my_fitted_DAG = as_fitted_DAG(psem)
# #   coef_plot( my_fitted_DAG, error_bar = "se")
# #   pdf(file=paste0(pathoutput, "/", modelname[semID], "SEM_DEG.pdf"), paper='a4r', width=200)
# #   plt <- plot(my_fitted_DAG, las = 2 )
# #   print(plt)
# #   dev.off()
# # }
# # 
# # 
# # CVconverged <- c(matrixCSV[1,], matrixCSV[2,], matrixCSVtot[1,], matrixCSVtot[2,])
# # CVconv <- data.frame(CVconverged)
# # write.csv(CVconv, paste0(pathoutput, "/CSVconverged", ".csv"))
# # 
# 
# 
# 
# 
# #* 
# #* COMPARE PHYLOGENETIC METHODS
# #* 
# # 
# 
# Grid = expand.grid( "OU" = c(FALSE,TRUE),
#                     "lambda" = c(FALSE,TRUE),
#                     "kappa" = c(FALSE,TRUE) )
# psem_grid = NULL
# for( i in 1:4){ #nrow(Grid)
#   data_CV <- dataext_traits
#   psem_grid[[i]] = phylosem(sem = as.character(modellist[[semID]]),
#                             data = data_CV,
#                             tree = P,
#                             family = c(rep("binomial", 3),  rep("fixed", 18) ),
#                             estimate_ou = Grid[i,'OU'],
#                             estimate_lambda = Grid[i,'lambda'],
#                             estimate_kappa = Grid[i,'kappa'],
#                             quiet = TRUE,
#                             covs = colnames(data_CV),
#                             newtonsteps = 1)
# }
# # Extract AIC for each model and rank-order by parsimony
# Grid$AIC = sapply( psem_grid, \(m) m$opt$AIC )
# Grid = Grid[order(Grid$AIC,decreasing=FALSE),]
# 
# write.csv2(as.data.frame(Grid), paste0(pathoutput, "/phylogeneticmethods.csv"))
