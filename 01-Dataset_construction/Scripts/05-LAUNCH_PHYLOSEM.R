###################### use dataset with phylosem : with all data from fishbase #######################
# needed  : output of the create_dataset_forphylosem_genus.R file, called : output_tot/dataset_traits_phylosem.csv
################################################################################################

#Load packages and files
path <- c("X:/PHYLOSEM")
path <- c("/home1/datahome/mbeneat/PHYLOSEM")
pathoutput <-  c("Z:/PHYLOSEM")
pathoutput <- c("/home1/scratch/mbeneat/PHYLOSEM/TLstdmeca")
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
#* # stdmeca
SEM = "
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


