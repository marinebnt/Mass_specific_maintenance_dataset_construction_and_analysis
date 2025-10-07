###################### use dataset with phylosem : with all data from fishbase #######################
# needed  : output of the create_dataset_forphylosem_genus.R file, called : output_tot/dataset_traits_phylosem.csv
################################################################################################

# Write paths and load dependencies

path <- paste0("01-Dataset_construction/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV")
pathoutput <- paste0( "01-Dataset_construction/Outputs/phylosem_alldatasets")
source(paste0("01-Dataset_construction/Scripts/00-Functions_for_phylosem.R"))
dir.create(pathoutput, recursive = T)

# load dataset
dataset        <- read.csv(paste0(path, "/output_tot_stdmorpho/dataset_phylosem.csv"), row.names = c(1))
dataset_traits <-  read.csv(paste0(path, "/output_tot_stdmorpho/dataset_traits.csv"), row.names = c(1))
rownames(dataset) <- stringr::str_replace(rownames(dataset), pattern = " ", replacement = "_")
rownames(dataset_traits) <- stringr::str_replace(rownames(dataset_traits), pattern = " ", replacement = "_")

if(length(which(dataset$Genus %in% c("Fundulus")))>0) {dataset_traits <- dataset_traits[-which(dataset$Genus %in% c("Fundulus")),]}
if(length(which(dataset$Genus %in% c("Fundulus")))>0) {dataset <- dataset[-which(dataset$Genus %in% c("Fundulus")),]}


# remove useless traits
dataset_traits <- dataset_traits[,-which(colnames(dataset_traits) %in% c("eps_m", "AspectRatio",
                                                                         "waterbrack", "watersaltwater", "lengthOffspring"))]

# replacing maturation traits with invariant traits
# WITH THE RATIO (INVARIANTS beverton & holt)
boxplot(exp(dataset_traits$Lm)/exp(dataset_traits$Loo))
dataset_traits$Lm <- dataset_traits$Lm-dataset_traits$Loo
boxplot(exp(dataset_traits$tm)*exp(dataset_traits$M))
dataset_traits$tm <- log(exp(dataset_traits$tm)*exp(dataset_traits$M))

if(any(is.nan(dataset_traits$tm))) {dataset_traits$tm[which(is.nan(dataset_traits$tm))] <- NA}
if(any(is.nan(dataset_traits$Lm))) {dataset_traits$Lm[which(is.nan(dataset_traits$Lm))] <- NA}


#*Define the SEM model
#*
#*
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
###





# phylogenetic tree
######
nb         <- dim(dataset)[1]
frm        = ~SuperClass/Class/Order/Family/Genus/Species
phylo      <- c()
phylo      = as.data.frame(dataset[which(names(dataset) %in% c("SuperClass", "Class", "Order", "Family", "Genus", "Species"))], stringsAsFactors = TRUE)
phylo$Species        = stringr::str_replace(phylo$Species, " ", "_")
phylo      <- mutate_if(phylo, is.character, as.factor)
phylo_tree_med             = as.phylo(x = frm, data = phylo, collapse = FALSE, use.labels = TRUE)
phylo_tree_med$edge.length = rep(1,length(phylo_tree_med$edge))

par(mfrow=c(1,1))
write.tree(phy = phylo_tree_med, file = paste0(path, "/P.tree"), append = FALSE, digits = 10, tree.names = FALSE) # format understood by the *ape* package
P = read.tree(paste0(path, "/P.tree"))
if (any(is.na(P))) {
  print("P has NA values - not allowed for")
} else {
  print("P looks OK") }



#*
#* RUN MODEL
#*

###########*************PARAMETERS
###########*

nbCV <- 10 # how many steps to cross-validation ?
semID=c(1) # if you have a list of SEM, specify on which one you want to run the analysis

listsem <- list(TLstdmeca) # list of SEM models
max     <- length(rownames(dataset)) 
IDNA_c_m<- which(!is.na(dataset_traits$c_m))
IDtot   <- length(dataset_traits$c_m)
max25   <- length(IDNA_c_m)/4
max   <- length(IDNA_c_m)/nbCV
max25tot   <- IDtot/4
maxtot   <- IDtot/nbCV
modelname <- c("TLstdmeca")
modellist <-list(TLstdmeca)
matrixCSV <- matrix("", 2,10)
matrixCSVtot <- matrix("", 2,1)
nameCV <- c("c_m", "all", "spe")

###########*
###########*


#*
#* RUN PHYLOSEM
#*

err <- matrix()
sem <- matrix()

for (semID in semID){
  data_CV <- dataset_traits
  psem = phylosem(sem = as.character(modellist[[semID]]),
                  data = data_CV,
                  tree = P,
                  family = c(rep("binomial", 3),  rep("fixed", 15) ),
                  estimate_ou = FALSE,
                  estimate_lambda = FALSE,
                  estimate_kappa = FALSE,
                  covs = colnames(data_CV)
  )
  matrixCSVtot[semID] <-  psem$opt$convergence
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


  if (file.exists(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))){
    file.remove(paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
  }
  if (file.exists(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))){
    file.remove(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
  }
  if (file.exists(paste0(pathoutput, "/", "semmodel", "_", ".csv"))){
    file.remove(paste0(pathoutput, "/", "semmodel", "_", ".csv"))
  }
  if (file.exists(paste0(pathoutput, "/", "MLE", "_", ".csv"))){
    file.remove(paste0(pathoutput, "/", "MLE", "_", ".csv"))
  }

  write.csv(semmodel, paste0(pathoutput, "/", "semmodel", "_", modelname[semID],".csv"))
  write.csv(coef, paste0(pathoutput, "/", "coefnostd_std", "_", modelname[semID],".csv"))
  write.csv(df, paste0(pathoutput, "/", "output", "_", modelname[semID], "psemFINALtot", ".csv"))
  # output4 : only OSMOSE values
  write.csv(df[IDosmose, c("c_m")], paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))

  #* PLOT SEM WITH CORRELATION
  #*
  my_fitted_DAG = as_fitted_DAG(psem)
  coef_plot( my_fitted_DAG, error_bar = "se")
  pdf(file=paste0(pathoutput, "/", modelname[semID], "SEM_DEG.pdf"), paper='a4r', width=200)
  plt <- plot(my_fitted_DAG, las = 2 )
  print(plt)
  dev.off()
}


CVconverged <- c(matrixCSV[1,], matrixCSV[2,], matrixCSVtot[1,], matrixCSVtot[2,])
CVconv <- data.frame(CVconverged)
write.csv(CVconv, paste0(pathoutput, "/CSVconverged", ".csv"))

save.image(file=paste0(pathoutput, "/", "/imageworkspace_phylosem.RData"))


 
# #*
# #* CROSS-VALIDATION RUN
# #*
# 
# 
# #RUN cross validation
# 
load(file=paste0(pathoutput, "/", "/imageworkspace_phylosem.RData"))
source(paste0("01-Dataset_construction/Scripts/00-Functions_for_phylosem.R"))
# 
# # cross-validation on all traits at the same time
sampletot           <- psem(semID = semID, trait=colnames(dataset_traits), nameCV = nameCV[2]) #CV on all
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))
# 
# # cross-validation on c_m only : genus by genus
samplec_mspe  <- psem(semID = semID, trait="c_m",nameCV =  nameCV[1])
save.image(file=paste0(pathoutput, "/", "imageworkspace_phylosem.RData"))


