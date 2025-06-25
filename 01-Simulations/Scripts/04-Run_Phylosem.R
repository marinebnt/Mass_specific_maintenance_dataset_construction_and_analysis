###################### use dataext with phylosem : with all data from fishbase #######################
# needed  : output of the create_dataext_forphylosem_genus.R file, called : output_tot/dataext_traits_phylosem.csv
################################################################################################

# Write paths and load dependencies

path <- paste0("01-Simulations/Outputs/dataset_creation_output/dataset_for_phylosem")
pathoutput <- paste0( "01-Simulations/Outputs/phylosem_output")
source(paste0("01-Simulations/Scripts/00-Functions_for_phylosem.R"))


# load dataset
dataset        <- read.csv(paste0(path, "/output_tot_stdmorpho/dataset_phylosemLOG.csv"))
rownames(dataset) <- dataset$X
dataset_traits <-  read.csv(paste0(path, "/output_tot_stdmorpho/dataset_traits_phylosemLOG.csv"))
rownames(dataset_traits) <- dataset_traits$X
if(length(grepl("_", dataset$X))==0){
  rownames(dataset) <- stringr::str_replace(pattern = " ", replacement = "_", dataset$X)
  rownames(dataset_traits) <- stringr::str_replace(pattern = " ", replacement = "_", dataset_traits$X)}
  
# remove useless traits
dataset        <- dataset[, -which(grepl("X", colnames(dataset)))]
dataset_traits <- dataset_traits[, -which(grepl("X", colnames(dataset_traits)))]
dataset_traits <- dataset_traits[,-which(colnames(dataset_traits) %in% c("eps_m", "AspectRatio",
                                                                         "waterbrack", "watersaltwater", "lengthOffspring"))]

# replacing maturation traits with invariant traits
# WITH THE RATIO (INVARIANTS beverton & holt)
boxplot(exp(dataset_traits$Lm)/exp(dataset_traits$Loo))
dataset_traits$Lm <- dataset_traits$Lm-dataset_traits$Loo
boxplot(exp(dataset_traits$tm)*(dataset_traits$M))
dataset_traits$tm <- log(exp(dataset_traits$tm)*dataset_traits$M)

# dataset names
dataset -> dataext
dataset_traits -> dataext_traits


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
nb         <- dim(dataext)[1]
frm        = ~SuperClass/Class/Order/Family/Genus/Species
phylo      <- c()
phylo      = as.data.frame(dataext[which(names(dataext) %in% c("SuperClass", "Class", "Order", "Family", "Genus", "Species"))], stringsAsFactors = TRUE)
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
#trait = c("c_m") # choose the trait you want to cross validate ****
nbCV <- 10
semID=c(1)

listsem <- list(TLstdmeca)
max     <- length(rownames(dataext))
IDNA_c_m<- which(!is.na(dataext_traits$c_m))
IDtot   <- length(dataext_traits$c_m)
max25   <- length(IDNA_c_m)/4
max33   <- length(IDNA_c_m)/nbCV
max25tot   <- IDtot/4
max33tot   <- IDtot/nbCV
modelname <- c("TLstdmeca")
modellist <-list(TLstdmeca)
matrixCSV <- matrix("", 2,10)
matrixCSVtot <- matrix("", 2,1)
nameCV <- c("c_m", "all", "spe")

###########*
###########*


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
#                "Max age", "K", "Infinite weight", "Mortality", "Trophic level", "Infinite length", "Fecundity",  "Peduncle depth", 
#                "Jaw length", "Body depth", "Body width", "Temperature")
# 
# 
# # 
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
# # # #*
# # # #* CROSS-VALIDATION CHECK & PLOTS
# # # #*
# 
# 
# 
# ncol=5
# nrow=4
# plot_checkphylosemdata(semID, trait=c("c_m", "K", "M", "Loo", "habitatpelagic", "habitatbenthopelagic", "habitatdemersal", "Woo"), name=nameCV[3],
#                        sample=list(samplec_mc_m$sample, sampleKspe$sample, sampleMspe$sample, sampleLoospe$sample,
#                                    samplehabitatpelagicspe$sample, samplehabitatbenthopelagicspe$sample, samplehabitatdemersalspe$sample,
#                                    sampleWoospe$sample),
#                        maxCV=list(samplec_mc_m$maxCV, sampleMspe$maxCV, sampleKspe$maxCV, sampleLoospe$maxCV,
#                                   samplehabitatpelagicspe$maxCV, samplehabitatbenthopelagicspe$maxCV, samplehabitatdemersalspe$maxCV,
#                                   sampleWoospe$maxCV), names_var)



#* 
#* COMPARE PHYLOGENETIC METHODS
#* 
# 

Grid = expand.grid( "OU" = c(FALSE,TRUE),
                    "lambda" = c(FALSE,TRUE),
                    "kappa" = c(FALSE,TRUE) )
psem_grid = NULL
for( i in 1:4){ #nrow(Grid)
  data_CV <- dataext_traits
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

write.csv2(as.data.frame(Grid), paste0(pathoutput, "/phylogeneticmethods.csv"))
