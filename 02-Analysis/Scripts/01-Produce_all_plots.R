#* 17/08/25
#* Beneat marine
#* This script loops over phylosem outputs folders. 
#* It creates one output folder per SEM. 
#* If the SEM is the chosen one, it runs all the analyses related, including Archetypal analysis. 
#* If previously produced Archetypal analysis missing : longer to run, and results might change slightly from publication. 
#* If the SEM is not the chosen one, it produces the cross-validation plot only, with the related outliers. 



# loop over all phylosem outputs (tested different SEMs)
supp      <- "_output"
paths_dir <-  paste0("01-Dataset_construction/Outputs/phylosem", supp)
list_dir  <- list.dirs(paths_dir, recursive = F, full.names = F)
list_dir  <- list_dir[-grep("jack", list_dir)]
kmax=3 

# functions to run to produce plots :
path_func <- "02-Analysis/Scripts/functions/"
functions <- list.files( "02-Analysis/Scripts/functions/")
functions <- functions[grep("plot", functions)]
functions <- functions[-c(4,8, 11)] # plots not used anymore : used to compare teleo and elasmo

for (dir in list_dir){ # loop over all phylosem outputs to analyse
  
  # What is the chosen model out of the tested SEM models ? 
  chosenSEM = c("TLstdmeca")
  model <- dir
  
  ######### SEM MODEL & associated plots names
  cat("\n", dir)
  if(model == "stdevol00"){
    semname = "Evolutionary (a)"
  }
  if(model == "stdevol0"){
    semname = "Evolutionary (b)"
  }
  if(model == "stdLoo"){
    semname = "Mechanistic 2(a)"
  }
  if(model == "stdmeca"){
    semname = "Mechanistic 1(a)"
  }
  if(model == "TLstdLoo"){
    semname = "Mechanistic 2(b)"
  }
  if(model == "TLstdmeca"){
    semname = "Mechanistic 1(b)"
  }
  if(model == "TLWoo"){
    semname = "Mechanistic 3(b)"
  }
  if(model == "TLWooWoo"){
    semname = "Mechanistic 4(b)"
  }
  if(model == "Woo2"){
    semname = "Mechanistic 3(a)"
  }
  if(model == "WooWoo2"){
    semname = "Mechanistic 4(a)"
  }
  
  
  # load R functions producing the plots
  path_image <- paste0("01-Dataset_construction/Outputs/phylosem", supp, "/", model, "/imageworkspace_phylosem.RData")
  load(path_image) #data needed for cross validation  #imageworkspaceEND.RData
  
  
  ######## What input and output folder ? 
  traits           <- c("Age.mat", "Age.max", "Mortality", "K", "Fecundity") # time related traits
  OUTPUT           <- paste0("Outputs_phylosem/", model)
  OUTPUT_phylo     <- paste0("Outputs/phylosem", supp, "/", model, "/phylosem")
  path_phylosem_out<- paste0("01-Dataset_construction/", OUTPUT_phylo) #, "/phylosem")
  path_analysis_out<- paste0("02-Analysis/", OUTPUT)
  path_plots       <- paste0("02-Analysis/", OUTPUT,"/plots/")
  path_phylosem_CV <- paste0("01-Dataset_construction/Outputs/phylosem", supp, "/", model, "/CrossValidation")
  path_phylosem_jack<- paste0("01-Dataset_construction/Outputs/phylosem", supp, "/", model)
  path_CV          <- paste0("02-Analysis/", OUTPUT,"/CrossValidation_error/")
  pathoutput_CV    <- paste0(path_analysis_out,"/plots")
  pathoutput_jack  <- paste0(pathoutput_CV,"/jackknife")
  source(paste0("02-Analysis/Scripts/00-Functions_for_analysis.R"))
  ########
  dir.create(path = path_CV, recursive =T)
  dir.create(path = pathoutput_jack, recursive =T)
  dir.create(paste0(pathoutput_CV, "/plot_CrossValidation"), recursive =T)
  

  
  ##### load input functions ######
  source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
  
  
  
  
  #### Prepare dataset #####
  IS_LOG_M=T
  # dataset output phylosem
  dataphylo <- read.csv(paste0(path_phylosem_out, "/output_SEMpsemFINALtot.csv"), row.names = c(2))
  # dataset input phylosem
  datagenus <- dataset
  # make both datasets clean
  datagenus$na <- is.na(datagenus$c_m)
  dataphylo$Species <- stringr::str_replace(rownames(dataphylo), "_", " ")
  dataphylo <- addhabitatcolumn(dataphylo)
  # Beverton-holt
  dataphylo$Lm <- ratioTOtrait(dataphylo$Lm, dataphylo$Loo, IS_OPERATOR_Loo=TRUE, IS_LOG_M=IS_LOG_M)
  dataphylo$tm <- ratioTOtrait(dataphylo$tm, dataphylo$M, IS_OPERATOR_Loo=FALSE, IS_LOG_M=IS_LOG_M) # 1/ M =~ longevity
  if(length(which(dataphylo$Lm>dataphylo$Loo))>0) {dataphylo <- dataphylo[-which(dataphylo$Lm>dataphylo$Loo),]}
  if(length(which(dataphylo$tm>dataphylo$tmax))>0) {dataphylo <- dataphylo[-which(dataphylo$tm>dataphylo$tmax),]}
  #Standardize
  dataphylo[, -which(colnames(dataphylo) %in% c("Species", "X", "label", "hhabtot", "osmose", "na"))] <-
    decostand(dataphylo[, -which(colnames(dataphylo) %in% c("Species", "X", "label", "hhabtot", "osmose", "na"))], "standardize")
  # merge dataset with species/genus/family names and phylosem output
  dataphylo <- left_join(dataphylo[,-which(colnames(dataphylo) == 'X')], datagenus[, c("Species", "na", "Genus", "Class", "Order", "Family", "SpecCode")])
  # prettier colnames
  colnames(dataphylo)[which(colnames(dataphylo) == "hhabtot")]<- "Habitat"
  colnames(dataphylo)[which(colnames(dataphylo) == "fecundity")]<- "Fecundity"
  colnames(dataphylo)[which(colnames(dataphylo) == "TLDiet")]<- "Trophic.lvl"
  colnames(dataphylo)[which(colnames(dataphylo) == "tmax")]<- "Age.max"
  colnames(dataphylo)[which(colnames(dataphylo) == "tm")]<- "Age.mat"
  colnames(dataphylo)[which(colnames(dataphylo) == "M")]<- "Mortality"
  colnames(dataphylo)[which(colnames(dataphylo) == "Woo")]<- "Weight.Inf"
  dataplot  <- data.frame("K"=dataphylo$K, "M"=dataphylo$Mortality, "Temperature"=dataphylo$Temperature,
                          dataphylo[, -which(colnames(dataphylo) %in%
                                               c("Class", "Order", "Family", "Genus", "Species", "SpecCodeode", "c_m", "T"))])
  dataacp   <- data.frame(dataphylo)

  
  ###### Run archetype analysis only for the selected SEM model #######
  if (model == chosenSEM){
    if (!file.exists(paste0(path_analysis_out ,"/IMAGE_AA_CONSTRUCTED.RData"))){
      cat("\nRunning Archetypal analysis.\n")
      source("02-Analysis/Scripts/functions/01-Archetypal_analysis.R")
    }
    else{
      cat("\nUsing previously produced Archetypal analysis outputs.\n")
      load(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_CONSTRUCTED.RData"))
    }
  }


  ##### loop to produce plots #####
  for (func in functions){
    if (func == c("plot_AA_elbow_criterion.R")){next} # this is really long, run only if needed
    if (model != chosenSEM && func == c("plot_cross_validations.R")){
      cat("=>", func, "\n")
      source(paste0(path_func, "/", func), echo = FALSE, print.eval = FALSE, verbose = FALSE)
    }
    if (model == c("TLstdmeca")){
      cat("=>", func, "\n")
      source(paste0(path_func, "/", func), echo = FALSE, print.eval = FALSE, verbose = FALSE)
    }
    else {
      next
    }
  }
}

