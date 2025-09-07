#* 17/08/25
#* Beneat marine
#* This script loops over phylosem outputs folders. 
#* It creates one output folder per SEM. 
#* If the SEM is the chosen one, it runs all the analyses related, including Archetypal analysis. 
#* If previously produced Archetypal analysis missing : longer to run, and results might change slightly from publication. 
#* If the SEM is not the chosen one, it produces the cross-validation plot only, with the related outliers. 



# loop over all phylosem outputs (tested different SEMs)
paths_dir <-  paste0("01-Dataset_construction/Outputs/phylosem")
list_dir <- list.dirs(paths_dir, recursive = F, full.names = F)

for (dir in list_dir){
  
  # What is the chosen model out of the tested SEM models ? 
  chosenSEM = c("TLstdmeca")
  
  ######### SEM MODEL & associated plots names
  model <- dir
  cat("\n", model)
  if(model == "stdevol0"){
    semname = "Evolutionary (a)"
  }
  if(model == "stdevol00"){
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
  path_func <- "02-Analysis/Scripts/functions/"
  functions <- list.files(path_func)
  load(paste0("01-Dataset_construction/Outputs/phylosem/", model, "/imageworkspace_phylosem.RData")) #data needed for cross validation  #imageworkspaceEND.RData
  
  
  ######## What input and output folder ? 
  traits = c("Age.mat", "Age.max", "Mortality", "K") # time related traits
  OUTPUT = paste0("Outputs/", model)
  OUTPUT_phylo = paste0("Outputs/phylosem/", model, "/phylosem")
  path_phylosem_out <- paste0(getwd(), "/01-Dataset_construction/", OUTPUT_phylo)
  path_plots <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")
  path_CV <- paste0(getwd(), "/01-Dataset_construction/", OUTPUT_phylo)
  pathoutput_CV <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")
  source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
  ########
  dir.create(path = path_CV, recursive =T)
  dir.create(paste0(pathoutput_CV, "/plot_CrossValidation"), recursive =T)
  
  # Run archetype analysis only for the selected SEM model
  if (model == chosenSEM){
    if (!file.exists(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_FOR_ANALYSIS.RData"))){
      cat("\nRunning Archetypal analysis.\n")
      source("02-Analysis/Scripts/01-Prepare_data_for_analysis.R")
    }
    else{
      cat("\nUsing previously produced Archetypal analysis outputs.\n")
    }
  }
  
  # load inputs 
  load(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_CONSTRUCTED.RData"))
  source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
  
  # correction, because loaded old names
  # path_output_genus <- stringr::str_replace(path_output_genus, "01-Simulations", replacement = "01-Dataset_construction")
  # path_phylosem_out <- stringr::str_replace(path_phylosem_out, "01-Simulations", replacement = "01-Dataset_construction")
  # OUTPUT = paste0("Outputs/", model) # time-related outputs
  # path_phylosem_out <- paste0(getwd(), "/02-Analysis/", OUTPUT)
  # path_plots <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")
  
  ## loop the plots
  for (func in functions){
    if (func == c("plot_AA_elbow_criterion.R")){next} # this is really long, run only if needed
    path_CV <- paste0(getwd(), "/01-Dataset_construction/", OUTPUT_phylo)
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
  # }
}
