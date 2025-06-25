



# loop over all phylosem outputs (tested different SEMs)
paths_dir <-  paste0("01-Simulations/Outputs/phylosem")
list_dir <- list.dirs(paths_dir, recursive = F, full.names = F)

for (dir in list_dir[6]){
  
  ######### SEM MODEL
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
  path_func <- "C:/Users/mbeneat/Documents/osmose/updated_parameters_osmose-med/Mass_specific_maintenance_dataset_building/02-Analysis/Scripts/functions"
  functions <- list.files(path_func)
  load(paste0("01-Simulations/Outputs/phylosem/", model, "/imageworkspaceEND.RData")) #data needed for cross validation
  
  
  ######## What input and output folder ? 
  # traits = c("Age.mat", "Age.max", "Mortality", "K", "Trophic.lvl", "Habitat")
  traits = c("Age.mat", "Age.max", "Mortality", "K")
  OUTPUT = paste0("Outputs_time/", model)
  OUTPUT_phylo = paste0("Outputs/phylosem/", model)
  path_phylosem_out <- paste0(getwd(), "/02-Analysis/", OUTPUT)
  path_plots <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")
  path_CV <- paste0(getwd(), "/01-Simulations/", OUTPUT_phylo)
  pathoutput_CV <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")
  source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
  ########
  dir.create(path = path_CV, recursive =T)
  dir.create(paste0(pathoutput_CV, "/plot_CrossValidation"), recursive =T)
  
  
  
  if (!file.exists(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_FOR_ANALYSIS.RData"))){
    cat("\nRunning Archetypal analysis.\n")
    source("02-Analysis/Scripts/01-Prepare_data_for_analysis.R")
  }
  else{
    cat("\nUsing previously produced Archetypal outputs.\n")
  }
  
  # load inputs 
  load(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_FOR_ANALYSIS.RData"))
  source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
  
  # loop the plots
  for (func in functions[3]){
    if (func == c("plot_AA_elbow_criterion.R")){next} # this is really long, run only if needed
    path_CV <- paste0(getwd(), "/01-Simulations/", OUTPUT_phylo)
    if (model != c("TLstdmeca") && func == c("plot_cross_validations.R")){
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
