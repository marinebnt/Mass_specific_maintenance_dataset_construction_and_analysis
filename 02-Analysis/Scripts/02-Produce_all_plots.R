



# loop over all phylosem outputs (tested different SEMs)
paths_dir <-  paste0("01-Dataset_construction/Outputs/phylosem")
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
  path_func <- "02-Analysis/Scripts/functions/"
  functions <- list.files(path_func)
  load(paste0("01-Dataset_construction/Outputs/phylosem/", model, "/imageworkspaceEND.RData")) #data needed for cross validation
  
  
  ######## What input and output folder ? 
  # traits = c("Age.mat", "Age.max", "Mortality", "K", "Trophic.lvl", "Habitat")
  traits = c("Age.mat", "Age.max", "Mortality", "K") # time related traits
  OUTPUT = paste0("Outputs_time/", model) # time-related outputs
  OUTPUT_phylo = paste0("Outputs/phylosem/", model)
  path_phylosem_out <- paste0(getwd(), "/02-Analysis/", OUTPUT)
  path_plots <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")
  path_CV <- paste0(getwd(), "/01-Dataset_construction/", OUTPUT_phylo)
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
    cat("\nUsing previously produced Archetypal analysis outputs.\n")
  }
  
  # load inputs 
  load(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_FOR_ANALYSIS.RData"))
  source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
  
  # correction, because loaded old names
  path_output_genus <- stringr::str_replace(path_output_genus, "01-Simulations", replacement = "01-Dataset_construction")
  path_phylosem_out <- stringr::str_replace(path_phylosem_out, "01-Simulations", replacement = "01-Dataset_construction")
  
  # loop the plots
  for (func in functions){
    if (func == c("plot_AA_elbow_criterion.R")){next} # this is really long, run only if needed
    path_CV <- paste0(getwd(), "/01-Dataset_construction/", OUTPUT_phylo)
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
