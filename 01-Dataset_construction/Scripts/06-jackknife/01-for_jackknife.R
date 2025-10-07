###################### use dataext with phylosem : with all data from fishbase #######################
# needed  : output of the create_dataext_forphylosem_genus.R file, called : output_tot/dataext_traits_phylosem.csv
################################################################################################

# Write paths and load dependencies

path <- paste0("01-Dataset_construction/Outputs/dataset_creation_output", supp, "/dataset_for_phylosem_NOUNITCV")
pathoutput <- paste0("01-Dataset_construction/Outputs/dataset_creation_output", supp, "/jackknife")
source(paste0("01-Dataset_construction/Scripts/00-Functions_for_phylosem.R"))
dir.create(pathoutput)

# load dataset
dataset        <- read.csv(paste0(path, "/output_tot_stdmorpho/eladataset_phylosem.csv"), row.names = c(1))
dataset_traits <-  read.csv(paste0(path, "/output_tot_stdmorpho/eladataset_traits.csv"), row.names = c(1))
rownames(dataset) <- stringr::str_replace(rownames(dataset), pattern = " ", replacement = "_")
rownames(dataset_traits) <- stringr::str_replace(rownames(dataset_traits), pattern = " ", replacement = "_")

# dataset names
dataset -> dataext
dataset_traits -> dataext_traits

# prepare_jacklist.R (run interactively on login node)
dat <- dataset[which(!is.na(dataset$c_m)),]  # or read.csv/other
group_col <- "Genus"  # or "SpecCode"
groups <- unique(as.character(dat[[group_col]]))
length(groups)  # number of jackknife replicates needed

saveRDS(groups, paste0(pathoutput, "/jack_list.rds")) # list of group ids to leave out
a <- readRDS(paste0(pathoutput, "/jack_list.rds"))
length(a)
