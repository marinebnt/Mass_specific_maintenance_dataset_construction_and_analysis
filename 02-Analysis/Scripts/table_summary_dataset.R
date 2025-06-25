# beneat marine 
# 14/10/24 
# To build a data frame of the summary of the observed data, before inference


# Convert summary to a dataframe (just an example if it's stored)
source(paste0("02-Analysis/Scripts/00-Functions_for_analysis.R"))
path_output_genus <- paste0("01-Simulations/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV/output_tot_stdmorpho")
datagenus <- read.csv(paste0(path_output_genus, "/dataset_phylosem.csv"))
datagenus_t <- read.csv(paste0(path_output_genus, "/dataset_traits_phylosem.csv"))

summary_data <- summary(datagenus)


# Create a nice clean table using kable
summary_table <- kable(summary_data, format = "html", table.attr = "class='table table-bordered'", caption = "Summary of Dataset") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                full_width = F, 
                position = "left") %>%
  scroll_box(height = "500px")


# Print the table
summary_table




# Further analysis to spot differences elasmobranchii and teleostei

# How many NAs are there ? 

elasmo <- datagenus_t[which(datagenus$Class=="Elasmobranchii"),]
teleo <- datagenus_t[which(datagenus$Class=="Teleostei"),]

whichX <- which(colnames(elasmo)==c("X"))
elasmo_na <- elasmo[,-whichX]
elasmo_na[which(!is.na(elasmo_na), arr.ind = T)] <- 0
elasmo_na[which(is.na(elasmo_na), arr.ind = T)] <- 1
telasmo <- table(rowSums(elasmo_na))
nelasmo <- nrow(elasmo)
p_telasmo <- telasmo*100/nelasmo

whichX <- which(colnames(teleo)==c("X"))
teleo_na <- teleo[,-whichX]
teleo_na[which(!is.na(teleo_na), arr.ind = T)] <- 0
teleo_na[which(is.na(teleo_na), arr.ind = T)] <- 1
tteleo <- table(rowSums(teleo_na))
nteleo <- nrow(teleo)
p_tteleo <- tteleo*100/nteleo


elasmo <- dataphylo[which(dataphylo$Class=="Elasmobranchii"),]
teleo <- dataphylo[which(dataphylo$Class=="Teleostei"),]

range(elasmo$c_m)
range(teleo$c_m)
