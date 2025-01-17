# This is to fill in the content of the table 



# Convert summary to a dataframe (just an example if it's stored)
library(knitr)
library(kableExtra)

source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
path_output_genus <- paste0(getwd(), "/01-Simulations/Outputs/dataset_creation_output/dataset_for_phylosem/output_tot_stdmorpho")
datagenus <- read.csv(paste0(path_output_genus, "/dataset_phylosemLOG.csv"))
summary_data <- summary(datagenus)


# Create a nice clean table using kable
summary_table <- kable(summary_data, format = "html", table.attr = "class='table table-bordered'", caption = "Summary of Dataset") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                full_width = F, 
                position = "left") %>%
  scroll_box(height = "500px")


# Print the table
summary_table

