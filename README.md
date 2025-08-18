### 2025 - Mass-specific dataset construction and data analysis
## Overview
This repository gathers the codes and data we used in the study. First is the dataset construction, and second is the data analysis. Codes and data generated and transformed in the study are separated in their respective folders.

## Authors & Contact details:
Marine BENEAT (Marine.Beneat@ifremer.fr),  
Alaia MORELL (amorell@ifremer.fr),  
Fabien MOULLEC (fabien.moullec@umontpellier.fr),  
Nicolas BARRIER (nicolas.barrier@ird.fr),  
Yunne SHIN (yunne-jai.shin@ird.fr),  
Bruno ERNANDE (Bruno.Ernande@ifremer.fr).  
Responsible of data collection, codes development and online upload: Marine BENEAT.

## Cite the code:
DOI

## Brief summary of the study


## Layout
The repository is split into three main directories: The first one is about the dataset construction, the second about data analysis, and the third one is the .csv file with the 18,129 species trait values estimated using the inference process. 
The folders 01- and 02- contain an Input folder (for 02-: the input data is the output of the folder 01-), a Script folder, and an Output folder.

The dataset construction folder (01-) has two steps. First is the collection of the temperature- and mass-specific respiration data. Second is the input of the mass-specific data and other traits values into the inference process. The inference process is used as input of the data analysis. 

The analysis folder (02-) both produces the cross-validation plots for all the SEM tested, and produces the analysis plots for the selected SEM model. 

The dataset folder (03-) is the raw output from the inference process using the selected SEM (mechanistic 1 (b)). Some of the trait values are log-transformed according to the method described in the corresponding paper. 
