What folder name corresponds to what model ? 
if you have a doubt, you can check the files SEMSEM_DEG.pdf and semmodel_SEM.csv in each folder

stdevol00   > Evolutionary
stdevol0    > Evolutionary with trophic level imapacted by all morphometric traits
stdmeca     > Mechanistic 1
TLstdmeca   > Mechanistic 1 with trophic level imapacted by all morphometric traits
stdLoo      > Mechanistic 2
TLstdLoo    > Mechanistic 2 with trophic level imapacted by all morphometric traits
Woo2        > Mechanistic 3
TLWoo       > Mechanistic 3 with trophic level imapacted by all morphometric traits
WooWoo2     > Mechanistic 4
TLWooWoo    > Mechanistic 4 with trophic level imapacted by all morphometric traits



######################
# File by file description
######################

>The Cross-validation outputs are numbered and each correspond to a trait the cross-validation is running on, or to the complete dateset ('allTOT'). 
ex : output1_SEMpsemFINALspehabitatbenthopelagic.csv, first round of the cross-validation, where the NAs are randomly replacing known habitatbenthopelagic trait binary values. These NAs are then inferred.  

>The comparison between observed and inferred trait values are in the csv files called "dataseterrorSEM_something_.csv". 

>The complete infered dataset is in 'output_SEMpsemFINALtot.csv'

>Preliminary plots are stored in the folders 'plot' and 'boxplot'. 

>coefnostd_std_SEM.csv and semmodel_SEM.csv are giving complementary information about the path coefficient values. SEMSEM_DEG.pdf is a rough plot showing the SEM. 

>imageworkspaceEND.RData gives many info about the whole process. For instance, for the cross-validation process, it informs on what trait values were set as NAs and in what order. 