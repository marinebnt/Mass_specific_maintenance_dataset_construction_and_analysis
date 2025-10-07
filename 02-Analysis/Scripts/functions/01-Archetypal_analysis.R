# beneat marine 
# 14/10/24 
############################
#* 
#* THIS SCRIPT GOAL IS TO RUN PCAs AND phylogenetic PCAs AND Archetypal Analysis
#* ITS OUTPUTS ARE USED AS INPUTS FOR THE PLOT CONSTRUCTIONS OF THE FOLLOWING SCRIPTS
#* 




############################## RUN Archetypal analysis ##########################

AAtot <- runAA(dataplot, traits, kmax)


############################## SAVE outputs ##########################

save.image(paste0(path_analysis_out, "/IMAGE_AA_FOR_ANALYSIS.RData"))

save(AAtot, dataphylo, dataacp, dataplot, datagenus, 
          file=paste0(path_analysis_out, "/IMAGE_AA_CONSTRUCTED.RData"))

