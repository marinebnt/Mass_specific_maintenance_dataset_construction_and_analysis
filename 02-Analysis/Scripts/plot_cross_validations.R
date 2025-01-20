# beneat marine 
# 14/10/24 
# Cross validation plots of the data inference


load(paste0(getwd(), "/01-Simulations/Outputs/phylosem_output/imageworkspaceEND.RData")) #data needed for cross validation
path_plots <- paste0(getwd(), "/02-Analysis/Outputs")
pathoutput_CV <- path_plots
path_CV <- paste0(getwd(), "/01-Simulations/Outputs/phylosem_output")
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))


semID=1


p <- plot_checkphylosemdata(semID, trait="c_m", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(MSM)")
p2 <- plot_checkphylosemdata(semID, trait="tm", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Maturation age)")
p3 <- plot_checkphylosemdata(semID, trait="Lm", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Maturation length)")
p4 <- plot_checkphylosemdata(semID, trait="tmax", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Max age)")
p5 <- plot_checkphylosemdata(semID, trait="K", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="K")
ph1 <- plot_checkphylosemdata(semID, trait="habitatbenthopelagic", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="Benthopelagic")
p6 <- plot_checkphylosemdata(semID, trait="Woo", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Infinity weight)")
p1 <- plot_checkphylosemdata(semID, trait="Loo", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Infinity length)")
p7 <- plot_checkphylosemdata(semID, trait="M", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="Mortality")
p8 <- plot_checkphylosemdata(semID, trait="TLDiet", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="Trophic level")
p9 <- plot_checkphylosemdata(semID, trait="fecundity", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Fecundity)")
ph2 <- plot_checkphylosemdata(semID, trait="habitatdemersal", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="Demersal")
p10 <- plot_checkphylosemdata(semID, trait="Min_caudalpeduncle_depth", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Peduncle depth)")
p11 <- plot_checkphylosemdata(semID, trait="Lower_jaw_length", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Jaw length)")
p12 <- plot_checkphylosemdata(semID, trait="Max_body_depth", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Body depth)")
p13 <- plot_checkphylosemdata(semID, trait="Max_body_width", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Body width)")
p14 <- plot_checkphylosemdata(semID, trait="Temperature", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(Temperature)")
ph3<- plot_checkphylosemdata(semID, trait="habitatpelagic", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="Pelagic ")

library(patchwork)


design <- "
  112222
  112222
  332222
"

final<-p+
 (p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14 + plot_layout(ncol = 4, nrow=4, guides = "collect")) + 
 ((ph1+ ph2 + ph3 + plot_layout(ncol=2, nrow=2, guides = "collect"))) + 
  plot_layout(guides = "collect", design=design)   # guide_area() + 

grDevices::pdf(file=paste0(pathoutput_CV, "/plot_CrossValidation/", "plotarrangedCV_clean.pdf"), height=9, width=15.37)
print(final)
dev.off()


