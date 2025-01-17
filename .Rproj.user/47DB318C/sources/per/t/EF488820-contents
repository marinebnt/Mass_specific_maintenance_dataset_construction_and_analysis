# beneat marine 
# 14/10/24 
# Cross validation plots of the data inference


setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/tests/repository_for_zenodo")
load(paste0(getwd(), "/01-Simulations/Outputs/phylosem_output/imageworkspaceEND.RData")) #data needed for cross validation
path_plots <- paste0(getwd(), "/02-Analysis/Outputs")
pathoutput_CV <- path_plots
path_CV <- paste0(getwd(), "/01-Simulations/Outputs/phylosem_output")
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))

semID=1

ncol = 3
nrow = 3
# plot_checkphylosemdata(semID, trait=c("c_m", "K", "M", "Loo", "habitatpelagic", "habitatbenthopelagic", "habitatdemersal", "Woo"), name=nameCV[2], 
                       # sample=list(sampletot$sample), 
                       # maxCV=list(sampletot$maxCV))
ncol = 6
nrow = 3
names_var <- c("log(MSM)", "log(Maturation age)", "log(Maturation length)", 
                             "log(Max age)", "K",
               "Benthopelagic", "log(Infinity weight)", "M", "Trophic level", "log(Infinitylength)", "log(Fecundity)", "Demersal",  
                             "log(Peduncle depth)", 
                             "log(Jaw length)", "log(Body depth)", "log(Body width)", "log(Temperature)", "Pelagic ")
traits = c("c_m", "tm", "Lm", "tmax", "K", "habitatbenthopelagic", 
           "Woo", "M", "TLDiet", "Loo", "fecundity", "habitatdemersal",                         
           "Min_caudalpeduncle_depth", "Lower_jaw_length", "Max_body_depth", "Max_body_width","Temperature", "habitatpelagic") 
plot_checkphylosemdata(semID, trait=traits, name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var, 
                       plot_layout)





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

grDevices::pdf(file=paste0(pathoutput_CV, "/plot/", modelname[sem], name, "plotarrangedCV_clean.pdf"), height=9, width=15.37)
print(final)
dev.off()



plot_checkphylosemdata(semID, trait=traits, name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var)semID=1
plot_checkphylosemdata(semID, trait=c("c_m"), name=nameCV[3], 
                       sample=list(samplec_mspe$sample), 
                       maxCV=list(samplec_mspe$maxCV), names_var=c("log(MSM)"))

plot_checkphylosemdata(semID, trait=c("Max_body_depth"), name=nameCV[2], 
                       sample=list(sampletot$sample), 
                       maxCV=list(sampletot$maxCV), names_var=c("log(age mat)"))


semmodel="TLstdmeca"
semID=1
plot_checkphylosemdata(semID, trait=c("habitatpelagic"), name=nameCV[2], 
                       sample=list(samplec_mspe$sample), 
                       maxCV=list(samplec_mspe$maxCV), names_var=c("Pelagic"))



