# beneat marine 
# 14/10/24 
# Cross validation plots of the data inference



semID=1


p <- plot_checkphylosemdata(semID = semID, trait="c_m", name=nameCV[2],
                       sample=list(sampletot$sample),
                       maxCV=list(sampletot$maxCV), names_var="log(RMR0)")
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
                       maxCV=list(sampletot$maxCV), names_var="log(K)")
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
                       maxCV=list(sampletot$maxCV), names_var="log(Mortality)")
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

first = (p)
second =  (p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14 + plot_layout(ncol = 4, nrow=4, guides = "collect", axis_titles = "collect"))
third =  ((ph1+ ph2 + ph3 + plot_layout(ncol=2, nrow=2, guides = "collect", axis_titles = "collect")))

final<- (first)+
    (second + plot_layout(tag_level = "new"))+
    (third + plot_layout(tag_level = "new"))+
  plot_layout(design=design)+
  plot_annotation(title = paste0("Cross-validation with ", semname, " SEM"),            
                  tag_levels = c("A", "1"))


grDevices::pdf(file=paste0(path_plots, "/plot_CrossValidation/", "plotarrangedCV_clean.pdf"), height=9, width=15.37)
print(final)
dev.off()









