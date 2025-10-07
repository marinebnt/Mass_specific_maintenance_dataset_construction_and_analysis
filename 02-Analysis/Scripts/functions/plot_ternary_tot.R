# beneat marine 
# 14/10/24 
# Ternary plots of the archetypal analysis



# veccol veccnames
AAAtot <- AAtot
AAAtot$alphas <- as.data.frame(AAAtot$alphas)[,order((as.data.frame(AAAtot$archetypes)$K), decreasing = F)]
AAAtot$archetypes <- as.data.frame(AAAtot$archetypes)[order((as.data.frame(AAAtot$archetypes)$K), decreasing = F),]

colAA <- c()
labelAA <- c()
pchvec <- c()
iopp <- which.max(as.data.frame(AAAtot$archetypes)$K) 
colAA[iopp] = "tomato"
labelAA[iopp] = "Fast"
pchvec[iopp] = c(23)
ieq <- which.min(as.data.frame(AAAtot$archetypes)$K) 
colAA[ieq] = "royalblue"
labelAA[ieq] = "Slow"
pchvec[ieq] = c(22)
iper <- c(1:3)[-c(ieq, iopp)]
colAA[iper] = "darkgreen"
labelAA[iper] ="Intermediate"
pchvec[iper] = c(21)

i_elasmo = 30
i_opp = 21

vecAA_eq = c("a","b","c")


# ternary plot with the gradient of mass specific maintenance values
beta_iv1_tot = dataplot[, which(colnames(dataplot) %in%  traits)]
t(parameters(AAAtot))     ## This shows the trait values of the selected archetypes
veccol_tot_barplot <- colAA
barplot(AAAtot, beta_iv1_tot, percentiles=TRUE, col.atypes = veccol_tot_barplot)
barplot(AAAtot, beta_iv1_tot, percentiles=TRUE)  # same but as a barplot
pcplot(AAAtot, beta_iv1_tot, atypes.col = veccol_tot_barplot, atypes.lty =c(1,2,3)) # same but as a plot showing each value
plot1 <- plotggtern(dataphylo=dataphylo, AA=AAAtot, vecAA_eq=vecAA_eq, bins=6)+
  labs(
    x = labelAA[1],       # Custom name for x-axis
    y = labelAA[2],       # Custom name for y-axis
    z = labelAA[3],        # Custom name for z-axis
  )+
  theme(
    tern.axis.title.T = element_text(size = 12, angle = 0, hjust = 0.5), # Top axis
    tern.axis.title.L = element_text(size = 12, angle = -60, hjust = 0.5), # Left axis
    tern.axis.title.R = element_text(size = 12, angle = -60, hjust = 1.25, vjust = -0.5), # Right axis
    plot.title = element_text(size=15, face="bold")
  )+
  labs(title="B")

plot1

pdf(paste0(path_plots, "/ternary_plot_tot_time.pdf"))
print(plot1)
dev.off()






# ternary plot with example species
cmstd <- dataphylo$c_m # standardized c_m values
alpha <- AAAtot$alphas # position of each species in the ternary plot. Same order as the original dataset (dataphylo). 
                       # To know which archetype is which, have a look at the order of the AAAtot$archetype
df <- cbind(alpha, cmstd)
df_mean <- aggregate(df[, 1:3], list(dataphylo$Order), mean)
# df_mean the family with the highest average value on each archetype is the most representative one of the archetype
maxAA1 <- df_mean$Group.1[which.max(df_mean$V1)] # unknown family, I prefere to pick another one with df_mean$Group.1[which(df_mean$V1>0.7)]
maxAA1 <-df_mean$Group.1[which(df_mean$V1>0.6)][1] # unknown family, I prefere to pick another one with df_mean$Group.1[which(df_mean$V1>0.7)]
# maxAA1 <- "Psettodidae" # pleuronectiformes
maxAA2 <- df_mean$Group.1[which.max(df_mean$V2)]
maxAA2 <- df_mean$Group.1[which(df_mean$V2>0.4)][1]
# maxAA2 <- "Spratelloididae" 
maxAA3 <- df_mean$Group.1[which.max(df_mean$V3)]
# maxAA3 <- "Cetorhinidae" 

a1 = AAAtot
genecolour <- rep('NA', nrow(df))
genecolour2 <- rep('NA', nrow(df))



genecolour[dataphylo$Family == 'Engraulidae'] <- 'brown3' #Gobiidae #Poeciliidae
# genecolour[dataphylo$Order == maxAA2 ] <- colAA[2]  #Gobiidae #Poeciliidae
# genecolour[dataphylo$Family == maxAA2 ] <- colAA[2]  #Gobiidae #Poeciliidae
# genecolour[dataphylo$Family == 'Engraulidae'] <- 'brown3' #Gobiidae #Poeciliidae
genecolour[dataphylo$Family == 'Sparidae'] <- '#00b159' # Mugilidae # *** Labridae *** # Sparidae #Trichiuridae # Gadidae #Clupeidae # Salmonidae #Scombridae
# genecolour[dataphylo$Family == 'Psettodidae'] <- '#00b159' # Mugilidae # *** Labridae *** # Sparidae #Trichiuridae # Gadidae #Clupeidae # Salmonidae #Scombridae
# genecolour[dataphylo$Order == maxAA1] <- colAA[1] # Mugilidae # *** Labridae *** # Sparidae #Trichiuridae # Gadidae #Clupeidae # Salmonidae #Scombridae
# genecolour[dataphylo$Family == maxAA1] <- colAA[1] # Mugilidae # *** Labridae *** # Sparidae #Trichiuridae # Gadidae #Clupeidae # Salmonidae #Scombridae
genecolour[dataphylo$Order == 'Siluriformes'] <- 'royalblue' #Squlidae #Ariidae
# genecolour[dataphylo$Order == maxAA3] <-colAA[3] #Squlidae #Ariidae
# genecolour[dataphylo$Family == maxAA3] <-colAA[3] #Squlidae #Ariidae
# genecolour2[dataphylo$Family == 'Blenniidae'] <- 'black'
# genecolour2[dataphylo$Family == 'Labridae'] <- 'black'
# genecolour2[dataphylo$Class == 'Elasmobranchii'] <- 'black'

genecolour_names <-c()
genecolour_names[2] <- "Sparidae"
genecolour_names[1] <- "Siluriformes"
genecolour_names[3] <- "Engraulidae"

colnames(df)[1:4] <-  c("x", "y", "z", "d")
df <- as.data.frame(df)
df$group <- genecolour
# here : the ith label of the ternary plots is associated to the element i of the constants_x/y/z vectors
constants_x <- c(1, 0.2, 0.05)
constants_y <- c(0.1, 0.7, 0.025)
constants_z <- c(0.25, 0.2, 0.25)
dfLabs <- data.frame(x = constants_x, y = constants_y, z = constants_z)

options(tern.arrow = arrow(type = "open", length = unit(.75, "cm")))

plot <- ggtern(df,aes(x=x,y=y,z=z)) + 
  geom_point(colour = genecolour, size=2.5, alpha=0.5)+
  # theme_classic()+
  labs(
    x = labelAA[1],       # Custom name for x-axis
    y = labelAA[2],       # Custom name for y-axis
    z = labelAA[3],    # Custom name for z-axis
    colour= "Archetype"
  )+
  theme(
    tern.axis.title.T = element_text(size = 12, angle = 0, hjust = 0.5), # Top axis
    tern.axis.title.L = element_text(size = 12, angle = -58, hjust = 0.5), # Left axis
    tern.axis.title.R = element_text(size = 12, angle = -60, hjust = 1.25, vjust = -0.5), # Right axis
    plot.title = element_text(size=15, face="bold")
  )+
  ggtitle("A")+
  theme_arrowlarge()+
  theme_custom(
    base_size = 12,
    base_family = "",
    tern.plot.background = NULL,
    tern.panel.background = NULL,
    col.T = colAA[2],
    col.L = colAA[1],
    col.R = colAA[3],
    col.grid.minor = "white"
  )+
  # theme_rgbw(base_size = 12, base_family = "")+
  theme(tern.axis.arrow = element_line(size = 3),
        tern.axis.text = element_text(size=15)) +
  scale_colour_manual(values = c( maxAA1 = colAA[ieq], 
                                  maxAA2 = colAA[iopp], maxAA3= colAA[ieq]))+
  geom_label(data = dfLabs, label = genecolour_names, alpha=0.4)


plot
pdf(paste0(path_plots, "/ternary_plot_tot_spe_archetypes_time.pdf"))
print(plot)
dev.off()


blank<-grid::grid.rect(gp=gpar(col="white"))
# arrangeGrob(plot, plot1, grid)
combined <- ggtern::grid.arrange(plot, plot1)
ggsave(combined, file=paste0(path_plots, "/ternary_combined_time.pdf"), width=7 , height=9)

