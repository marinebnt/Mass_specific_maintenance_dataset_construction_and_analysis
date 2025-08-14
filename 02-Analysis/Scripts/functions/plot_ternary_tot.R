# beneat marine 
# 14/10/24 
# Ternary plots of the archetypal analysis


# without log_k
# colAA = c("darkgreen", "royalblue",  "tomato")
# labelAA = c("Periodic", "Equilibrium", "Opportunistic")
# i_elasmo = 2 
# i_opp = 6

# with log_k

# veccol veccnames
colAA <- c()
labelAA <- c()
pchvec <- c()
iopp <- which.max(as.data.frame(AAtot$archetypes)$K) 
colAA[iopp] = "tomato"
labelAA[iopp] = "Fast"
pchvec[iopp] = c(23)
ieq <- which.min(as.data.frame(AAtot$archetypes)$K) 
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

AAAtot <- AAtot


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
cmstd <- dataphylo$c_m
alpha <- AAAtot$alphas
df <- cbind(alpha, cmstd)

a1 = AAAtot
genecolour <- rep('NA', nrow(df))
genecolour2 <- rep('NA', nrow(df))

genecolour[dataphylo$Family == 'Gobiidae'] <- 'brown3' #Gobiidae #Poeciliidae
genecolour[dataphylo$Family == 'Mugilidae'] <- '#00b159' # Mugilidae # *** Labridae *** # Sparidae #Trichiuridae # Gadidae #Clupeidae # Salmonidae #Scombridae
genecolour[dataphylo$Family == 'Sebastidae'] <- 'royalblue' #Squlidae #Ariidae
genecolour2[dataphylo$Family == 'Blenniidae'] <- 'black'
genecolour2[dataphylo$Family == 'Labridae'] <- 'black'
genecolour2[dataphylo$Class == 'Elasmobranchii'] <- 'black'

genecolour_names <-c()
genecolour_names[iopp] <- c("Gobiidae")
genecolour_names[iper] <- c("Mugilidae")
genecolour_names[ieq] <- c("Sebastidae")

colnames(df) <-  c("x", "y", "z", "d")
df <- as.data.frame(df)
df$group <- genecolour
constants_x <- c(1, 0.1, 0.05)
constants_y <- c(0.1, 0.7, 0.01)
constants_z <- c(0.25, 0.1, 0.25)
dfLabs <- data.frame(x = constants_x, y = constants_y, z = constants_z)

plot <- ggtern(df,aes(x=x,y=y,z=z)) + 
  geom_point(colour = genecolour, size=2.5, alpha=0.5)+
  theme_classic()+
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
  scale_colour_manual(values = c( "Sebastidae" = "royalblue", 
                                  "Gobiidae" = "tomato", "Labridae"="#00b159"))+
  geom_label(data = dfLabs, label = genecolour_names, alpha=0.4)
  # geom_text(aes(label = genecolour_names), 
  #           hjust = -0.2, 
  #           vjust = -1.5, 
  #           size = 4.5, 
  #           fontface = "italic",  
  #           color = genecolour2)
plot
pdf(paste0(path_plots, "/ternary_plot_tot_spe_archetypes_time.pdf"))
print(plot)
dev.off()

# par(mfrow = c(2, 2), mar = rep(0.5, 4))
# plot
# plot1
# dev.off()
