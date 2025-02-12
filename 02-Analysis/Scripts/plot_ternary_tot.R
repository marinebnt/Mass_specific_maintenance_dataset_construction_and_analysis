# beneat marine 
# 14/10/24 
# Ternary plots of the archetypal analysis


OUTPUT = "Outputs_1000"
OUTPUT = "Outputs_try_again"
OUTPUT="Outputs_try_again"

colAA = c("royalblue", "darkgreen", "tomato")
labelAA = c("Opportunistic", "Equilibrium", "Periodic")

load(paste0(getwd(), "/02-Analysis/", OUTPUT,"/IMAGE_AA_FOR_ANALYSIS.RData"))
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))
path_phylosem_out <- paste0(getwd(), "/02-Analysis/", OUTPUT)
path_plots <- paste0(getwd(), "/02-Analysis/", OUTPUT,"/plots")

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
    x = "Opportunistic",       # Custom name for x-axis
    y = "Equilibrium",       # Custom name for y-axis
    z = "Periodic",        # Custom name for z-axis
  )+
  theme(
    tern.axis.title.T = element_text(size = 12, angle = 0, hjust = 0.5), # Top axis
    tern.axis.title.L = element_text(size = 12, angle = -60, hjust = 0.5), # Left axis
    tern.axis.title.R = element_text(size = 12, angle = -60, hjust = 1.25, vjust = -0.5) # Right axis
  )
plot1

ggsave(plot1, file = paste0(path_plots, "/ternary_plot_tot_time.pdf"), width=8.02, height=5.1 ) 

# ternary plot with example species
cmstd <- dataphylo$c_m
alpha <- AAAtot$alphas
df <- cbind(alpha, cmstd)

a1 = AAAtot
genecolour <- rep('NA', nrow(df))
genecolour2 <- rep('NA', nrow(df))

genecolour[dataphylo$Family == 'Gobiidae'] <- 'brown3' 
genecolour[dataphylo$Family == 'Sparidae'] <- '#00b159'
genecolour[dataphylo$Family == 'Squalidae'] <- 'royalblue' #Squlidae
genecolour2[dataphylo$Family == 'Gobiidae'] <- 'black'
genecolour2[dataphylo$Family == 'Sparidae'] <- 'black'
genecolour2[dataphylo$Family == 'Squalidae'] <- 'black'

genecolour_names <- rep(' ', length(genecolour))
genecolour_names[which(genecolour == 'brown3')[21]] <- c("Gobiidae")
genecolour_names[which(genecolour == '#00b159')[2]] <- c("Sparidae")
genecolour_names[which(genecolour == 'royalblue')[30]] <- c("Squalidae")
colnames(df) <-  c("x", "y", "z", "d")
df <- as.data.frame(df)
df$group <- genecolour

plot <- ggtern(df,aes(x=x,y=y,z=z)) + 
  geom_point(colour = genecolour, size=2.5)+
  theme_classic()+
  labs(
    x = "Opportunistic",       # Custom name for x-axis
    y = "Equilibrium",       # Custom name for y-axis
    z = "Periodic",        # Custom name for z-axis
    colour= "Archetype"
  )+
  theme(
    tern.axis.title.T = element_text(size = 12, angle = 0, hjust = 0.5), # Top axis
    tern.axis.title.L = element_text(size = 12, angle = -58, hjust = 0.5), # Left axis
    tern.axis.title.R = element_text(size = 12, angle = -60, hjust = 1.25, vjust = -0.5) # Right axis
  )+
  scale_colour_manual(values = c( "Elasmobranchii" = "royalblue", 
                                  "Gobiidae" = "tomato", "Trichiuridae"="#00b159"))+
  geom_text(aes(label = genecolour_names), 
            hjust = -0.2, 
            vjust = -1.5, 
            size = 4.5, 
            fontface = "italic",  
            color = genecolour2)
plot
ggsave(plot, file = paste0(path_plots, "/ternary_plot_tot_spe_archetypes_time.pdf"), width=8.02, height=5.1) #, width = 80, height = 180, dpi = 600, units = "mm"

par(mfrow = c(2, 2), mar = rep(0.5, 4))
plot
plot1
