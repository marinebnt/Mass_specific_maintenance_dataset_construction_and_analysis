# beneat marine 
# 14/10/24 
# Ternary plots of the archetypal analysis

path_phylosem_out <- paste0(getwd(), "/02-Analysis/Outputs")
path_plots <- paste0(getwd(), "/02-Analysis/Outputs/plots")

load(paste0(path_phylosem_out, "/IMAGE_AA_FOR_ANALYSIS.RData")) # archetypal analysis and pca outputs
source(paste0(getwd(), "/02-Analysis/Scripts/00-Functions_for_analysis.R"))

vecAA_eq = c("a","b","c")

AAAtot <- AAtot


# ternary plot with the gradient of mass specific maintenance values
beta_iv1_tot = dataplot[, which(colnames(dataplot) %in%  traits)]
t(parameters(AAAtot))     ## This shows the trait values of the selected archetypes
veccol_tot_barplot <- c("royalblue", "darkgreen", "tomato")
barplot(AAAtot, beta_iv1_tot, percentiles=TRUE, col.atypes = veccol_tot_barplot)
barplot(AAAtot, beta_iv1_tot, percentiles=TRUE)  # same but as a barplot
pcplot(AAAtot, beta_iv1_tot, atypes.col = veccol_tot_barplot, atypes.lty =c(1,2,3)) # same but as a plot showing each value
plot1 <- plotggtern(dataphylo=dataphylo, AA=AAAtot, vecAA_eq=vecAA_eq, bins=6)+
  labs(
    x = "Equilibrium",       # Custom name for x-axis
    y = "Periodic",       # Custom name for y-axis
    z = "Opportunistic",        # Custom name for z-axis
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
genecolour[dataphylo$Family == 'Squalidae'] <- 'royalblue'
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
    x = "Equilibrium",       # Custom name for x-axis
    y = "Periodic",       # Custom name for y-axis
    z = "Opportunistic",        # Custom name for z-axis
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
