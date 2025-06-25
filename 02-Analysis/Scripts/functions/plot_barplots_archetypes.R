# beneat marine 
# 14/10/24 
# to represent the barplots for the MSM paper :
# 2 analysis : 
#   MSM value elasmo//teleo
#   Barplots of the archetypal analysis (AAA) trait values


#INITIALIZATION


AAAelasmo <- AAelasmo
AAAteleo <- AAteleo
AAAtot <- AAtot

# veccol veccnames
colAAteleo <- c()
labelAAteleo <- c()
iopp <- which.max(as.data.frame(AAAteleo$archetypes)$K) 
colAAteleo[iopp] = "tomato"
labelAAteleo[iopp] = "Opportunistic"
ieq <- which.min(as.data.frame(AAAteleo$archetypes)$K) 
colAAteleo[ieq] = "royalblue"
labelAAteleo[ieq] = "Equilibrium"
iper <- c(1:3)[-c(ieq, iopp)]
colAAteleo[iper] = "darkgreen"
labelAAteleo[iper] ="Periodic"

colAAtot <- c()
labelAAtot <- c()
iopp <- which.max(as.data.frame(AAAtot$archetypes)$K) 
colAAtot[iopp] = "tomato"
labelAAtot[iopp] = "Opportunistic"
ieq <- which.min(as.data.frame(AAAtot$archetypes)$K) 
colAAtot[ieq] = "royalblue"
labelAAtot[ieq] = "Equilibrium"
iper <- c(1:3)[-c(ieq, iopp)]
colAAtot[iper] = "darkgreen"
labelAAtot[iper] ="Periodic"

colAAelasmo <- c()
labelAAelasmo <- c()
iopp <- which.max(as.data.frame(AAAelasmo$archetypes)$K) 
colAAelasmo[iopp] = "tomato"
labelAAelasmo[iopp] = "Opportunistic"
ieq <- which.min(as.data.frame(AAAelasmo$archetypes)$K) 
colAAelasmo[ieq] = "royalblue"
labelAAelasmo[ieq] = "Equilibrium"
iper <- c(1:3)[-c(ieq, iopp)]
colAAelasmo[iper] = "darkgreen"
labelAAelasmo[iper] ="Periodic"


colnames(AAAelasmo$archetypes)[which(colnames(AAAelasmo$archetypes) == "fecundity")]<- "Fecundity"
colnames(AAAteleo$archetypes)[which(colnames(AAAteleo$archetypes) == "fecundity")]<- "Fecundity"
colnames(AAAtot$archetypes)[which(colnames(AAAtot$archetypes) == "fecundity")]<- "Fecundity"



# compare RMR0 ELASMO TELEO
sd_teleo <- sd(dataphylo[which(dataphylo$Class == "Teleostei"),"c_m"])
sd_elas <- sd(dataphylo[which(dataphylo$Class == "Elasmobranchii"),"c_m"])
mean_teleo <- mean(dataphylo[which(dataphylo$Class == "Teleostei"),"c_m"])
mean_elas <- mean(dataphylo[which(dataphylo$Class == "Elasmobranchii"),"c_m"])


data <- data.frame(
  Class=c("Teleostei", "Elasmobranchii"),
  RMR0=c(mean_teleo, mean_elas),
  sd=c(sd_teleo, sd_elas)
)

stat.test <- dataphylo %>%
  t_test(c_m ~ Class) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(fun = "mean_sd", x = "Class", dodge = 0.8) 
stat.test

lm(dataphylo$c_m ~ dataphylo$Class)->a
# plot(a)  # control of the diagnostic plots 

bp <- ggbarplot(
  data, x = "Class", y = "RMR0", add = "sd", 
  fill= "Class", palette = c("grey", "black"),
  position = position_dodge(0.8)
)
bp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0.01,
  bracket.nudge.y = -1.05)+
  xlab(" ")+
  theme(legend.position = 'none')




# Barplot TOT
rownames(AAAtot$archetypes)<- labelAAtot
archetype_traits_long <- as.data.frame(AAAtot$archetypes) %>%
  tibble::rownames_to_column(var = "Archetype") %>%
  pivot_longer(cols = -Archetype, names_to = "Trait", values_to = "Value")

custom_colors <- c("Equilibrium" = "royalblue", 
                   "Periodic" = "#00b159", 
                   "Opportunistic" = "tomato")


plottot <- ggplot(archetype_traits_long, aes(x = Trait, y = Value, fill = Archetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.3, size=4) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Complete dataset archetypes trait values",
       x = "Trait",
       y = "Trait value",
       fill = "Archetype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plottot, file = paste0(path_plots, "/barplot_tot.pdf"), width=5.63, height=4.03) #, width = 80, height = 180, dpi = 600, units = "mm"


# Barplot ELASMO
rownames(AAAelasmo$archetypes)<- labelAAelasmo # c("AA1", "AA2", "AA3")
archetype_traits_long <- as.data.frame(AAAelasmo$archetypes) %>%
  tibble::rownames_to_column(var = "Archetype") %>%
  pivot_longer(cols = -Archetype, names_to = "Trait", values_to = "Value")

custom_colors <- c("Equilibrium" = "royalblue", 
                   "Periodic" = "#00b159", 
                   "Opportunistic" = "tomato")


plot <- ggplot(archetype_traits_long, aes(x = Trait, y = Value, fill = Archetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.3, size=4) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Elasmobranchs species archetypes trait values",
       x = "Trait",
       y = "Trait value",
       fill = "Archetype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot, file = paste0(path_plots, "/barplot_elasmo.pdf"), width=5.63, height=4.03) #, width = 80, height = 180, dpi = 600, units = "mm"



# Barplot TELEO
rownames(AAAteleo$archetypes)<- labelAAteleo
archetype_traits_long <- as.data.frame(AAAteleo$archetypes) %>%
  tibble::rownames_to_column(var = "Archetype") %>%
  pivot_longer(cols = -Archetype, names_to = "Trait", values_to = "Value")

custom_colors <- c("Equilibrium" = "royalblue", 
                   "Periodic" = "#00b159", 
                   "Opportunistic" = "tomato")


plotteleo <- ggplot(archetype_traits_long, aes(x = Trait, y = Value, fill = Archetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.3, size=4) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Teleost species archetypes trait values",
       x = "Trait",
       y = "Trait value",
       fill = "Archetype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plotteleo, file = paste0(path_plots, "/barplot_teleo.pdf"), width=5.63, height=4.03) 


# save combined barplots 
combined <- ggarrange(plottot, plotteleo, plot, ncol = 1, nrow=3, common.legend = T)
ggsave(combined, file=paste0(path_plots, "/barplot_combined_time.pdf"), width=7 , height=12)

