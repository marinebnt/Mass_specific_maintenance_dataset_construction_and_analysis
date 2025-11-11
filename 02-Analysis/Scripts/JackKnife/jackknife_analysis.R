library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(tidyr)

########
##### AFTER LAUNCHING JACKKNIFE : COLLECT OUTPUT ########
#########


# 
# files_jack_dataset <- list.files("01-Dataset_construction/Outputs/phylosem_final/jackknife", recursive = T, pattern = "output_SEMpsemFINALtot", full.names = T)
# files_jack_coef <- list.files("01-Dataset_construction/Outputs/phylosem_final/jackknife", recursive = T, pattern = "coef", full.names = T)
# 
# 
# for (fi in 1:length(files_jack_dataset)){
#   cat(fi, "\n")
#   # sem coefficients
#   fi_coef <- files_jack_coef[fi]
# 
#   name <- stringr::str_remove(fi_coef, "/phylosem/coefnostd_std_SEM.csv")
#   name <- stringr::str_remove(name, "01-Dataset_construction/Outputs/phylosem_final/jackknife/")
# 
#   data_i <- read.csv(fi_coef, row.names = c(2))
# 
#   if( fi == 1){dataset_coef <- data.frame("path" = rownames(data_i),data_i[, c("EstimateNotStd")],
#                                           data_i[, c("EstimateStd")])
#   } else {dataset_coef <- full_join(dataset_coef, data.frame("path" = rownames(data_i), data_i[, c("EstimateNotStd")],
#                                                           data_i[, c("EstimateStd")]))}
#   colnames(dataset_coef)[c(ncol(dataset_coef)-1):ncol(dataset_coef)]<- c(paste0(name, "_EstimateNotStd"), paste0(name, "_EstimateStd"))
# 
# 
#   # dataset
#   fi_dataset <- files_jack_dataset[fi]
#   name <- stringr::str_remove(fi_dataset, "/phylosem/output_SEMpsemFINALtot.csv")
#   name <- stringr::str_remove(name, "01-Dataset_construction/Outputs/phylosem_final/jackknife/")
#   data_i <- read.csv(fi_dataset, row.names = c(2))
#   for (tra in colnames(data_i)[-c(1,ncol(data_i))]){
#     if( fi == 1){data_tempo <- data.frame("species" = rownames(data_i), data_i[, tra])
#     } else {data_tempo <- full_join(get(paste0("dataset_", tra)), data.frame("species" = rownames(data_i), data_i[, tra]))}
#     colnames(data_tempo)[ncol(data_tempo)]<- c(paste0(name, "_", tra))
#     assign(paste0("dataset_", tra), data_tempo)
#   }
# }

pathoutput = "02-Analysis/Outputs_jackknife/Outputs_jackknife_only_c_m/"
# dir.create(pathoutput)
# output <- list(dataset_coef, dataset_c_m, dataset_fecundity, dataset_habitatbenthopelagic, dataset_habitatdemersal, dataset_habitatpelagic,
#                dataset_K, dataset_M, dataset_Lm, dataset_Loo, dataset_Woo, dataset_tmax, dataset_tm, dataset_TLDiet, dataset_Temperature,
#                dataset_Min_caudalpeduncle_depth, dataset_Max_body_width, dataset_Max_body_depth, dataset_Lower_jaw_length)
# save(output, file=paste0(pathoutput, "/data_jackknife.RData"))

load(file=paste0(pathoutput, "/data_jackknife.RData"))
dataset_coef <- output[[1]]
dataset_c_m <- output[[2]]
dataset_fecundity<- output[[3]]
dataset_habitatbenthopelagic <- output[[4]]
dataset_habitatdemersal <- output[[5]]
dataset_habitatpelagic <- output[[6]]
dataset_K <- output[[7]]
dataset_M <- output[[8]]
dataset_Lm <- output[[9]]
dataset_Loo<- output[[10]]
dataset_Woo<- output[[11]]
dataset_tmax <- output[[12]]
dataset_tm <- output[[13]]
dataset_TLDiet <- output[[14]]
dataset_Temperature <- output[[15]]
dataset_Min_caudalpeduncle_depth <- output[[16]]
dataset_Max_body_width<-output[[17]]
dataset_Max_body_depth<-output[[18]]
dataset_Lower_jaw_length<-output[[19]]
load("01-Dataset_construction/Outputs/phylosem_output/TLstdmeca/imageworkspace_phylosem.RData")



#####################
###### COMPARE INFERRED each trait  #########
#####################

# Compare with full model (no jackknife)
trait_sd       <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=nrow(df)))
trait_bias     <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=nrow(df)))
trait_unbiased <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=nrow(df)))
trait_perc_diff<- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=nrow(df)))
trait_mean     <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=nrow(df)))
trait_relative_std_erro <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=nrow(df)))
min                     <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
max                     <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
mean                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
median                  <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
minrel                     <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
maxrel                     <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
meanrel                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
medianrel                  <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q10                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q25                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q75                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q90                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q10rel                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q25rel                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q75rel                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))
q90rel                    <- data.frame(matrix(-99, ncol=ncol(dataset_traits), nrow=1))

colnames(trait_sd)       <- colnames(dataset_traits)
colnames(trait_bias)     <- colnames(dataset_traits)
colnames(trait_unbiased) <- colnames(dataset_traits)
colnames(trait_perc_diff)<- colnames(dataset_traits)
colnames(trait_mean)     <- colnames(dataset_traits)
colnames(min)            <- colnames(dataset_traits)
colnames(max)            <- colnames(dataset_traits)
colnames(mean)           <- colnames(dataset_traits)
colnames(median)         <- colnames(dataset_traits)
colnames(minrel)            <- colnames(dataset_traits)
colnames(maxrel)            <- colnames(dataset_traits)
colnames(meanrel)           <- colnames(dataset_traits)
colnames(medianrel)         <- colnames(dataset_traits)
colnames(q10)            <- colnames(dataset_traits)
colnames(q25)            <- colnames(dataset_traits)
colnames(q75)           <- colnames(dataset_traits)
colnames(q90)         <- colnames(dataset_traits)
colnames(q10rel)            <- colnames(dataset_traits)
colnames(q25rel)            <- colnames(dataset_traits)
colnames(q75rel)           <- colnames(dataset_traits)
colnames(q90rel)         <- colnames(dataset_traits)
colnames(trait_relative_std_erro)<- colnames(dataset_traits)
rownames(trait_relative_std_erro)<- rownames(dataset_traits)
rownames(trait_sd)        <- rownames(dataset_traits)
rownames(trait_bias)      <- rownames(dataset_traits)
rownames(trait_unbiased)  <- rownames(dataset_traits)
rownames(trait_perc_diff) <- rownames(dataset_traits)
rownames(trait_mean)      <- rownames(dataset_traits)

data_origin <- dataset
full_params <- df
n = nrow(data_origin)
for (trai in colnames(df)[-c(1,ncol(df))]){
  
# for (trai in colnames(df)[5]){
  cat(trai)
  params <- get(paste0("dataset_", trai))
  
  # data frame
  data_compare_trait <- left_join(data.frame("genus"=data_origin[, which(colnames(data_origin) %in% c("Genus"))], 
                                       "species"=stringr::str_replace(data_origin[, which(colnames(data_origin) %in% c("Species"))], " ", "_")), params)
  data_compare_trait <- full_join(data_compare_trait, data.frame("species"=stringr::str_replace(data_origin$Species, " ", "_"), full_params[,trai]))
  rownames(data_compare_trait) <- data_compare_trait$species
  
  if (trai == "c_m"){
    data_compare_c_m <- data_compare_trait
  }
  
  # std error of the jackknife : 
  mean_trai         <- apply(data_compare_trait[,-c(1,2,ncol(data_compare_trait))], 1, function(x) mean(as.numeric(x)))
  bias              <- (n-1)*(as.numeric(mean_trai)-full_params[, trai])
  val_unbiased      <- full_params[, trai] - bias
  sd_sp             <- apply(data_compare_trait[,-c(1:2, ncol(data_compare_trait))], 1, function(x) (((length(x)-1)^2/(length(x)))*var(as.numeric(x)))^(1/2))
  mean_sd           <- mean(sd_sp)
  min_sd            <- min(sd_sp)
  max_sd            <- max(sd_sp)
  median_sd         <- median(sd_sp)
  percent_diff      <- 100*(full_params[, trai] - as.numeric(mean_trai))/full_params[, trai]
  relative_std_erro <- sd_sp/ifelse(full_params[, trai]==0, 0.00001, abs(full_params[, trai]))*100
  mean_relsd           <- mean(relative_std_erro)
  min_relsd            <- min(relative_std_erro)
  max_relsd            <- max(relative_std_erro)
  median_relsd         <- median(relative_std_erro)
  quant_i           <- quantile(sd_sp, probs = c(0.1,0.25,0.75,0.9))
  quant_i_rel       <- quantile(relative_std_erro, probs = c(0.1,0.25,0.75,0.9))
  
  
  meanrel[, trai]         <- mean_relsd
  maxrel[, trai]          <- max_relsd
  minrel[, trai]          <- min_relsd
  medianrel[, trai]       <- median_relsd  
  mean[, trai]         <- mean_sd
  min[, trai]          <- min_sd
  max[, trai]          <- max_sd
  median[, trai]       <- median_sd
  trait_mean[, trai]   <- mean_trai
  trait_sd[,trai]      <- sd_sp
  trait_bias[,trai]    <- bias
  trait_unbiased[,trai]<- val_unbiased
  trait_perc_diff[,trai]<- percent_diff
  trait_relative_std_erro[,trai] <- relative_std_erro
  
  q10[, trai]                    <- quant_i["10%"]
  q25[, trai]                    <- quant_i["25%"]
  q75[, trai]                    <- quant_i["75%"]
  q90[, trai]                    <- quant_i["90%"]
  q10rel[, trai]                    <- quant_i_rel["10%"]
  q25rel[, trai]                    <- quant_i_rel["25%"]
  q75rel[, trai]                    <- quant_i_rel["75%"]
  q90rel[, trai]                    <- quant_i_rel["90%"]
}


path = "02-Analysis/Scripts/JackKnife"
pathoutput = "02-Analysis/Outputs_jackknife/Outputs_jackknife_only_c_m"

write.csv(data_compare_c_m, paste0(pathoutput, "/data_compare_c_m.csv"))

jack_genus <- readRDS(paste0(path, "/jack_list.rds"))
obs <- c()
jak <- c()
sp  <- c()
for (jack_i in jack_genus){
  id_rep_jack <- grep(jack_i, colnames(data_compare_c_m))
  id_spe_jack <- grep(jack_i, data_compare_c_m$genus)
  id_rep_tot  <- grep(jack_i, dataset$Genus)
  
  obs <- c(obs, dataset$c_m[id_rep_tot])
  jak <- c(jak, data_compare_c_m[id_spe_jack, id_rep_jack])
  sp  <- c(sp, dataset$Species[id_rep_tot])
}
table_jack_obs <- data.frame(obs, jak, sp)
table_jack_obs <- na.omit(table_jack_obs)

plot_j <- ggplot(table_jack_obs, aes(y=jak, x=obs))+
  geom_point()+
  geom_text_repel(aes(label=sp),
                  size = 3.5,
                  max.overlaps = 15,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color = "grey50")+
  geom_abline(aes(intercept=0, slope=1), col="red")+
  labs(title="Jackknife vs Inferred Value (RMR0 parameter)", y="Jackknife value", x="Observed value")+
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )
plot_j

pdf(paste0(pathoutput, "/jaccknife_rmr0.pdf"))
print(plot_j)
dev.off()



which(is.na(trait_mean[, trai]))   -> mean_trai
trait_sd[,trai]      <- sd_sp
trait_bias[,trai]    <- bias
trait_unbiased[,trai]<- val_unbiased
trait_perc_diff[,trai]<- percent_diff
trait_relative_std_erro[,trai] <- relative_std_erro

write.csv(trait_sd, paste0(pathoutput, "/trait_sd.csv"))
write.csv(trait_bias, paste0(pathoutput, "/trait_bias.csv"))
write.csv(trait_unbiased, paste0(pathoutput, "/trait_unbiased.csv"))
write.csv(trait_perc_diff, paste0(pathoutput, "/trait_perc_diff.csv"))
write.csv(trait_relative_std_erro, paste0(pathoutput, "/trait_relative_std_erro.csv"))




#####################
###### COMPARE PATH COEFFICIENTS #########
#####################
path = "02-Analysis/Scripts/JackKnife"
pathoutput = "02-Analysis/Outputs_jackknife/Outputs_jackknife_only_c_m"

params <- dataset_coef

# Compare with full model
params$path <- stringr::str_replace(params$path, " -> ", "->")
params$path <- stringr::str_replace(params$path, "-> ", "->")
params$path <- stringr::str_replace(params$path, " ->", "->")

full_params <- coef
colnames(coef)[1]<- "path"
n = nrow(full_params)
full_params$path <- stringr::str_replace(full_params$path, " -> ", "->")
full_params$path <- stringr::str_replace(full_params$path, "-> ", "->")
full_params$path <- stringr::str_replace(full_params$path, " ->", "->")

# data frame
data_compare_PC <- full_join(full_params, params)

write.csv(data_compare_PC, paste0(pathoutput, "/data_compare_PC.csv"))


# std error of the jackknife : 
notstd <- data_compare_PC[, grep("EstimateNotStd", colnames(data_compare_PC))]
std <- data_compare_PC[, grep("EstimateStd", colnames(data_compare_PC))]
write.csv(notstd, paste0(pathoutput, "/data_compare_PCnotstd.csv"))
write.csv(std, paste0(pathoutput, "/data_compare_PCstd.csv"))


sd_sp_NotStd <- c()
bias_NotStd  <- c()
val_unbiased_NotStd <- c()

mean_NotStd         <- apply(notstd[,-c(1,2,ncol(notstd))], 1, function(x) mean(as.numeric(x)))
bias_NotStd         <- (n-1)*(as.numeric(mean_NotStd)-full_params[, "EstimateNotStd"])
val_unbiased_NotStd <- full_params[, "EstimateNotStd"] - bias_NotStd
sd_sp_NotStd        <- apply(notstd[,-c(1)], 1, function(x) (((length(x)-1)^2/(length(x)))*var(as.numeric(x)))^(1/2))
min_sd_NotStd       <- min(sd_sp_NotStd)
max_sd_NotStd       <- max(sd_sp_NotStd)
median_sd_NotStd    <- median(sd_sp_NotStd)
mean_sd_NotStd      <- mean(sd_sp_NotStd)
percent_diff_NotStd <- 100*(full_params[, "EstimateNotStd"] - as.numeric(mean_NotStd))/full_params[, "EstimateNotStd"]
relative_std_erro_NotStd <- sd_sp_NotStd/abs(full_params[, "EstimateNotStd"])*100
# relative_std_erro_NotStd <- sd_sp_NotStd/abs(mean(full_params[, "EstimateNotStd"]))*100
min_relsd_NotStd       <- min(relative_std_erro_NotStd)
max_relsd_NotStd       <- max(relative_std_erro_NotStd)
median_relsd_NotStd    <- median(relative_std_erro_NotStd)
mean_relsd_NotStd      <- mean(relative_std_erro_NotStd)
quant_NotStd           <- quantile(sd_sp_NotStd, probs = c(0.1,0.25,0.75,0.9))
quant_rel_NotStd       <- quantile(relative_std_erro_NotStd, probs = c(0.1,0.25,0.75,0.9))
q10_NotStd          <- quant_NotStd["10%"]
q25_NotStd          <- quant_NotStd["25%"]
q75_NotStd          <- quant_NotStd["75%"]
q90_NotStd          <- quant_NotStd["90%"]
q10_rel_NotStd      <- quant_rel_NotStd["10%"]
q25_rel_NotStd      <- quant_rel_NotStd["25%"]
q75_rel_NotStd      <- quant_rel_NotStd["75%"]
q90_rel_NotStd      <- quant_rel_NotStd["90%"]

NotStd <- data.frame(sd_sp_NotStd, "PC_original"=full_params[, "EstimateNotStd"], "relative_standard_erro"=relative_std_erro_NotStd)
rownames(NotStd) <- data_compare_PC$path


sd_sp_Std <- c()
bias_Std  <- c()
val_unbiased_Std <- c()

mean_Std         <- apply(std[,-c(1,2,ncol(std))], 1, function(x) mean(as.numeric(x)))
bias_Std         <- (n-1)*(as.numeric(mean_Std)-full_params[, "EstimateStd"])
val_unbiased_Std <- full_params[, "EstimateStd"] - bias_Std
sd_sp_Std        <- apply(std[,-c(1)], 1, function(x) (((length(x)-1)^2/(length(x)))*var(as.numeric(x)))^(1/2))
min_sd_Std       <- min(sd_sp_Std)
max_sd_Std       <- max(sd_sp_Std)
median_sd_Std    <- median(sd_sp_Std)
mean_sd_Std      <- mean(sd_sp_Std)
percent_diff_Std <- 100*(full_params[, "EstimateStd"] - as.numeric(mean_Std))/full_params[, "EstimateStd"]
relative_std_erro_Std <- sd_sp_Std/abs(full_params[, "EstimateStd"])*100
# relative_std_erro_Std <- sd_sp_Std/abs(mean(full_params[, "EstimateStd"]))*100
min_relsd_Std       <- min(relative_std_erro_Std)
max_relsd_Std       <- max(relative_std_erro_Std)
median_relsd_Std    <- median(relative_std_erro_Std)
mean_relsd_Std      <- mean(relative_std_erro_Std)
quant_Std           <- quantile(sd_sp_Std, probs = c(0.1,0.25,0.75,0.9))
quant_rel_Std       <- quantile(relative_std_erro_Std, probs = c(0.1,0.25,0.75,0.9))
q10_Std          <- quant_Std["10%"]
q25_Std          <- quant_Std["25%"]
q75_Std          <- quant_Std["75%"]
q90_Std          <- quant_Std["90%"]
q10_rel_Std      <- quant_rel_Std["10%"]
q25_rel_Std      <- quant_rel_Std["25%"]
q75_rel_Std      <- quant_rel_Std["75%"]
q90_rel_Std      <- quant_rel_Std["90%"]

Std <- data.frame(sd_sp_Std, "PC_original"=full_params[, "EstimateStd"], 
                  "relative_standard_erro"=relative_std_erro_Std)
rownames(Std) <- data_compare_PC$path


write.csv(Std, paste0(pathoutput, "/pcstd.csv"))
write.csv(NotStd, paste0(pathoutput, "/pcnostd.csv"))





mmm_data_frame <- (data.frame(
           "Min_sd"=c(as.numeric(min), min_sd_Std, min_sd_NotStd),
           "Max_sd"=c(as.numeric(max), max_sd_Std, max_sd_NotStd),
           "Mean_sd"=c(as.numeric(mean), mean_sd_Std, mean_sd_NotStd),
           "Median_sd"=c(as.numeric(median), median_sd_Std, median_sd_NotStd)))
rownames(mmm_data_frame) <- c(names(median), "Standardised path coefficients", "Not standardised path coefficients")
write.csv(mmm_data_frame, paste0(pathoutput, "/Tab_max_mean_med_sd.csv"))

mmmrel_data_frame <- (data.frame(
           "Min_relative_sd"=c(as.numeric(minrel), min_relsd_Std, min_relsd_NotStd),
           "Max_relative_sd"=c(as.numeric(maxrel), max_relsd_Std, max_relsd_NotStd),
           "Mean_relative_sd"=c(as.numeric(meanrel), mean_relsd_Std, mean_relsd_NotStd),
           "Median_relative_sd"=c(as.numeric(medianrel), median_relsd_Std, median_relsd_NotStd)))
rownames(mmmrel_data_frame) <- c(names(median), "Standardised path coefficients", "Not standardised path coefficients")
write.csv(mmmrel_data_frame, paste0(pathoutput, "/Tab_max_mean_med_RELATIVEsd.csv"))


q_data_frame <- (data.frame(
  "Q10"=c(as.numeric(q10), q10_Std, q10_NotStd),
  "Q25"=c(as.numeric(q25), q25_Std, q25_NotStd),
  "Median"=c(as.numeric(median), median_sd_Std, median_sd_NotStd),
  "Q75"=c(as.numeric(q75), q75_Std, q75_NotStd),
  "Q90"=c(as.numeric(q90), q75_Std, q75_NotStd)))
rownames(q_data_frame) <- c(names(median), "Standardised path coefficients", "Not standardised path coefficients")
write.csv(q_data_frame, paste0(pathoutput, "/Tab_sd_quantiles.csv"))

q_rel_data_frame <- (data.frame(
  "Q10"=c(as.numeric(q10rel), q10_rel_Std, q10_rel_NotStd),
  "Q25"=c(as.numeric(q25rel), q25_rel_Std, q25_rel_NotStd),
  "Median"=c(as.numeric(medianrel), median_relsd_Std, median_relsd_NotStd),
  "Q75"=c(as.numeric(q75rel), q75_rel_Std, q75_rel_NotStd),
  "Q90"=c(as.numeric(q90rel), q90_rel_Std, q90_rel_NotStd)))
rownames(q_rel_data_frame) <- c(names(median), "Standardised path coefficients", "Not standardised path coefficients")
write.csv(q_rel_data_frame, paste0(pathoutput, "/Tab_relsd_quantiles.csv"))




###
#### BOXPLOT std errors #########
###
library(tidyr)
library(tidyverse)


trait_rse<-trait_relative_std_erro
summary(trait_rse)
trait_rse<-abs(trait_rse[,-c(1, 2, 3, 11)])
summary(trait_rse)
quantile(trait_rse$c_m,probs=c(0.10,0.25,0.50,0.75,0.90))

# update trait names
colnames(trait_rse) <- c("RMR0", "Maturation age invariant", 
                         "Maturation length invariant",
                         "Maximum age", "Growth coefficient", "Infinite weight",
                         "Mortality rate", "Infinite length", "Absolute fecundity", 
                         "Min caudal peduncle depth", "Lower jaw length", 
                         "Max body depth", "Max body width", "Temperature")

# Long format
trait_rse_long <- trait_rse[,!(names(trait_rse) %in% c("habitatbenthopelagic","habitatdemersal","habitatpelagic","TLDiet"))] %>% pivot_longer(cols = everything(), names_to = "variable", values_to = "valeur")

#Processing 0
trait_rse_long <- trait_rse_long %>%
  mutate(valeur_adj = 100*(valeur + min(valeur[valeur > 0]) / 2)/100)

# Median calculation and variable reordering
median_order <- trait_rse_long %>%
  group_by(variable) %>%
  summarize(mediane = median(valeur_adj, na.rm = TRUE)) %>%
  arrange(mediane) %>%
  pull(variable)

trait_rse_long$variable <- factor(trait_rse_long$variable, levels = median_order)

#trait_rse_long <- trait_rse_long %>%
#  mutate(variable = fct_reorder(variable, valeur, .fun = median, .desc = FALSE))

# Boxplot with log10 scale and graduation control
p <- ggplot(trait_rse_long, aes(x = variable, y = valeur_adj)) +
  geom_boxplot(
    fill = "white",
    outlier.shape = 21,
    outlier.fill = "grey51",
    outlier.colour = "grey51",
    outlier.alpha = 0.5
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 5, 10, 30, 100, 1000, 10000),
    labels = scales::comma
  ) +
  labs(
    x = "",
    y = "RSE % (log10)"
  ) +
  theme_bw(base_size = 14) +  # 👈 base text size (affects everything)
  theme(
    axis.title = element_text(size = 17),        # Axis titles
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),  # X tick labels
    axis.text.y = element_text(size = 12),                      # Y tick labels
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 1)
  )

p  

pdf(paste0(pathoutput, "/boxplot_RSE.pdf"), width=10, height = 8)
print(p)
dev.off()






#   
# 
# data_origin <- dataset
# full_params <- df
# for (pc in rownames(data_compare_PC)){
#   
#   # std error of the jackknife : 
#   sd_sp <- c()
#   bias <- c()
#   val_unbiased <- c()
#   for (sp in 1:nrow(data_compare_trait)){
#     mean_trai <- mean(as.numeric(data_compare_trait[sp,-c(1,2,ncol(data_compare_trait))]), na.rm = T)
#     jack_sp   <- as.numeric(data_compare_trait[sp,-c(1:2, ncol(data_compare_trait))])
#     sd_sp     <- c(sd_sp, (n-1/n)*sum((jack_sp-mean_trai)^2)^(1/2))
#     bias      <- c(bias, (n-1)*(mean_trai-full_params[sp, trai]))
#     val_unbiased <- c(val_unbiased, full_params[sp, trai] - bias)
#   }
#   
#   trait_sd[,trai]      <- sd_sp
#   trait_bias[,trai]    <- bias
#   trait_unbiased[,trai]<- val_unbiased
#   
#     # std error of the jackknife : 
#   sd_sp <- c()
#   bias <- c()
#   val_unbiased <- c()
#   for (sp in 1:nrow(data_compare_trait)){
#     mean_trai <- mean(as.numeric(data_compare_trait[sp,-c(1,2,ncol(data_compare_trait))]), na.rm = T)
#     jack_sp   <- as.numeric(data_compare_trait[sp,-c(1:2, ncol(data_compare_trait))])
#     sd_sp     <- c(sd_sp, (n-1/n)*sum((jack_sp-mean_trai)^2)^(1/2))
#     bias      <- c(bias, (n-1)*(mean_trai-full_params[sp, trai]))
#     val_unbiased <- c(val_unbiased, full_params[sp, trai] - bias)
#   }
#   
#   trait_sd[,trai]      <- sd_sp
#   trait_bias[,trai]    <- bias
#   trait_unbiased[,trai]<- val_unbiased
#   
#   
#   
#   # plnstd <-list()
#   # k=1
#   # for (i in colnames(notstd[,c(-1, -2, -ncol(notstd))])){
#   #   p<-ggplot(notstd, aes_string(x = colnames(notstd)[ncol(notstd)], y = i))+
#   #     geom_point()+
#   #     labs(x = "Inferred dataset", y = "Original dataset",
#   #          title = paste0(i, "\nJackknife influence on PC not standardized" ))+
#   #     geom_text_repel(aes(label=data_compare_PC$Parameter))+
#   #     geom_abline(aes(intercept = 0,slope = 1), col = "red")+
#   #     theme_bw()
#   #   k=k+1
#   #   plnstd[[k]] <- p
#   # }
#   # plstd <-list()
#   # k=1
#   # for (i in colnames(std[,c(-1, -2, -ncol(std))])){
#   #   p<-ggplot(std, aes((std[, ncol(std)]), (std[, i])))+
#   #     geom_point()+
#   #     labs(x = "Inferred dataset", y = "Original dataset",
#   #          title = paste0(i, "\nJackknife influence on PC not standardized" ))+
#   #     geom_text_repel(aes(label=data_compare_PC$Parameter))+
#   #     geom_abline(aes(intercept = 0,slope = 1), col = "red")+
#   #     theme_bw()
#   #   k=k+1
#   #   plstd[[k]] <- p
#   # }
#   # pdf(paste0("PLOT_jackknife_PCnotstd",pc,".pdf"))
#   # for (k in 1:length(plnstd)) print(plnstd[[k]])
#   # dev.off()
#   # 
#   # pdf(paste0("PLOT_jackknife_PCstd",pc,".pdf"))
#   # for (k in 1:length(plstd)) print(plstd[[k]])
#   # dev.off()
# }
# 
# 
# 
# 
# 
