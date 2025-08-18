# beneat marine 
# 14/10/24 
# Make sure that the dataset we use makes sense
# Question: can we have in the same dataset routine, standard and NA metabolic rates from Fishbase? 



#Load packages
library(tidyverse)
library(rstatix)
library(ggpubr)
library(lme4)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(performance)  # for R²
library(partR2)
library(car)
library(merDeriv)
library(nlme)

#####
path <- paste0("01-Dataset_construction/Scripts")
pathoutput <- paste0("01-Dataset_construction/Outputs/dataset_creation_output")

####
oxdata <- read.csv(paste0(pathoutput, "/dataset_oxygen.csv"))
spetot <- read.csv(paste0(pathoutput, "/dataset_spetot.csv"))
#####################

oxdata$MetabolicLevel[which(is.na(oxdata$MetabolicLevel))] <- "else"
oxdata$MetabolicLevel[which(oxdata$MetabolicLevel=="NA")] <- "else"
oxdata$OxygenCons <- log(oxdata$OxygenCons)



# 1- Should we include 'NA' metabolisms in the dataset ? YES

boxplot(oxdata$OxygenCons~oxdata$MetabolicLevel)
oxdata %>%
  group_by(MetabolicLevel) %>%
  get_summary_stats(OxygenCons, type = "mean_sd")

res.aov <- oxdata %>% anova_test(OxygenCons ~ MetabolicLevel)
res.aov

pwc <- oxdata %>%
  pairwise_t_test(OxygenCons ~ MetabolicLevel, p.adjust.method = "bonferroni")
pwc


pairwise.wilcox.test(oxdata$OxygenCons, oxdata$MetabolicLevel)



hist((oxdata$OxygenCons[which(oxdata$MetabolicLevel == "standard")]))
hist((oxdata$OxygenCons[which(oxdata$MetabolicLevel == "active")]))
hist((oxdata$OxygenCons[which(oxdata$MetabolicLevel == "else")]))
hist((oxdata$OxygenCons[which(oxdata$MetabolicLevel == "routine")]))




# 2- Comparing what happens to the slope and intercept when we use the dataset with only standard or routine data ? 

load(file = "01-Dataset_construction/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV/LMER.RData")
model <- lmm_ox_S
path <- paste0("01-Dataset_construction/Scripts")
pathoutput <- paste0("01-Dataset_construction/Outputs/dataset_creation_output")

# visualise different metabolisms 
library(ggplot2)
ggplot(fortify.merMod(model3), aes(temp, ox, color=mod)) +
  stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=.fitted), fun.y=mean, geom="line")

###
oxdata <- read.csv(paste0(pathoutput, "/dataset_oxygen.csv"))
spetot <- read.csv(paste0(pathoutput, "/dataset_spetot.csv"))

# linear model to infer eps_m and c_m
ox <- oxdata
residlm <- c()
residlmm <- c()
meanvarlm <- c()
meanvarlmm <- c()
outlierslm <- c()
outlierslmm <- c()

#adapt dataset to the linear model linearization + create data frame
ox$weight_g_gamma    <- (ox$Weight*10^3)^(3/4)
ox$temp_inv_k        <- -1/(ox$Temperature + 273.5)
ox$ox_div_weight_log <-  log(ox$OxygenCons/ox$weight_g_gamma)

data_ox_fr           <- data.frame(ox$weight_g_gamma, ox$temp_inv_k,
                                   ox$ox_div_weight_log, ox$DemersPelag, ox$Genus, oxdata$MetabolicLevel)
colnames(data_ox_fr) <- c("weig", "temp", "ox", "hab", "ge", "mod")

# scale the numeric variables
icol     <- which(colnames(data_ox_fr)%in% c("ox"))
torm     <- which(colnames(data_ox_fr)%in% c("weig", "ge"))
char_col <- which(sapply((data_ox_fr[1,]), is.numeric)==0)
p.order  <- c(icol,(1:ncol(data_ox_fr))[-c(icol, as.numeric(char_col), torm)])
m        <- colMeans(data_ox_fr[p.order])
s        <- apply(data_ox_fr[p.order],2,sd)

data_ox_fr[, p.order] <- sapply(data_ox_fr[, p.order], scale)
data_ox_fr$mod[which(is.na(data_ox_fr$mod))] <- "else"
data_ox_fr$mod[which(data_ox_fr$mod=="NA")] <- "else"
data_ox_fr$mod <- as.factor(data_ox_fr$mod)

data_ox_fr_rout <- data_ox_fr[which(data_ox_fr$mod == "routine"),]
data_ox_fr_std <- data_ox_fr[which(data_ox_fr$mod == "standard"),]


# => compare models #
lmer(data=data_ox_fr, ox ~  temp + (temp|ge)) -> model1 # normal
lmer(data=data_ox_fr, ox ~  temp + mod + (temp|ge)) -> model2 # normal
lmer(data=data_ox_fr, ox ~  temp * mod + (temp|ge)) -> model3 # normal

(summary(model2))
anova(model3, test.statistic = "F")

lmer(data=data_ox_fr, ox ~  temp + (temp|ge)) -> model1 # normal
lmer(data=data_ox_fr, ox ~  temp + mod + (temp|ge)) -> model2 
lmer(data=data_ox_fr, ox ~  temp * mod + (temp|ge)) -> model3
lmer(data=data_ox_fr, ox ~  temp:mod + (temp|ge)) -> model4
anova(model1)
anova(model2)
anova(model3)
anova(model1,model2) # delta AIC < 2 => 
                     # the metabolic rate does not change the results significatively, 
                     # so we can have both resting and standard metabolic rate in our dataset
anova(model3,model2)
anova(model4,model1)

# other analysis 
library(emmeans)
emmeans(model2, pairwise ~ mod)
emtrends(model2, pairwise ~ mod, var = "temp")
emmeans(model3, pairwise ~ mod)
emtrends(model3, pairwise ~ mod, var = "temp")
emmeans(model4, pairwise ~ mod)
emtrends(model4, pairwise ~ mod, var = "temp")
emm <- emmeans(model3, ~ mod)
plot(emm)


hist((data_ox_fr$ox[which(data_ox_fr$mod == "standard")]))
hist((data_ox_fr$ox[which(data_ox_fr$mod == "else")]))
hist((data_ox_fr$ox[which(data_ox_fr$mod == "routine")]))





