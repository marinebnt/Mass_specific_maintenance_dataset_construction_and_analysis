##########################
## TO USE PHYLOSEM I NEED A FEW THINGS :
## - 1 - CREATE A DATASET WITH THE OXYGEN CONSUMPTION PARAMETERS
## -> 2 - DEDUCE A RELATIONSHIP BETWEEN OXYGEN DATA AND WEIGHT AND TEMPERATURE
## - 3 - CREATE A SEPERATE DATASET WITH DATA EXTRACTED FROM FISHBASE
#        AND ADD C_M AND EPS_M DATA TO THE DATASET 
## - 4 - RUN PHYLOSEM
##########################

#Load packages
library(rfishbase)
library(dplyr)
library(car)
library(stringr)
library(ape)
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(ggpubr)
library(lme4)
library(optimx)

#####
setwd("C:/Users/mbeneat/Documents/osmose/parameterizing_ev-osmose-med/repository_for_zenodo")
path <- paste0(getwd(), "/01-Simulations/Scripts")
pathoutput <- paste0(getwd(), "/01-Simulations/Outputs/dataset_creation_output")

####
oxdata <- read.csv(paste0(pathoutput, "/dataset_oxygen.csv"))
spetot <- read.csv(paste0(pathoutput, "/dataset_spetot.csv"))
#####################
######LINEAR MODEL###
#####################

# linear model to infer eps_m and c_m 

ox <- oxdata

residlm <- c()
residlmm <- c()
meanvarlm <- c()
meanvarlmm <- c()
outlierslm <- c()
outlierslmm <- c()


#************ to compare the residuals when changing the dataset
#################  ***** comparison of the model residuals when the number of measures per genus is i = 3, 5, 10, 15
# for (i in 1:4){
# ox <- listoxdataforresid[[i]]
# dim(ox)




#adapt dataset to the linear model linearization + create data frame
ox$weight_g_gamma    <- (ox$Weight*10^3)^(3/4)
ox$temp_inv_k        <- -1/(ox$Temperature + 273.5)
ox$ox_div_weight_log <-  log(ox$OxygenCons/ox$weight_g_gamma)
data_ox_fr           <- data.frame(ox$weight_g_gamma, ox$temp_inv_k,
    ox$ox_div_weight_log, ox$DemersPelag, ox$Genus)
colnames(data_ox_fr) <- c("weig", "temp", "ox", "hab", "ge")

# create a function that will unscale the data when it got scaled
# adapted to lm and lmm 
# https://stackoverflow.com/questions/23642111/how-to-unscale-the-coefficients-from-an-lmer-model-fitted-with-a-scaled-respon
rescale.coefs <- function(beta,mu,sigma, repnb, coefname_) {
   if (repnb>1){ #lmm
        beta3 <- matrix(-99, dim(beta$ge)[1], dim(beta$ge)[2])}
   else { #lm
        beta3 <- c()}
   for (i in seq(repnb)){
      if (repnb>1){
        beta2 <- beta$ge[i,]}
      else {
        beta2 <- beta}
      temp_ID <- grep("tem", coefname_)
      beta2[temp_ID] <- sigma[-temp_ID]*beta2[temp_ID]/sigma[temp_ID]
      beta2[-c(temp_ID)]  <- 
        sigma[-temp_ID]*beta2[-temp_ID]+mu[-temp_ID]*(c(1:length(beta2))[-temp_ID]==1)-
        beta2[temp_ID]*mu[temp_ID]
      if (repnb>1){
        # beta2[-c(temp_ID)] <- beta2[-c(temp_ID)] - 
        #   rep(as.numeric(beta2[temp_ID]*mu[temp_ID]), length(beta2[-temp_ID]))
        beta3[i,] <- as.numeric(beta2)}
      else {
        # beta2[-c(temp_ID)] <- beta2[-c(temp_ID)] - beta2[temp_ID]*mu[temp_ID]
        beta3 <- beta2}
    }
    beta3
}

# linear models without scaling
lm_ox    <- lm(data=data_ox_fr, ox ~ ge+temp+temp:ge)
Anova(lm_ox)
#lmm_ox   <- lmer(data=data_ox_fr, ox ~ hab*temp+(temp|ge)) #does not run because of scales
coefname <- names(coef(lm_ox))
coefnamelm <- coefname
b1       <- coef(lm_ox)

# scale the numeric variables
icol     <- which(colnames(data_ox_fr)%in% c("ox"))
torm     <- which(colnames(data_ox_fr)%in% c("weig", "ge"))
char_col <- which(sapply((data_ox_fr[1,]), is.numeric)==0)
p.order  <- c(icol,(1:ncol(data_ox_fr))[-c(icol, as.numeric(char_col), torm)])
m        <- colMeans(data_ox_fr[p.order])
s        <- apply(data_ox_fr[p.order],2,sd)
data_ox_fr[, p.order] <- sapply(data_ox_fr[, p.order], scale)
library(doBy)
c <- scaleBy(.~., data=data_ox_fr[,c("ox","temp", "ge")], scale=FALSE)
c <- purrr::list_rbind(c)
d=c
# c=data_ox_fr

# linear models with the scaled data
lm_ox_S  <- lm(data=data_ox_fr, ox ~ ge+temp+temp:ge + 0)
lmm_ox_S <- lmer(data=data_ox_fr, ox ~  temp + (temp|ge))
lattice::dotplot(ranef(lmm_ox_S))
par(mfrow=c(2,2))
lattice::qqmath(lmm_ox_S)
plot(lmm_ox_S,
       sqrt(abs(resid(.)))~fitted(.),
       type=c("p","smooth"), col.line=1)
plot(lmm_ox_S, rstudent(.) ~ hatvalues(.))
plot(lmm_ox_S, type=c("p","smooth"), col.line=1)
lmm_ox_S_2 <- lmer(data=c, ox ~ temp+(temp|ge))
lattice::dotplot(ranef(lmm_ox_S_2))
par(mfrow=c(2,2))
lattice::qqmath(lmm_ox_S_2)
plot(lmm_ox_S_2,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1)
plot(lmm_ox_S_2, rstudent(.) ~ hatvalues(.))
plot(lmm_ox_S_2, type=c("p","smooth"), col.line=1)
coefname_m <- colnames(coef(lmm_ox_S)$ge)
b1S      <- coef(lm_ox_S)
b1S_m    <- coef(lmm_ox_S)
nb_ge <- length(table(data_ox_fr$ge))
nb_ha <- 0

############################
#PLOT Prediction
############################
# training set
dat_train <- data_ox_fr %>%
  group_by(ge) %>%
  slice(head(row_number(), 3)) %>%
  ungroup()
# testing set
dat_test <- data_ox_fr %>%
  group_by(ge) %>%
  slice(-head(row_number(), 3)) %>%
  ungroup()
## Fit mixed model
fit_lmer2 <- lmer(data=dat_train, ox ~  temp + (temp|ge))
summary(fit_lmer2)
# Predict on training set
train_preds  <- merTools::predictInterval(fit_lmer2, newdata = dat_train, n.sims = 100, returnSims = TRUE, seed = 657, level = 0.9) %>%
  as.data.frame()
dat_train <- dat_train %>% bind_cols(train_preds)
dat_train$group <- "train"
# Predict on test set with 90% prediction intervals
test_preds  <- merTools::predictInterval(fit_lmer2, newdata = dat_test, n.sims = 100, returnSims = TRUE, seed = 657, level = 0.9) %>%
  as.data.frame()
dat_test <- dat_test %>% bind_cols(test_preds)
dat_test$group <- "test"
## Combine the data together
combined_dat <- bind_rows(dat_train, dat_test) %>%
  arrange(ox, temp)
## Plot the time series of predictions and observed data
combined_dat %>%
  mutate(group = factor(group, levels = c("train", "test"))) %>%
  ggplot(aes(x = temp, y = ox)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              fill = "light grey",
              alpha = 0.8) +
  geom_line(aes(y = fit),
            col = "red",
            size = 1) +
  geom_point(aes(fill = group),
             size = 3,
             shape = 21) +
  geom_line() +
  facet_wrap(~ge) +
  theme(strip.background = element_rect(fill = "black"),
        strip.text = element_text(face = "bold", color = "white"),
        legend.position = "top") +
  labs(x = "-1/Temperature",
       y = "log(Oxygen/Weight^0.75)",
       title = "Predict validation scale globally")
boxplot(combined_dat$fit~combined_dat$group)
## Plot the time series of predictions and observed data
# training set
dat_train <- c %>%
  group_by(ge) %>%
  slice(head(row_number(), 3)) %>%
  ungroup()
# testing set
dat_test <- c %>%
  group_by(ge) %>%
  slice(-head(row_number(), 3)) %>%
  ungroup()
## Fit mixed model
fit_lmer2 <- lmer(data=dat_train, ox ~  temp + (temp|ge))
summary(fit_lmer2)
# Predict on training set
train_preds  <- merTools::predictInterval(fit_lmer2, newdata = dat_train, n.sims = 100, returnSims = TRUE, seed = 657, level = 0.9) %>%
  as.data.frame()
dat_train <- dat_train %>% bind_cols(train_preds)
dat_train$group <- "train"
# Predict on test set with 90% prediction intervals
test_preds  <- merTools::predictInterval(fit_lmer2, newdata = dat_test, n.sims = 100, returnSims = TRUE, seed = 657, level = 0.9) %>%
  as.data.frame()
dat_test <- dat_test %>% bind_cols(test_preds)
dat_test$group <- "test"
## Combine the data together
combined_dat <- bind_rows(dat_train, dat_test) %>%
  arrange(ox, temp)
## Plot the time series of predictions and observed data
combined_dat %>%
  mutate(group = factor(group, levels = c("train", "test"))) %>%
  ggplot(aes(x = temp, y = ox)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "light grey",
              alpha = 0.8) +
  geom_line(aes(y = fit),
            col = "red",
            size = 1) +
  geom_point(aes(fill = group),
             size = 3,
             shape = 21) +
  geom_line() +
  facet_wrap(~ge) +
  theme(strip.background = element_rect(fill = "black"),
        strip.text = element_text(face = "bold", color = "white"),
        legend.position = "top") +
  labs(x = "-1/Temperature",
       y = "log(Oxygen/Weight^0.75)",
       title = "Predict validation scale grouped by genus")
par(mfrow=c(1,1))
boxplot(combined_dat$fit~combined_dat$group)
#############



# rescale the output coefficients from the model
coef_R <- b1S
coef_R   <- rescale.coefs(b1S, as.numeric(c(rep(m[1], nb_ge+nb_ha), rep(m[2], nb_ge+nb_ha))),
 as.numeric(c(rep(s[1], nb_ge+nb_ha), rep(s[2], nb_ge + nb_ha))), repnb=1, coefname_=coefname)
coef_R_m <- b1S_m$ge
coef_R_m <- rescale.coefs(b1S_m, as.numeric(c(m[1], m[2])),
 as.numeric(c(s[1], s[2])), repnb=nb_ge, coefname=coefname_m)
all.equal(b1,coef_R) # compare the outputs from coefficient of scaled and unscaled model. If TRUE = perfect.
#all.equal(b1S_m,coef_R_m) # same but won't work because lmm with unscaled data isn't running


# reassamble the coefficients from lmm
coef_R_m <- as.data.frame(coef_R_m)
colnames(coef_R_m) <- colnames(b1S_m$ge)
rownames(coef_R_m) <- rownames(b1S_m$ge)
coef_R_m # check the values

# plot(lmm_genOFF)
plot(resid(lmm_ox_S)~fitted((lmm_ox_S)))
abline(h=0, col="red")
qqnorm(resid(lmm_ox_S))
qqline(resid(lmm_ox_S), col="black")
plot(lmm_ox_S,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1)
plot(lmm_ox_S, rstudent(.) ~ hatvalues(.))
#plot(lmm_genOFF2)
#Anova(lmm_ox_S)

residmesurelm <- resid(lm_ox_S)
residmesurelmm <- resid(lmm_ox_S)

# residlm[i] <- mean(residmesurelm)
# residlmm[i] <- mean(residmesurelmm)
# meanvarlm[i] <- var(residmesurelm)
# meanvarlmm[i] <- var(residmesurelmm)
# outlierslm[i] <- max(boxplot.stats(residmesurelm)$out)
# outlierslmm[i] <- max(boxplot.stats(residmesurelmm)$out)
# }

#plot(resid(lmm_ox_S)~ox$Genus)

###############################
#### extract the linear model data 
##############################

# First : linear model 
coeftot  <- coef_R
coefname <- coefnamelm
boltzmann<- 8.62*10^-5 
gamma = 0.75
habitat_genus    <- table(oxdata$Genus, oxdata$DemersPelag) 
id_habitat_genus <- rownames(habitat_genus)

c_m   <- rep(-99, length(spetot$Species))
eps_m <- rep(-99, length(spetot$Species))
for (i  in  seq_along(spetot$Species)){ # in 249){ #
  id_genus    <- which(id_habitat_genus == spetot$Genus[i])
  id_habitat  <- names(which(habitat_genus[id_genus,]>0))
  
  if (spetot$DemersPelag[i] %in% id_habitat){
    genspe   <- paste0("ge", spetot$Genus[i])
    #demerspe <- paste0("hab", spetot$DemersPelag[i])
    
    intercept <- coeftot[1]
    tem       <- coeftot[which(coefname == c("temp"))]
    
    gen =0 ; dem = 0 ; temgen = 0 ; temdemer = 0
    coef <- coeftot[which(coefname == genspe)]
    if (length(which(coefname == genspe))>0)   {
      if(!is.na(coef)){gen    <-coef}
      else {gen <- 0}}
    else {
      if (genspe == "geAcipenser") {gen <- 0}
      else {gen       <- NA}
    }
    # coef <- coeftot[which(coefname == demerspe)]
    # if (length(which(coefname == demerspe))>0) {
    #   if(!is.na(coef)){dem    <-coef}
    #   else {dem <- 0}}
    # else {
    #   if (demerspe == "habbenthopelagic") {dem <- 0}
    #   else {dem       <- NA}
    # }
    coef <- coeftot[which(coefname == paste0(genspe, ":", "temp"))]
    if (length(which(coefname == paste0(genspe, ":", "temp")))>0)  {
      if(!is.na(coef)){temgen    <-coef}
      else {temgen <- 0}}
    else {
      if (genspe == "geAcipenser") {temgen <- 0}
      else {temgen       <- NA}
    }
    # coef <- coeftot[which(coefname == paste0(demerspe, ":", "temp"))]
    # if (length(which(coefname == paste0(demerspe,  ":", "temp")))>0) {
    #   if(!is.na(coef)){temdemer    <-coef}
    #   else {temdemer <- 0}}
    # else {
    #   if (demerspe == "habbenthopelagic") {temdemer <- 0}
    #   else {temdemer       <- NA}
    # }
    c_m[i]        <- intercept + gen #+ dem
    eps_m[i]      <- (temgen + tem + temdemer)*boltzmann
  }
  else {
    c_m[i]        <- NA
    eps_m[i]      <- NA
  }
  
}
spetotm <- cbind(spetot, c_m, eps_m)




# Second : the mixed linear model
coeftot  <- coef_R_m
coefname <- rownames(coef_R_m)
boltzmann<- 8.62*10^-5
habitat_genus    <- table(oxdata$Genus, oxdata$DemersPelag) 
id_habitat_genus <- rownames(habitat_genus)

c_mm   <- rep(-99, length(spetot$Species))
eps_mm <- rep(-99, length(spetot$Species))
for (i  in seq_along(spetot$Species)){ # in 249){ #
  id_genus    <- which(id_habitat_genus == spetot$Genus[i])
  id_habitat  <- names(which(habitat_genus[id_genus,]>0))
  
  if (spetot$DemersPelag[i] %in% id_habitat){
    genspe   <- spetot$Genus[i]
    demerspe <- paste0("hab", spetot$DemersPelag[i])
    
    gen =0 ; temp =0 ; dem = 0 ; temhab =0
    coef <- coefname == genspe
    ID <- which(coef)
    if (sum(coef)>0)   {
      ID  <- which(coef)
      gen <- coeftot$`(Intercept)`[ID]
      temp<- coeftot$temp[ID]
      # if(demerspe=="habbenthopelagic"){dem = 0}
      # else {dem <- coeftot[ID,demerspe]}
      # if(demerspe=="habbenthopelagic"){temhab = 0}
      # else {temhab <- coeftot[ID,paste0(demerspe, ":temp")]}
    }
    else {
      gen       <- NA
      dem       <- NA
      temp      <- NA
      temhab    <- NA}
    c_mm[i]   <- gen #+ dem
    eps_mm[i] <- (temp )*boltzmann #+ temhab
  }
  else {
    c_mm[i]     <- NA
    eps_mm[i]   <- NA
  }
}
spetotmm <- cbind(spetot, c_mm, eps_mm)
head(spetotmm)

datajoined <- left_join(oxdata, spetotmm, by="Species")
write.csv(datajoined, paste0(pathoutput, "/dataset_oxygen_completed.csv"))


# plot c_m ~ c_mm and eps_m ~eps_mm
par(mfrow=c(2,1))

c_mm_na <- na.omit(c_mm)
c_m_na <- na.omit(c_m)
#plot(c_m_na~c_mm_na)
locm <- loess(c_m_na~c_mm_na, span=1, degree=1)
lox <- seq(min(c_mm_na), max(c_mm_na), length.out=50)
loy <- predict(locm, lox)
lines(lox, loy, col="red")
abline(a=0, b=1)

c_mm_na <- na.omit(c_mm)
c_m_na <- na.omit(c_m)
plot(log(c_m_na)~log(c_mm_na))
locm <- loess(log(c_m_na)~log(c_mm_na), span=1, degree=1)
lox <- seq(min(log(c_mm_na)), max(log(c_mm_na)), length.out=50)
loy <- predict(locm, lox)
lines(lox, loy, col="red")
abline(a=0, b=1)

eps_mm_na <- na.omit(eps_mm)
eps_m_na <- na.omit(eps_m)
plot(eps_m_na~eps_mm_na)
loepsm <- loess(eps_m_na~eps_mm_na, span=1, degree=1)
lox <- seq(min(eps_mm_na), max(eps_mm_na), length.out=50)
loy <- predict(loepsm, lox)
lines(lox, loy, col="red")
abline(a=0, b=1)

# plot the linear models residuals according to the genus
genustohighlight_eps_m <- spetotm$Genus[which(abs(spetotm$eps_m)>4)]
genustohighlight_c_m <- spetotm$Genus[which(abs(spetotm$c_m)>200)]
#Anova(lmm_genOFF2)
par(mfrow=c(2,1))
dataframe_residlm<-data.frame(c(rep("lmm_genOFF", length(data_ox_fr$ge)),
                                rep("lm_genOFF", length(data_ox_fr$ge))),
                      c(resid(lmm_ox_S), resid(lm_ox_S)), 
                      c(data_ox_fr$ge, data_ox_fr$ge))
colnames(dataframe_residlm) <- c("ref", "lm", "genus")
a <- ggplot(dataframe_residlm[which(dataframe_residlm$ref=="lm_genOFF"),], aes(x=genus, y=lm))+
     geom_boxplot()+
     geom_point(aes(col=(ifelse(dataframe_residlm[which(dataframe_residlm$ref=="lm_genOFF"),]$genus %in% genustohighlight_c_m, 
      "c_m",   ifelse(dataframe_residlm[which(dataframe_residlm$ref=="lm_genOFF"),]$genus  %in% genustohighlight_eps_m, 
      "eps_m","not outlier")))))+
      labs(title="residuals of the lm model", color = "Outliers of :")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
b <- ggplot(dataframe_residlm[which(dataframe_residlm$ref=="lmm_genOFF"),], aes(x=genus, y=lm))+
     geom_boxplot()+
     geom_point(
      aes(col=(ifelse(dataframe_residlm[which(dataframe_residlm$ref=="lmm_genOFF"),]$genus %in% genustohighlight_c_m, 
      "c_m",   ifelse(dataframe_residlm[which(dataframe_residlm$ref=="lmm_genOFF"),]$genus %in% genustohighlight_eps_m, 
      "eps_m","not outlier")))))+
      labs(title="residuals of the lmm model", color = "Outliers of :")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

d<-ggarrange(a, b, ncol = 1, nrow = 2, labels = "AUTO",  common.legend = TRUE, legend="right")
if(!dir.exists(paste0(pathoutput, "/plot_02"))){dir.create(paste0(pathoutput, "/plot_02"), showWarnings = F)}
ggsave(d, filename=paste0(pathoutput, "/plot_02/", "boxplot_genus_lm.jpg"), width=c(15), height =c(30))


#* CHANGE C_M AND EPS_M UNIT
#*
#******explanation : ****************
#* at first, c_m extracted with fishbase is ln(c_m)
#* and because ConsoOx ~ c_m*weight*f(temp), and ConsoOx in mg of O2 * kg^-1 * h^-1
#* then c_m in log(mg of O2 * kg^{-1} * h^-1)
#* the goal is to have c_m in g * g^{-beta} year^-1
#* Following previous conversions : we have c_m in log(mg of O2 * kg^{-1-beta} * h^-1) ie log(mg of O2 * kg^{-beta} * h^-1)

########
#___________________________________Unit conversion of the linear model outputs
########
#semID = 2
#c_mosmo   <- read.csv(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
#eps_mosmo <- (6.78*10^(-6) - 2.146*10^(-6)*c_mosmo$c_m)*(-1) # *According to the linear model . (-1) because in Arrhenius' formula : (-eps_m/(k_b*T))
#* log(mg of O2 * g^{-beta} * h^-1) -> mg of O2 * g^{-beta} * h^-1
spetotm$c_m <- exp(spetotm$c_m)
#* mg of O2 * g^{-beta} * h^-1  -> mmol of O2 * g^{-beta} * h^-1 : divided by 32
#* (molar mass) 
spetotm$c_m <- spetotm$c_m/32
#* mmol of O2 * g^{-beta} * h^-1 -> J * g^{-beta} * h^-1 : * 434
#* (Clarke & Johnston 1999)
spetotm$c_m <- spetotm$c_m*434


#___________________________________ Parenthesis
# compare c_m and eps_m values with the ones from Gillooly (_g) et al 2001
# to compare with Gillooly : needs to be Watt/g^3/4
#* J * g^{-beta} * h^-1 -> J * g^{-beta} * s^-1  : /3600
cm_compare_with_g <- log(spetotm$c_m/3600)
# considering that the oxygen consumption is in mg of O2/h in my dataset
c_m_g       <- 14.47 # W/g^3/4
eps_m_gfish <- 5.02  # K/1000
eps_m_gmax  <- 0.74  # eV 
eps_m_gmin  <- 0.41  # eV  

"Values for c_m"
head(na.omit(cm_compare_with_g))
range(na.omit(cm_compare_with_g))
mean(na.omit(cm_compare_with_g))
median(na.omit(cm_compare_with_g))
c_m_g

"Values for eps_m"
head(na.omit(eps_m)) 
range(na.omit(eps_m)) 
mean(na.omit(eps_m)) 
median(na.omit(eps_m)) 
eps_m_gmax  <- 0.74  # eV 
eps_m_gmin  <- 0.41  # eV  
#___________________________________Restart the unit conversion
#* J * g^{-beta} * h^-1 -> kg * g^{-beta} * h^-1 : ????
#* (Holt & Jorgensen 2014) : energy density+
EDsoma <- 4.62*10^6 #in J/kg
EDgona <- 6.93*10^6
spetotm$c_m <- spetotm$c_m/mean(EDsoma, EDgona)
#* kg * g^{-beta} * h^-1 -> g * g^{-beta} * h^-1 
spetotm$c_m <- spetotm$c_m*10^3
#* g * g^{-beta} * h^-1 -> g * g^{-beta} * y^-1
spetotm$c_m <- spetotm$c_m*(24*365)

# c_meps_mFINAL_afterconversion <- data.frame(c_mosmo$label, c_mosmo$c_m, eps_mosmo)
# colnames(c_meps_mFINAL_afterconversion) <- c("species", "c_m (g * g^{-1-beta} * y^-1)", "eps_m")
# write.csv(c_meps_mFINAL_afterconversion, paste0(pathoutput, "/", "c_meps_mFINAL_afterconversion.csv"))
# 



cor(na.omit(spetotm$eps_m), log(na.omit(spetotm$c_m)))
cor(na.omit(spetotmm$eps_m), log(na.omit(spetotmm$c_m)))
cor(na.omit(spetotm$eps_m), na.omit(spetotm$c_m))
cor(na.omit(spetotmm$eps_m), na.omit(spetotmm$c_m))
lmepsmcm <- lm(eps_mm~log(c_mm), data=spetotmm)


########
#___________________________________Unit conversion of the mixed linear model outputs
########
#semID = 2
#c_mosmo   <- read.csv(paste0(pathoutput, "/", "c_mosmose", "_", modelname[semID], "psemFINALtot", ".csv"))
#eps_mosmo <- (6.78*10^(-6) - 2.146*10^(-6)*c_mosmo$c_m)*(-1) # *According to the linear model . (-1) because in Arrhenius' formula : (-eps_m/(k_b*T))
#* log(mg of O2 * g^{-beta} * h^-1) -> mg of O2 * g^{-beta} * h^-1
spetotmm$c_mm <- exp(spetotmm$c_mm)
#* mg of O2 * g^{-beta} * h^-1  -> mmol of O2 * g^{-beta} * h^-1 : divided by 32
#* (molar mass) 
spetotmm$c_mm <- spetotmm$c_mm/32
#* mmol of O2 * g^{-beta} * h^-1 -> J * g^{-beta} * h^-1 : * 434
#* (Clarke & Johnston 1999)
spetotmm$c_mm <- spetotmm$c_mm*434

#___________________________________ Parenthesis
# compare c_m and eps_m values with the ones from Gillooly (_g) et al 2001
# to compare with Gillooly : needs to be Watt/g^3/4
#* J * g^{-beta} * h^-1 -> J * g^{-beta} * s^-1  : /(60*60)
cmm_compare_with_g <- log(spetotmm$c_mm/3600 )
# considering that the oxygen consumption is in mg of O2/h in my dataset
c_m_g       <- 14.47 # W/g^3/4
eps_m_gfish <- 5.02  # K/1000
eps_m_gmax  <- 0.74  # eV
eps_m_gmin  <- 0.41  # eV

"Values for c_m"
head(na.omit(cmm_compare_with_g))
range(na.omit(cmm_compare_with_g))
mean(na.omit(cmm_compare_with_g))
median(na.omit(cmm_compare_with_g))
c_m_g

"Values for eps_m"
head(na.omit(eps_mm)) 
range(na.omit(eps_mm)) 
mean(na.omit(eps_mm)) 
median(na.omit(eps_mm)) 
eps_m_gmax  <- 0.74  # eV 
eps_m_gmin  <- 0.41  # eV  

#___________________________________Restart the unit conversion
#* J * g^{-beta} * h^-1 -> kg * g^{-beta} * h^-1 : ????
#* (Holt & Jorgensen 2014) : energy density
EDsoma <- 4.62*10^6 #in J/kg
EDgona <- 6.93*10^6
spetotmm$c_mm <- spetotmm$c_mm/mean(EDsoma, EDgona)
#* kg * g^{-beta} * h^-1 -> g * g^{-beta} * h^-1 :
spetotmm$c_mm <- spetotmm$c_mm*10^3
#* g * g^{-beta} * h^-1 -> g * g^{-beta} * y^-1
spetotmm$c_mm <- spetotmm$c_mm*(24*365)




cor(na.omit(spetotm$eps_m), log(na.omit(spetotm$c_m)))
cor(na.omit(spetotmm$eps_m), log(na.omit(spetotmm$c_m)))
cor(na.omit(spetotm$eps_m), na.omit(spetotm$c_m))
cor(na.omit(spetotmm$eps_m), na.omit(spetotmm$c_m))
lmepsmcm <- lm(eps_mm~log(c_mm), data=spetotmm)


#####################
#####################
## plot(eps_m and c_m)
#####################
#####################

eps_malaia <- c(0.595682223, 0.333531516, 0.47339626, 0.268933469, 0.173632518, 0.37530696, 0.417973605,
                0.2588154, 0.477614462, 0.307223157, 0.246413775, 0.298152628, 0.458920549, 0.308743327,
                0.248132689, 0.404345929)
c_malaia   <- c(1.49E+11, 5034832.882, 540388517.7, 256579.0172, 5816.281889, 15368755.3, 52809960.9, 
                274469.8036, 2066760649, 1961447.834, 103677.1587, 1357998.53, 286079471.6, 1474805.094,
                174241.9494, 31010952.85) 

dataframe_boxplot_alaia <- data.frame(c(rep("lm", length(na.omit(spetotm$c_m))), rep("lmm", length(na.omit(spetotmm$c_mm))), rep("alaia", length(eps_malaia))),
                                            c(log(na.omit(spetotm$c_m)), log(na.omit(spetotmm$c_mm)), log(c_malaia)),
                                            c(log(na.omit(spetotm$eps_m)), log(na.omit(spetotmm$eps_mm)), log(eps_malaia)))
colnames(dataframe_boxplot_alaia) <- c("ref", "logc_m", "logeps_m")
par(mfrow=c(1,1))
ggplot(dataframe_boxplot_alaia, aes(x=ref, y=logc_m)) + 
  geom_boxplot()+
  geom_signif(comparisons = list(c("alaia", "lm"), c("lm", "lmm"), 
                                 c("alaia", "lmm")), 
              map_signif_level=TRUE)
ggplot(dataframe_boxplot_alaia, aes(x=ref, y=logeps_m)) + 
  geom_boxplot()+
  geom_signif(comparisons = list(c("alaia", "lm"), c("lm", "lmm"), 
                                 c("alaia", "lmm")), 
              map_signif_level=TRUE)
range(na.omit(c_malaia))
range(na.omit(spetotm$c_m))
range(na.omit(spetotmm$c_mm))
range(na.omit(eps_malaia))
range(na.omit(spetotm$eps_m))
range(na.omit(spetotmm$eps_mm))

par(mfrow=c(3,1))
hist(log(na.omit(spetotm$c_m)))
hist(log(na.omit(spetotmm$c_mm)))
hist(log(c_malaia))
par(mfrow=c(3,1))
hist(log(na.omit(spetotm$eps_m)))
hist(log(na.omit(spetotmm$eps_mm)))
hist(log(eps_malaia))


# correlation estimation  : carreful : set the correlation estimation after exp(log(c_m)) in the unit conversion process 
cor(na.omit(spetotm$eps_m), log(na.omit(spetotm$c_m)))
cor(na.omit(spetotmm$eps_m), log(na.omit(spetotmm$c_m)))
cor(na.omit(spetotm$eps_m), na.omit(spetotm$c_m))
cor(na.omit(spetotmm$eps_m), na.omit(spetotmm$c_m))
lmepsmcm <- lm(eps_mm~log(c_mm), data=spetotmm)


par(mfrow=c(2,2))
a <- which(!is.na(spetotm$c_m))
b <- ggplot(oxdata, aes(x=Class, y=log(OxygenCons))) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Elasmobranchii", "Teleostei"), c("Chondrostei", "Elasmobranchii"), 
                                 c("Petromyzonti", "Elasmobranchii"), c("Petromyzonti", "Teleostei"),
                                 c("Chondrostei", "Teleostei"), c("Chondrostei", "Petromyzonti")), 
              map_signif_level=TRUE)
c <- ggplot(spetotm[a,], aes(x=Class, y=log(c_m))) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Elasmobranchii", "Teleostei"), c("Chondrostei", "Elasmobranchii"), 
                                 c("Petromyzonti", "Elasmobranchii"), c("Petromyzonti", "Teleostei"),
                                 c("Chondrostei", "Teleostei"), c("Chondrostei", "Petromyzonti")),
              map_signif_level=TRUE)
d <- ggplot(spetotm[a,], aes(x=Class, y=eps_m)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Elasmobranchii", "Teleostei"), c("Chondrostei", "Elasmobranchii"), 
                                 c("Petromyzonti", "Elasmobranchii"), c("Petromyzonti", "Teleostei"),
                                 c("Chondrostei", "Teleostei"), c("Chondrostei", "Petromyzonti")),
              map_signif_level=TRUE)
e <- ggarrange(d,b,c)
e


par(mfrow=c(2,2))
a <- which(!is.na(spetotm$c_m))
ggplot(data=spetotm[a,], aes(x=eps_m, y=log(c_m), label=Genus, col=DemersPelag))+
  geom_point()+
  theme_classic()
# +
#   geom_label_repel() +
#   theme_classic()
ggplot(data=spetotm) +
  geom_point(aes(x=eps_m, y=log(c_m)), color=ifelse(spetot$Genus == "Oreochromis", 'red', 'black'))
# ggplot(data=spetotmm) +
#   geom_point(aes(x=eps_mm, y=c_mm, col=DemersPelag), color=ifelse(spetot$Genus == "Oreochromis", 'red', 'black'))
ggplot(data=spetotm) +
  geom_point(aes(x=log(c_m), y=eps_m, col=DemersPelag))
ggplot(data=spetotmm) +
  geom_point(aes(x=log(c_mm), y=eps_mm, col=DemersPelag))


##############
# write a csv file with this output
##############


spetotmm$c_m   <- spetotmm$c_mm
spetotmm$eps_m <- spetotmm$eps_mm
spetotmm       <- spetotmm[, -which(colnames( spetotmm) %in% c("c_mm", "eps_mm"))]
write.csv(spetotmm, paste0(pathoutput, "/spetot_fishabse_c_m_eps_mLMnoHABTREF.csv"))


