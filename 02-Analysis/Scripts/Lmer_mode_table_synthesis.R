#* marine beneat
#* for the paper : need to summarize the linear model into a table. 

load(file = "01-Dataset_construction/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV/LMER.RData")

model <- lmm_ox_S

library(lme4)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(performance)  # for R²
library(partR2)
library(car)
library(merDeriv)


lmer(data=data_ox_fr, ox ~  temp + (temp|ge)) -> model1


cat("Step 1: Estimate column and Std error of the Fixed effect")
summary(model1)
# from this I have : The estimate of the fixed effect + The estimate of the random effect (Variance column)
# it also gives us the standard error of the fixed effect. The std error of the random effect has to be calculated. 


cat("Step 2: stat, df and p-values of the fixed effect")
# we are going to compare the model with the fixed effect, and without it
#* REFIT : we are running a likelihood test to compare the models
#* The likelihood estimator is more precise for fixed effect when we use ML instead of REML (REFIT=T)
lmer(data=data_ox_fr, ox ~  1 + (temp|ge))->model4
anova(model1, model4, refit=T)


cat("Step 3: stat, df and p-value of the random effect")
# we are comparing the model with the random effects and without it 
lmer(data=data_ox_fr, ox ~  temp + (1|ge))->model2
lmer(data=data_ox_fr, ox ~  temp + (0+temp|ge))->model3
#* REFIT : we are running a likelihood test to compare the models
#* The likelihood estimator is more precise for random effect when we use REML instead of ML (REFIT=F)
anova(model1, model2, refit=F)
anova(model1, model3, refit=F)


cat("Step 4: Estimate the standard error of the random effect")
#* 1 : are the random effects symetrically distributed ? 
#* if yes we can estimate the std error, otherwise we will keep the confidence intervals (IC) as results
confint(model)->cf # Confidence interval (IC)
data.frame(VarCorr(model))[,5]->col5 #these are the estimates of the different parts of the models
col5_bis <- col5 # reorder the estimates in the same order as the confidence intervals (see their ranges)
col5_bis[2] <- col5[3]
col5_bis[3] <- col5[2]
cf[1:4,]-col5_bis # the IC - estimates of the parameters of the model = centered IC to see if the IC are symetrical

#* 2 : the std error can be estimated because the confidence intervals are overall symetrical
#* method : https://stackoverflow.com/questions/31694812/standard-error-of-variance-component-from-the-output-of-lmer 
vv <- vcov(model, full=TRUE) # the variance covariance matrix of the model corresponds to the essian matrix of the model
                             # 1/essian matrix = is equivalent to the secondary derivative of the matrix, 
                             # which gives the slope of the log-likelihood distributions of the estimate
                             # : the wider the distribution, the less precise is the estimate
sqrt(diag(vv))
colnames(vv) # associate the names of the different parts of the model to the std errors




cat("The you can fill in the table by hand, following the comments 
    of the script to know what parameter corresponds to what.")
