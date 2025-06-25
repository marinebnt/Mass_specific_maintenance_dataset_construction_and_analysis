load(file = "01-Simulations/Outputs/dataset_creation_output/dataset_for_phylosem_NOUNITCV/LMER.RData")

model <- lmm_ox_S

library(lme4)
library(broom.mixed)
library(performance)
library(dplyr)

fixed_effects <- broom.mixed::tidy(lmm_ox_S, effects = "fixed")

fixed_effects

random_effects <- broom.mixed::tidy(lmm_ox_S, effects = "ran_pars") %>%
  select(group, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)
random_effects <- random_sd[,c(1,grep("sd", colnames(random_sd)))]

r2_vals <- performance::r2(lmm_ox_S)

r2_vals

summary_table <- tibble(
  Term = c("FIXED effect", "Intercept", "Temperature", "RANDOM effect", "Intercept:Genus", "Temperature:Genus"),
  Estimate = c("", fixed_effects$estimate, "", random_effects$`sd__(Intercept)`[random_effects$group == "ge"], random_effects$sd__temp[random_effects$group == "ge"]),
  SD_Intercept = random_effects$estimate[grepl("(Intercept)", random_effects$term) & (random_effects$group == "ge")],
  SD_Slope = random_effects$estimate[grepl("temp", random_effects$term)],
  SD_Residual = random_effects$estimate[random_effects$group == "Residual"],
  R2_marginal = r2_vals$R2_marginal,
  R2_conditional = r2_vals$R2_conditional
)

summary_table <- summary_table %>%
  mutate(SlopeDirection = ifelse(Term == "temp" & Estimate > 0, "+", "-"))






fixed_effects <- fixef(lmm_ox_S)
ranefs <- ranef(lmm_ox_S)$ge  # per-group intercepts and slopes
ranefs$ge <- rownames(ranefs) # add group identifier
ranefs <- ranefs %>%
  rename(
    rand_intercept = `(Intercept)`,
    rand_slope = temp
  ) %>%
  mutate(
    group_intercept = fixed_effects["(Intercept)"] + rand_intercept,
    group_slope     = fixed_effects["temp"] + rand_slope
  )
summary_table <- ranefs %>%
  select(ge, group_intercept, group_slope) %>%
  rename(Group = ge,
         Intercept = group_intercept,
         Slope = group_slope)




library(lme4)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(performance)  # for R²

# Fit model (assuming it's already done)
# lmm_ox_S <- lmer(ox ~ temp + (temp | ge), data = data_ox_fr)

# Fixed effects
fixed_effects <- fixef(lmm_ox_S)

# Random effects (BLUPs)
ranefs <- ranef(lmm_ox_S)$ge %>%
  tibble::rownames_to_column("Group") %>%
  rename(rand_intercept = `(Intercept)`, rand_slope = temp) %>%
  mutate(
    Intercept = fixed_effects["(Intercept)"] + rand_intercept,
    Slope     = fixed_effects["temp"] + rand_slope
  )

# Standard deviations of random effects
random_sd <- broom.mixed::tidy(lmm_ox_S, effects = "ran_pars") %>%
  select(group, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)
random_sd <- random_sd[,c(1,grep("sd", colnames(random_sd)))]

# Residual SD
resid_sd <- as.data.frame(VarCorr(lmm_ox_S)) %>%
  filter(is.na(var1) & is.na(var2)) %>%
  pull(sdcor)

# R² values
r2_vals <- performance::r2(lmm_ox_S)

# Add the same SDs and R²s to each row
summary_table <- ranefs %>%
  mutate(
    SD_Intercept = as.numeric(random_sd[which(random_sd$group == "ge"), grepl("(Intercept)", colnames(random_sd))]),
    SD_Slope     = as.numeric(random_sd[which(random_sd$group == "ge"), grepl("temp", colnames(random_sd))]),
    SD_Residual  = resid_sd,
    R2_marginal  = r2_vals$R2_marginal,
    R2_conditional = r2_vals$R2_conditional
  ) %>%
  select(Group, Intercept, Slope, SD_Intercept, SD_Slope, SD_Residual, R2_marginal, R2_conditional)







library(lme4)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(performance)  # for R²


fixed_effects <- fixef(lmm_ox_S)
blups <- broom.mixed::tidy(lmm_ox_S, effects = "ran_vals") %>%
  filter(group == "ge") %>%
  select(group, level, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(Group = level, Rand_Intercept = `(Intercept)`, Rand_Slope = temp)

blups <- blups %>%
  mutate(
    Intercept = fixed_effects["(Intercept)"] + Rand_Intercept,
    Slope     = fixed_effects["temp"] + Rand_Slope
  )

random_sd <- broom.mixed::tidy(lmm_ox_S, effects = "ran_pars") %>%
  select(group, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)
random_sd <- random_sd[,c(1,grep("sd", colnames(random_sd)))]


sd_intercept <- as.numeric(random_sd[which(random_sd$group == "ge"), grepl("(Intercept)", colnames(random_sd))])
sd_slope     <- as.numeric(random_sd[which(random_sd$group == "ge"), grepl("temp", colnames(random_sd))])
sd_resid     <- as.numeric(random_sd[which(random_sd$group == "Residual"), grepl("Observation", colnames(random_sd))])

r2_vals <- performance::r2(lmm_ox_S)

summary_table <- blups %>%
  mutate(
    SD_Intercept = sd_intercept,
    SD_Slope     = sd_slope,
    SD_Residual  = sd_resid,
    R2_marginal  = r2_vals$R2_marginal,
    R2_conditional = r2_vals$R2_conditional
  ) %>%
  select(Group, Intercept, Slope, Rand_Intercept, Rand_Slope,
         SD_Intercept, SD_Slope, SD_Residual, R2_marginal, R2_conditional)









#################
# Marine's version

library(lme4)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(performance)  # for R²
library(partR2)

# table summary slope and intercept per genus
fixed_effects <- fixef(lmm_ox_S)
std_error_rand <- broom.mixed::tidy(lmm_ox_S, effects = "ran_vals") %>%
  filter(group == "ge") %>%
  select(group, level, term, std.error) %>%
  pivot_wider(names_from = term, values_from = std.error) %>%
  rename(Group = level, SD_intercept = `(Intercept)`, SD_slope = temp)

table <- std_error_rand %>%
  mutate(
    Intercept = coef_R_m$`(Intercept)`,
    Slope     = coef_R_m$temp
  )

write.csv(file="01-Simulations/Outputs/dataset_creation_output/linear_model/Estimates_by_genus.csv", table[, -which(colnames(table)==c("group"))], row.names = FALSE)

# table summary fixed and random effects
  # collecting parameters one by one
random_sd <- broom.mixed::tidy(lmm_ox_S, effects = "ran_pars") %>%
  select(group, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)
random_sd <- random_sd[,c(1,grep("sd", colnames(random_sd)))]

fixed_sd <- broom.mixed::tidy(lmm_ox_S, effects = "fixed")

rand_sd_intercept <- as.numeric(random_sd[which(random_sd$group == "ge"), grepl("(Intercept)", colnames(random_sd))])
rand_sd_slope     <- as.numeric(random_sd[which(random_sd$group == "ge"), grepl("temp", colnames(random_sd))])
rand_sd_resid     <- as.numeric(random_sd[which(random_sd$group == "Residual"), grepl("Observation", colnames(random_sd))])

fixed_sd_intercept <- fixed_sd$std.error[fixed_sd$term=="(Intercept)"]
fixed_sd_slope <- fixed_sd$std.error[fixed_sd$term=="temp"]
fixed_estimate_intercept_sd <- fixed_sd$estimate[fixed_sd$term=="(Intercept)"]
fixed_estimate_slope_sd <- fixed_sd$estimate[fixed_sd$term=="temp"]
fixed_statistic_intercept <- fixed_sd$statistic[fixed_sd$term=="(Intercept)"]
fixed_statistic_slope <- fixed_sd$statistic[fixed_sd$term=="temp"]

  # to get the partial R² for each element of the model
# r2_all <- partR2(
#   lmm_ox_S,
#   partvars = "temp",
#   data = data_ox_fr,
#   nboot = 500,
#   R2_type = "conditional"    # So it includes random effects
# )


r2_vals <- performance::r2(lmm_ox_S)


table_summary <- data.frame(c("Fixed effect", "Intercept", "Slope", "Random effect", "Intercept", "Slope", "Residuals"), 
           c("", fixed_estimate_intercept_sd, fixed_estimate_slope_sd, "", "", "", ""),
           c("", fixed_sd_intercept, fixed_sd_slope, "", rand_sd_intercept, rand_sd_slope, rand_sd_resid),
           c("", fixed_statistic_intercept, fixed_statistic_slope, "", "", "", ""))
colnames(table_summary) <- c("", "Estimate", "SD", "p-value")


write.csv(file="01-Simulations/Outputs/dataset_creation_output/linear_model/Estimates_fixed_random.csv", table_summary, row.names = FALSE)
