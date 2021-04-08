library("tidyverse")
library("corrplot")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())

### PREPARE RESULTS ###
# ----
ssom_coefs = bind_rows(lapply(species_final$spec_id, FUN = function(spec){
  if(!any(grepl(spec, list.files("Data/models_fit/SSOM/run_nov_24/")))){return(NULL)}
  load(list.files("Data/models_fit/SSOM/run_nov_24/", pattern = as.character(spec), full.names = T))
  tmp_summary = jags_samples$summary[,c("mean")] %>% 
    as_tibble(rownames = "param") %>% 
    mutate(param = replace(param, list = grepl("alpha_coef", param), values = c("elev", "elev_sq", "day", "day_sq"))) %>% 
    mutate(param = replace(param, list = grepl("beta_coef", param), values = c("ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq"))) %>% 
    pivot_wider(names_from = param, values_from = c(value)) %>% 
    add_column(spec_id = spec, .before = "alpha_null")
})) %>% select(-fit, -fit_pr, -deviance)

# ----
msom1 = readRDS("Data/models_fit/MSOM/MSOM_1_run1.RDS")
msom1_coefs = cbind(species_final$spec_id, msom1$mean$alpha_null, msom1$mean$alpha_coef, msom1$mean$beta_null, msom1$mean$beta_coef) 
colnames(msom1_coefs) = c("spec_id","alpha_null", "elev", "elev_sq", "day", "day_sq", "beta_null", "ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq")
msom1_coefs = as_tibble(msom1_coefs)

# ----
msom2 = readRDS("Data/models_fit/MSOM/MSOM_2_run1.RDS")
msom2_coefs = cbind(species_final$spec_id, msom2$mean$alpha_null, msom2$mean$alpha_coef, msom2$mean$beta_null, msom2$mean$beta_coef) 
colnames(msom2_coefs) = c("spec_id","alpha_null", "elev", "elev_sq", "day", "day_sq", "beta_null", "ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq")
msom2_coefs = as_tibble(msom2_coefs)

# ----
msom3 = readRDS("Data/models_fit/MSOM/MSOM_3_run1.RDS")
msom3_coefs = cbind(species_final$spec_id, msom3$mean$alpha_null, msom3$mean$alpha_coef, msom3$mean$beta_null[,1], msom3$mean$beta_coef)
colnames(msom3_coefs) = c("spec_id","alpha_null", "elev", "elev_sq", "day", "day_sq", "beta_null", "ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq")
msom3_coefs = as_tibble(msom3_coefs)

### PLOT ###
plot_correlations = function(var){
  tmp_df = tibble("SSOM" = pull(ssom_coefs,var), "MSOM1" = pull(msom1_coefs,var), "MSOM2" = pull(msom2_coefs,var), "MSOM3" = pull(msom3_coefs,var)) %>% 
    mutate_if(is.character, as.numeric)
  pairs(~., tmp_df)
}

plot_correlations("alpha_null")
plot_correlations("elev")
plot_correlations("elev_sq")
plot_correlations("day")
plot_correlations("day_sq")
plot_correlations("ddeg0")
plot_correlations("bio_12")
plot_correlations("rad")
plot_correlations("asp")
plot_correlations("slp")

### Trait effects ###
load("Data/traits_final.RData")

# SSOM
ssom_trait_eff = lm(ssom_coefs$alpha_null ~ ., data = traits_final[,-1])
ssom_trait_sum = data.frame(summary(ssom_trait_eff)$coefficients[-1,1:2])
colnames(ssom_trait_sum) = c("estimate", "std_error")
ssom_trait_sum$model = "SSOM"
ssom_trait_sum$var = rownames(ssom_trait_sum)

# MSOM1
msom1_trait_eff = lm(msom1_coefs$alpha_null ~ ., data = traits_final[,-1])
msom1_trait_sum = data.frame(summary(msom1_trait_eff)$coefficients[-1,1:2])
colnames(msom1_trait_sum) = c("estimate", "std_error")
msom1_trait_sum$model = "MSOM1"
msom1_trait_sum$var = rownames(ssom_trait_sum)

# MSOM2
msom2_trait_mean = c(msom2$mean$alpha_traits_num[1:4], msom2$mean$alpha_color_bottom[-1], 
                     msom2$mean$alpha_color_top[-1], msom2$mean$alpha_traits_num[5:8])
msom2_trait_sd = c(msom2$sd$alpha_traits_num[1:4], msom2$sd$alpha_color_bottom[-1], 
                     msom2$sd$alpha_color_top[-1], msom2$sd$alpha_traits_num[5:8])
msom2_trait_sum = data.frame(estimate = msom2_trait_mean, std_error = msom2_trait_sd,
                             model = "MSOM2", var = rownames(ssom_trait_sum))

# MSOM3
msom3_trait_mean = c(msom3$mean$alpha_traits_num[1:4], msom3$mean$alpha_color_bottom[-1], 
                     msom3$mean$alpha_color_top[-1], msom3$mean$alpha_traits_num[5:8])
msom3_trait_sd = c(msom3$sd$alpha_traits_num[1:4], msom3$sd$alpha_color_bottom[-1], 
                   msom3$sd$alpha_color_top[-1], msom3$sd$alpha_traits_num[5:8])
msom3_trait_sum = data.frame(estimate = msom3_trait_mean, std_error = msom3_trait_sd,
                             model = "MSOM3", var = rownames(ssom_trait_sum))

# Plot trait effects
plot_df = bind_rows(ssom_trait_sum, msom1_trait_sum, msom2_trait_sum, msom3_trait_sum) %>%
  remove_rownames() %>% 
  mutate(model = fct_relevel(model, "SSOM"))

ggplot(plot_df, aes(x = estimate, y = var, xmin = estimate-std_error, xmax = estimate+std_error, color = model)) +
  geom_pointrange(position = position_dodge(0.5)) +
  theme_light() +


