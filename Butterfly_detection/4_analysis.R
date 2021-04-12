library("tidyverse")
library("MCMCvis")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())

load("Data/species_final.RData")
load("Data/traits_final.RData")
load("Data/models_fit/MSOM/MSOM_3_run1.RData")

# ------------------------------------------------------------------------------------- #
#### organize MSOM model estimates ####
antilogit = function(x){exp(x) / (1+exp(x))}
model_estimates = Reduce(cbind, list(jags_samples$mean$alpha_null, 
                                     jags_samples$mean$alpha_coef,
                                     jags_samples$mean$beta_null[,1], 
                                     jags_samples$mean$beta_coef)) %>% 
  as_tibble() %>% 
  setNames(., nm = c("alpha_null", "elev", "elev_sq", "day", "day_sq", "beta_null", "ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq")) %>% 
  add_column(spec_id = species_final$spec_id, .before = "alpha_null")

detection_df = right_join(species_final, model_estimates, by = "spec_id") %>% 
  left_join(traits_final) %>% 
  mutate(alpha_null = antilogit(alpha_null)) %>% 
  select(spec_id, alpha_null, Vol_min:mean_lgt_top) 

ggplot(detection_df, aes(x = Vol_min, y = alpha_null)) + geom_point()
ggplot(detection_df, aes(x = WIn, y = alpha_null)) + geom_point()
ggplot(detection_df, aes(x = HSI, y = alpha_null)) + geom_point()
ggplot(detection_df, aes(x = FMo_Average, y = alpha_null)) + geom_point()
ggplot(detection_df, aes(x = main_color_top, y = alpha_null)) + geom_boxplot()
ggplot(detection_df, aes(x = main_color_bottom, y = alpha_null)) + geom_boxplot()
ggplot(detection_df, aes(x = mean_sat_top, y = alpha_null)) + geom_point()
ggplot(detection_df, aes(x = mean_sat_bottom, y = alpha_null)) + geom_point()
ggplot(detection_df, aes(x = mean_lgt_top, y = alpha_null)) + geom_point()
ggplot(detection_df, aes(x = mean_lgt_bottom, y = alpha_null)) + geom_point()

MCMCplot(jags_samples, params = 'alpha_traits_num', ref_ovl = TRUE)
