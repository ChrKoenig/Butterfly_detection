library("tidyverse")
library("jagsUI")
library("coda")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())

load("Data/species_final.RData")
load("Data/traits_final.RData")

# ------------------------------------------------------------------------------------- #
#### organize SSOM model estimates ####
antilogit = function(x){exp(x) / (1+exp(x))}
SSOM_estimates = bind_rows(lapply(species_final$spec_id, FUN = function(spec){
  print(spec)
  if(!any(grepl(spec, list.files("Data/models_fit/SSOM/run_nov_19/")))){return(NULL)}
  load(list.files("Data/models_fit/SSOM/run_nov_19/", pattern = as.character(spec), full.names = T))
  tmp_summary = jags_samples$summary[,c("mean", "sd")] %>% 
    as_tibble(rownames = "param") %>% 
    mutate(param = replace(param, list = grepl("alpha_coef", param), values = c("elev", "elev_sq", "day", "day_sq"))) %>% 
    mutate(param = replace(param, list = grepl("beta_coef", param), values = c("ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq"))) %>% 
    pivot_wider(names_from = param, values_from = c(mean, sd))
}))

detection_df = right_join(species_final, model_estimates) %>% 
  left_join(traits_final) %>% 
  mutate(mean_p_conv = antilogit(mean_p_mean))

for(t in colnames(traits_final)[-1]){
  plot(pull(detection_df, "mean_p_conv") ~ pull(detection_df, t))
}

model_1 = glm(mean_p_conv ~ Vol_min, family = "binomial", data = detection_df)
summary(model_1)

# TODO add random effects for families

