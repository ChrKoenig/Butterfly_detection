library("tidyverse")
library("rjags")
library("coda")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())

load("Data/species_final.RData")
load("Data/traits_final.RData")
load("Data/colors_extr.RData")

# ------------------------------------------------------------------------------------- #
#### organize model estimates ####
antilogit = function(x){exp(x) / (1+exp(x))}
model_estimates = bind_rows(lapply(species_final$spec_id, FUN = function(spec){
  load(list.files("Data/model_fits/", pattern = as.character(spec), full.names = T))
  tmp_summary = summary(jags_samples)$statistics
  alpha_names = c("elev", "elev_sq", "day", "day_sq")
  beta_names = c("ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq")
  prob_names = c("mean_p", "mean_psi")
  
  tmp_means = setNames(tmp_summary[,"Mean"], paste0(c(alpha_names, beta_names, prob_names), "_mean"))
  tmp_sd = setNames(tmp_summary[,"SD"], paste0(c(alpha_names, beta_names, prob_names), "_sd"))
  return(cbind(spec_id = as.character(spec), as.data.frame(t(c(tmp_means, tmp_sd)))))
}))


detection_df = right_join(species_final, model_estimates) %>% 
  left_join(traits_final) %>% 
  mutate(mean_p_conv = antilogit(mean_p_mean))

model_1 = glm(mean_p_conv ~ generations + color + saturation + brightness + body_area, family = "binomial", data = detection_df)
summary(model_1)

# TODO add random effects for families

