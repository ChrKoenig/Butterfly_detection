library("tidyverse")
library("jagsUI")
library("coda")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())
load("Data/models_fit/MSOM/MSOM_fit_test.RData")

# ------------------------------------------------------------------------------------- #
summary(jags_samples)
samples = sample(1:105, 5)
                 
#### Convergence ####
jagsUI::traceplot(jags_samples, parameters = "alpha_null_sp[samples]") # poor convergence
jagsUI::traceplot(jags_samples, parameters = "alpha_elev[samples]") # good convergence
jagsUI::traceplot(jags_samples, parameters = "alpha_elev_sq[samples]") # good convergence
jagsUI::traceplot(jags_samples, parameters = "alpha_day[samples]") # good convergence
jagsUI::traceplot(jags_samples, parameters = "alpha_day_sq[samples]") # good convergence
jagsUI::traceplot(jags_samples, parameters = "alpha_traits_num") # very poor convergence
jagsUI::traceplot(jags_samples, parameters = "alpha_color_top") # very poor convergence
jagsUI::traceplot(jags_samples, parameters = "alpha_color_bottom") # very poor convergence
jagsUI::traceplot(jags_samples, parameters = "mu_beta_null[samples]") # good convergence
jagsUI::traceplot(jags_samples, parameters = "mu_beta_coef") # good convergence

#### Parameter estimates ####
whiskerplot(jags_samples, parameters = "alpha_null_sp[samples]")
whiskerplot(jags_samples, parameters = "mu_beta_null[samples]")
whiskerplot(jags_samples, parameters = "alpha_traits_num")
whiskerplot(jags_samples, parameters = "alpha_color_top")
whiskerplot(jags_samples, parameters = "alpha_color_bottom")


