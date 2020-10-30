library("tidyverse")
library("rjags")
library("coda")
library("MCMCvis")
library("mcmcplots")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())
load("Data/species_final.RData")

# ------------------------------------------------------------------------------------- #
#### Convergence ####
load(sample(list.files("Data/model_fits/", full.names = T), 1)) # sample model for random species
load("Data/model_fits/31133_MCMC.RData")

MCMCtrace(jags_samples, params = "alpha_brightness")
MCMCsummary(jags_samples, excl = c("chi2", "chi2_pr"), round = 2, Rhat = T, n.eff = T)

effectiveSize(jags_samples)

#### Model fit ####
# Bayesian p-value
load("Data/models_test/31232_MCMC.RData")
samples_fit = MCMCchains(jags_samples, "fit")
samples_fit_pr = MCMCchains(jags_samples, "fit_pr")
sum(samples_fit > samples_fit_pr) / length(samples_fit)

# AUC posterior
y_pr = MCMCchains(jags_samples, "y_pr")

#### Parameter estimates ####
MCMCplot(jags_samples, params = c("mu_p", "alpha_1"), labels = c("p", "elevation", "elevation_sq", "day", "day_sq"))
MCMCplot(jags_samples, params = c("mu_psi", "beta_1"), labels = c("psi", "ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq"))

#### Response curves ####
plot_response = function(spec_id, variable, process){
  # Prepare model results
  load(paste0("Data/model_fits/", spec_id, "_MCMC.RData"))
  post_means = data.frame(summary(jags_samples)$statistics)
  post_means$var_name = c("elevation", "elevation_sq", "day", "day_sq", "ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq","mu_p","mu_psi")
  coefs = post_means %>% 
    rownames_to_column(var = "row_name") %>% 
    mutate(process = case_when(grepl("alpha", row_name) ~ "detection",
                               grepl("beta", row_name) ~ "state",
                               TRUE ~ "probability")) %>% 
    dplyr::filter(process == !!process | process == "probability")
  
  # Get coefficients
  if(process == "state"){b0 = coefs$Mean[coefs$row_name == "mu_psi"]}
  if(process == "detection"){b0 = coefs$Mean[coefs$row_name == "mu_p"]}
  b1 = coefs$Mean[coefs$var_name == variable]
  b2 = coefs$Mean[coefs$var_name == paste0(variable, "_sq")]
  
  # Rescale variable to original range
  load("Data/sample_sites.RData")
  if(variable != "day"){
    var_orig = sample_sites@data[,variable]
  } else {
    var_orig = 1:365
  }  
  var_lims = (range(var_orig) - mean(var_orig)) / sd(var_orig)
  
  # Plot
  curve(plogis(b0 + b1*x + b2*x*x), var_lims[1], var_lims[2], col = "red", ylim = c(0, 1), axes = F, xlim = c(-2.5,2.5), main = spec,
        xlab = variable, ylab = ifelse(process == "state", "psi (probability of occurrence)", "p (probability of detection)"))
  axis(1, at = c(-2,-1,-0,1,2), labels = round(c(mean(var_orig)-2*sd(var_orig), mean(var_orig)-sd(var_orig), mean(var_orig), 
                                                 mean(var_orig)+sd(var_orig), mean(var_orig)+2*sd(var_orig))))
  axis(2, las = 1)
}

for(spec in species_final$spec_id){
  plot_response(spec, "rad", "state")  
  Sys.sleep(0.1)
}

