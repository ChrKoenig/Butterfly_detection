library("tidyverse")
library("coda")       # mcmc/coda objects
library("runjags")    # JAGS wrapper
library("bayesplot")  # plotting mcmc/coda objects

setwd("~/ownCloud/Projects/Berlin/06_Butterfly_detection/")
rm(list=ls())
load("Data/species_final.RData")

# ------------------------------------------------------- #
#                Convergence Diagnostics               ####
# ------------------------------------------------------- #
burnin = 4000
thin = 100

## MSOM1 ####
MSOM1 = readRDS("~/Data/Butterfly_models_fit/MSOM1/MSOM1_15.RDS")      # Full MCMC chain without burnin and thinning
MSOM1_mcmc = MSOM1$mcmc[seq(from = burnin, to = 14000, by = thin),,]   # remove burnin, apply thinning
MSOM1_summary = summary(as.mcmc.list(lapply(MSOM1_mcmc, as.mcmc)))     # create summary
save(MSOM1_summary, file = "Data/MSOM1_summary.RData")

# Look at random species
spec_index = sample(1:105, 1)                                                        # Re-run multiple times   
bayesplot::mcmc_trace(MSOM1_mcmc, regex_pars = paste0("\\[", spec_index, "[^0-9]"))  # Re-run multiple times
bayesplot::mcmc_areas(MSOM1_mcmc, regex_pars = paste0("\\[", spec_index, "[^0-9]"))  # Re-run multiple times

# Look at some critical parameters
bayesplot::mcmc_trace(MSOM1_mcmc, regex_pars = "alpha_null")
bayesplot::mcmc_trace(MSOM1_mcmc, regex_pars = "beta_null")
bayesplot::mcmc_trace(MSOM1$mcmc, regex_pars = "beta_coef\\[[0-9]*,6\\]")

# --> Generally very good convergence

## MSOM2 ####
MSOM2 = readRDS("~/Data/Butterfly_models_fit/MSOM2/MSOM2_15.RDS")      # Full MCMC chain without burnin and thinning
MSOM2_mcmc = MSOM2$mcmc[seq(from = burnin, to = 14000, by = thin),,]   # remove burnin, apply thinning
MSOM2_summary = summary(as.mcmc.list(lapply(MSOM2_mcmc, as.mcmc)))     # create summary
save(MSOM2_summary, file = "Data/MSOM2_summary.RData")

# Look at random species
spec_index = sample(1:105, 1)
bayesplot::mcmc_trace(MSOM2_mcmc, regex_pars = paste0("\\[", spec_index, "[^0-9]"))
bayesplot::mcmc_areas(MSOM2_mcmc, regex_pars = paste0("\\[", spec_index, "[^0-9]"))

# Look at some critical parameters
bayesplot::mcmc_trace(MSOM2_mcmc, regex_pars = "alpha_coef_traits")
bayesplot::mcmc_areas(MSOM2_mcmc, regex_pars = "alpha_coef_traits")

bayesplot::mcmc_trace(MSOM2_mcmc, regex_pars = "alpha_coef_env\\[[0-9]*,1\\]")
bayesplot::mcmc_trace(MSOM2_mcmc, regex_pars = "beta_null")

# --> Generally very good convergence
# --> High uncertainty on some parameters

## MSOM3 ####
MSOM3 = readRDS("~/Data/Butterfly_models_fit/MSOM3/MSOM3_15.RDS")      # Full MCMC chain without burnin and thinning
MSOM3_mcmc = MSOM3$mcmc[seq(from = burnin, to = 14000, by = thin),,]   # remove burnin, apply thinning
MSOM3_summary = summary(as.mcmc.list(lapply(MSOM3_mcmc, as.mcmc)))     # create summary
save(MSOM3_summary, file = "Data/MSOM3_summary.RData")

# Look at random species
spec_index = sample(1:105, 1)
bayesplot::mcmc_trace(MSOM3_mcmc, regex_pars = paste0("\\[", spec_index, "[^0-9]"))
bayesplot::mcmc_areas(MSOM3_mcmc, regex_pars = paste0("\\[", spec_index, "[^0-9]"))

# Look at some critical parameters
bayesplot::mcmc_trace(MSOM3_mcmc, regex_pars = "alpha_coef_l2") # Trait effects on detection response to env. variables
bayesplot::mcmc_trace(MSOM3_mcmc, regex_pars = "beta_null")

# --> Okay convergence, but not on all hierarchical effects

# ------------------------------------------------------- #
#                      Response Plots                  ####
# ------------------------------------------------------- #
# Define plotting function
# TODO use posterior samples instead of posterior mean to also reflect uncertainty in model predictions
plot_response = function(model_summary, spec_id, variable, process, same_plot = F){
  # Get index of spec_id
  spec_index = which(species_final$spec_id == spec_id)
  spec_name = species_final$species[species_final$spec_id == spec_id]
  
  # Check model structure (MSOM1, MSOM2, MSOM3)
  smry = model_summary$statistics
  param_names = str_extract(rownames(smry), "[^\\[]+") %>% unique()
  if("alpha_coef_l1" %in% param_names){
    model = "MSOM3"
  } else if("alpha_coef_traits" %in% param_names){
    model = "MSOM2"
  } else {
    model = "MSOM1"
  }
  
  # Extract intercept + coefficients
  if(process == "detection"){
    traits_design_matrix = readRDS("Data/traits_design_matrix.RDS")
    intercept = switch(model,
                       "MSOM1" = smry[paste0("alpha_null[", spec_index, "]"), "Mean"],
                       "MSOM2" = smry["alpha_null", "Mean"] + smry[paste0("alpha_coef_traits[", 1:13, "]"), "Mean"] %*% as.vector(traits_design_matrix[spec_index,]),  
                       "MSOM3" = smry[paste0("alpha_null_l1[", spec_index, "]"), "Mean"])
    coefs = switch(model,
                   "MSOM1" = smry[paste0("alpha_coef_env[", spec_index, ",", 1:4, "]"), "Mean"],
                   "MSOM2" = smry[paste0("alpha_coef_env[", spec_index, ",", 1:4, "]"), "Mean"],
                   "MSOM3" = smry[paste0("alpha_coef_l1[", spec_index, ",", 1:4, "]"), "Mean"])
    names(coefs) = c("elev", "elev_sq", "day", "day_sq")
    
  } else if (process == "state"){
    intercept = smry[paste0("beta_null[", spec_index, "]"), "Mean"]
    coefs = smry[paste0("beta_coef[", spec_index, ",", 1:10, "]"), "Mean"]
    names(coefs) = c("ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq")
  } else {
    stop("unknown process")
  }

  # Prepare for plotting
  if(!variable %in% names(coefs)){
    stop("unknown variable")
  } else {
    coefs = coefs[str_detect(names(coefs), variable)]
  }
  
  # Determine variable range for plotting
  load("Data/sample_sites.RData")
  if(variable != "day"){
    var_orig = dplyr::select(sample_sites@data, contains(variable))[,1]
  } else {
    var_orig = 1:365
  }  
  var_lims = (range(var_orig) - mean(var_orig)) / sd(var_orig)
  
  # Plot response curve
  curve(plogis(intercept + coefs[1]*x + coefs[2]*x*x),
        from = var_lims[1], to = var_lims[2], ylim = c(0, 1), axes = F, xlim = c(-3,3), main = spec_name,
        xlab = variable, ylab = ifelse(process == "state", "psi (probability of occurrence)", "p (probability of detection)"))
  axis(1, at = c(-3,-2,-1,-0,1,2,3), labels = round(c(mean(var_orig)-3*sd(var_orig), mean(var_orig)-2*sd(var_orig), mean(var_orig)-sd(var_orig), mean(var_orig),
                                                      mean(var_orig)+sd(var_orig), mean(var_orig)+2*sd(var_orig), mean(var_orig)+3*sd(var_orig))))
  axis(2, las = 1)
  abline(v = var_lims[1])
  abline(v = var_lims[2])
}

# Load summaries
load("Data/MSOM1_summary.RData")
load("Data/MSOM2_summary.RData")
load("Data/MSOM3_summary.RData")

# Plot selected variables
for(spec_id in species_final$spec_id){
  plot_response(MSOM1_summary, spec_id, "elev", "detection")
  Sys.sleep(0.25)
  plot_response(MSOM2_summary, spec_id, "elev", "detection")
  Sys.sleep(0.25)
  plot_response(MSOM3_summary, spec_id, "elev", "detection")
  Sys.sleep(0.5)
}

for(spec_id in species_final$spec_id){
  plot_response(MSOM1_summary, spec_id, "day", "detection")
  Sys.sleep(0.25)
  plot_response(MSOM2_summary, spec_id, "day", "detection")
  Sys.sleep(0.25)
  plot_response(MSOM3_summary, spec_id, "day", "detection")
  Sys.sleep(0.5)
}

for(spec in species_final$spec_id){
  plot_response(MSOM1_summary, spec, "bio_12", "state")
  Sys.sleep(0.25)
  plot_response(MSOM2_summary, spec, "bio_12", "state")
  Sys.sleep(0.25)
  plot_response(MSOM3_summary, spec, "bio_12", "state")
  Sys.sleep(0.5)
}

# --> Estimated responses broadly similar among models for most species
# --> MSOM1 and MSOM3 generally more alike, while MSOM2 is often somewhat different/shifted compared to the other two models 
# --> More variation in estimated responses in detection submodel than in occupancy submodel