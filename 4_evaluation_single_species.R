library("tidyverse")
library("jagsUI")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())
load("Data/species_final.RData")

# ------------------------------------------------------------------------------------- #
#### Convergence ####
load("Data/models_fit/SSOM/run_nov_24/31227_MCMC.RData") 

# look at p and psi of random model
load(sample(list.files("Data/models_fit/SSOM/run_nov_24", full.names = T), 1)) 
jagsUI::traceplot(jags_samples, parameters = c("alpha_null", "beta_null"))

# create convergence summary table
convergence_table = bind_rows(lapply(species_final$spec_id, FUN = function(spec){
  if(!any(grepl(spec, list.files("Data/models_fit/SSOM/run_nov_24/")))){return(NULL)}
  load(list.files("Data/models_fit/SSOM/run_nov_24/", pattern = as.character(spec), full.names = T))
  as.data.frame(t(sapply(jags_samples$Rhat, function(x) all(x < 1.1)))) %>% 
    add_column(spec_id = spec, .before = "alpha_null")
}))

#### Response curves ####
plot_response = function(spec, variable, process){
  # Prepare model results
  load(paste0("Data/models_fit/SSOM/run_nov_24/", spec, "_MCMC.RData"))
  if(process == "detection"){
    intercept = jags_samples$sims.list$alpha_null 
    coefs = jags_samples$sims.list$alpha_coef
    colnames(coefs) = c("elev", "elev_sq", "day", "day_sq")
  } else if (process == "state"){
    intercept = jags_samples$sims.list$beta_null
    coefs = jags_samples$sims.list$beta_coef
    colnames(coefs) = c("ddeg0","bio_12","rad","asp","slp","ddeg0_sq","bio_12_sq","rad_sq","asp_sq","slp_sq")
  } else {
    stop("unknown process")
  }
  
  if(!variable %in% colnames(coefs)){
    stop("unknown variable")
  }
  
  # Rescale variable to original range
  load("Data/sample_sites.RData")
  if(variable != "day"){
    var_orig = dplyr::select(sample_sites@data, contains(variable))[,1]
  } else {
    var_orig = 1:365
  }  
  var_lims = (range(var_orig) - mean(var_orig)) / sd(var_orig)
  
  # Plot
  for(i in 1:nrow(coefs)){
    if(i == 1){
      curve(plogis(intercept[1] + coefs[1,variable]*x + coefs[1,paste0(variable, "_sq")]*x*x), 
            from = var_lims[1], to = var_lims[2], col = "#FF000005", ylim = c(0, 1), axes = F, xlim = c(-3,3), main = spec,
            xlab = variable, ylab = ifelse(process == "state", "psi (probability of occurrence)", "p (probability of detection)"))
      axis(1, at = c(-3,-2,-1,-0,1,2,3), labels = round(c(mean(var_orig)-3*sd(var_orig), mean(var_orig)-2*sd(var_orig), mean(var_orig)-sd(var_orig), mean(var_orig), 
           mean(var_orig)+sd(var_orig), mean(var_orig)+2*sd(var_orig), mean(var_orig)+3*sd(var_orig))))
      axis(2, las = 1)
      abline(v = var_lims[1])
      abline(v = var_lims[2])
    } else (
      curve(plogis(intercept[i] + coefs[i,variable]*x + coefs[i,paste0(variable, "_sq")]*x*x),  from = var_lims[1], to = var_lims[2], add = T, col = "#FF000005")
    )
  }
}
for(spec in species_final$spec_id){
  plot_response(spec, "elev", "detection")
  Sys.sleep(0.1)
}

for(spec in species_final$spec_id){
  plot_response(spec, "slp", "state")
  Sys.sleep(0.1)
}
