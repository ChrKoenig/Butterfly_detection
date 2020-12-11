library("tidyverse")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())
load("Data/species_final.RData")
load("Data/models_fit/MSOM/MSOM_3_run1.RData")

# ------------------------------------------------------------------------------------- #
summary(jags_samples)

#### Convergence ####
convergence_table = bind_cols(lapply(jags_samples$Rhat, function(x){
  if(length(dim(x)) == 0){
    return(rep(x < 1.1, 105))
  } else if(length(dim(x)) == 1){
    return(x < 1.1)
  } else if(length(dim(x)) > 1){
    apply(x, 1, function(r){all(r < 1.1)})
  }
})) %>% add_column(spec_id = species_final$spec_id, .before = "alpha_null")

table(convergence_table$alpha_null) # species with rhat > 1.1 on alpha_null

#### Convergence ####
plot_response_MSOM = function(spec, variable, process){
  # Prepare model results
  load("Data/models_fit/MSOM/MSOM_1_run1.RData")
  load("Data/species_final.RData")
  if(!spec %in% species_final$spec_id){
    stop("unknown species")
  } else {
    i = which(species_final$spec_id == spec)
  }
  if(process == "detection"){
    intercept = jags_samples$sims.list$alpha_null[,i] 
    coefs = jags_samples$sims.list$alpha_coef[,i,]
    colnames(coefs) = c("elev", "elev_sq", "day", "day_sq")
  } else if (process == "state"){
    if(length(dim(jags_samples$sims.list$beta_null)) == 2){
      intercept = jags_samples$sims.list$beta_null[,i]
    } else {
      intercept = rowMeans(jags_samples$sims.list$beta_null[,i,])
    }
    coefs = jags_samples$sims.list$beta_coef[,i,]
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
  plot_response_MSOM(spec, "elev", "detection")
  Sys.sleep(0.1)
}

for(spec in species_final$spec_id){
  plot_response_MSOM(spec, "day", "detection")
  Sys.sleep(0.1)
}

for(spec in species_final$spec_id){
  plot_response_MSOM(spec, "rad", "state")
  Sys.sleep(0.1)
}

for(spec in species_final$spec_id){
  plot_response_MSOM(spec, "ddeg0", "state")
  Sys.sleep(0.1)
}

sp = sample(species_final$spec_id, 1)
plot_response_SSOM(sp, "ddeg0", "state")
plot_response_MSOM(sp,"ddeg0", "state")

