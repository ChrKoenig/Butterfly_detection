library("tidyverse")
library("jagsUI")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())
load("Data/species_final.RData")
MSOM1 = readRDS("Data/models_fit/MSOM/MSOM_1_testrun.RDS")
MSOM2 = readRDS("Data/models_fit/MSOM/MSOM_2_testrun.RDS")
MSOM3 = readRDS("Data/models_fit/MSOM/MSOM_3_testrun.RDS")
MSOM4 = readRDS("Data/models_fit/MSOM/MSOM_4_testrun.RDS")

# ------------------------------------------------------------------------------------- #
summary(MSOM1)
summary(MSOM2)
summary(MSOM3)
summary(MSOM4)

#### Convergence ####
convergence_table = function(model){
  bind_cols(lapply(model$Rhat, function(x){
    if(length(dim(x)) == 0){
      return(rep(x < 1.1 & x > 0.9, 105))
    } else if(length(dim(x)) == 1 & length(x) == 105){
      return(x < 1.1 & x > 0.9)
    } else if(length(dim(x)) == 1 & length(x) != 105){ # random effects 
      return(rep(all(x < 1.1 & x > 0.9), 105))
    } else if(length(dim(x)) > 1){
      apply(x, 1, function(r){all(r < 1.1 & r > 0.9)})
    }
  })) %>% add_column(spec_id = species_final$spec_id, .before = 1)
}

### Overview
view(convergence_table(MSOM1))
view(convergence_table(MSOM2))
view(convergence_table(MSOM3))
view(convergence_table(MSOM4))

#### Traceplots
# MSOM1
traceplot(MSOM1, parameters = sprintf("alpha_null[c(%s)]", paste(sample(length(MSOM1$mean$alpha_null), 5), collapse = ",")))
traceplot(MSOM1, parameters = sprintf("alpha_coef_env[c(%s)]", paste(sample(length(MSOM1$mean$alpha_coef_env), 5), collapse = ",")))
traceplot(MSOM1, parameters = sprintf("beta_coef[c(%s)]", paste(sample(length(MSOM1$mean$beta_coef), 5), collapse = ",")))
traceplot(MSOM1, parameters = sprintf("beta_null[c(%s)]", paste(sample(length(MSOM1$mean$beta_null), 5), collapse = ",")))
# Overall good convergence in MSOM1

#MSOM2
traceplot(MSOM2, parameters = sprintf("alpha_null[c(%s)]", paste(sample(length(MSOM2$mean$alpha_null), 5), collapse = ",")))
traceplot(MSOM2, parameters = sprintf("alpha_coef_traits[c(%s)]", paste(sample(length(MSOM2$mean$alpha_coef_traits), 5), collapse = ",")))
traceplot(MSOM2, parameters = sprintf("alpha_coef_env[c(%s)]", paste(sample(length(MSOM2$mean$alpha_coef_env), 5), collapse = ",")))
traceplot(MSOM2, parameters = sprintf("beta_coef[c(%s)]", paste(sample(length(MSOM2$mean$beta_coef), 5), collapse = ",")))
traceplot(MSOM2, parameters = sprintf("beta_null[c(%s)]", paste(sample(length(MSOM2$mean$beta_null), 5), collapse = ",")))

#MSOM3
traceplot(MSOM3, parameters = sprintf("alpha_null[c(%s)]", paste(sample(length(MSOM3$mean$alpha_null), 5), collapse = ",")))
traceplot(MSOM3, "alpha_coef_traits")
traceplot(MSOM3, parameters = sprintf("alpha_coef_env[c(%s)]", paste(sample(length(MSOM3$mean$alpha_coef_env), 5), collapse = ",")))
traceplot(MSOM3, parameters = sprintf("beta_coef[c(%s)]", paste(sample(length(MSOM3$mean$beta_coef), 5), collapse = ",")))
traceplot(MSOM3, parameters = sprintf("beta_null[c(%s)]", paste(sample(length(MSOM3$mean$beta_null), 5), collapse = ",")))
# Generally good convergence except for color base levels

#MSOM4
traceplot(MSOM4, parameters = sprintf("alpha_null_l1[c(%s)]", paste(sample(length(MSOM4$mean$alpha_null_l1), 5), collapse = ",")))
traceplot(MSOM4, parameters = sprintf("alpha_coef_l1[c(%s)]", paste(sample(length(MSOM4$mean$alpha_coef_l1), 5), collapse = ",")))
traceplot(MSOM4, parameters = "alpha_null_l2") 
traceplot(MSOM4, parameters = "alpha_coef_l2") 
traceplot(MSOM4, parameters = sprintf("beta_null[c(%s)]", paste(sample(length(MSOM4$mean$beta_null), 5), collapse = ",")))
traceplot(MSOM4, parameters = sprintf("beta_coef[c(%s)]", paste(sample(length(MSOM4$mean$beta_coef), 5), collapse = ",")))
# Convergence problems in level 2 parameters 

#### Response plots ####
plot_response_MSOM = function(model = c("MSOM_1", "MSOM_2", "MSOM_3", "MSOM_4"), spec, variable, process, same_plot = F){
  # Get species index
  load("Data/species_final.RData")
  traits_design_matrix = readRDS("Data/traits_design_matrix.RDS")

  if(!spec %in% species_final$spec_id){
    stop("unknown species")
  } else {
    i = which(species_final$spec_id == spec)
  }
  
  # Extract model parameters
  model_fit = readRDS(paste0("Data/models_fit/MSOM/", model, "_run3.RDS"))
  n_samples = model_fit$mcmc.info$n.samples
  if(process == "detection"){
    intercept = switch(model,
      "MSOM_1" = model_fit$sims.list$alpha_null[,i],
      "MSOM_2" = model_fit$sims.list$alpha_null + model_fit$sims.list$alpha_coef_traits %*% as.vector(traits_design_matrix[i,]),  
      "MSOM_3" = model_fit$sims.list$mu_alpha_null + model_fit$sims.list$alpha_coef_traits %*% as.vector(traits_design_matrix[i,])
    )
    coefs = model_fit$sims.list$alpha_coef_env[,i,]
    colnames(coefs) = c("elev", "elev_sq", "day", "day_sq")
  } else if (process == "state"){
    if(length(dim(model_fit$sims.list$beta_null)) == 2){
      intercept = model_fit$sims.list$beta_null[,i]
    } else {
      intercept = rowMeans(model_fit$sims.list$beta_null[,i,])
    }
    coefs = model_fit$sims.list$beta_coef[,i,]
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
plot_response_MSOM_hrc = function(spec, variable, process){
  # Get species index
  load("Data/species_final.RData")
  traits_design_matrix = readRDS("Data/traits_design_matrix.RDS")
  
  if(!spec %in% species_final$spec_id){
    stop("unknown species")
  } else {
    i = which(species_final$spec_id == spec)
  }
  
  # Extract model parameters
  model_fit = readRDS(paste0("Data/models_fit/MSOM/MSOM_4_testrun.RDS"))
  n_samples = model_fit$mcmc.info$n.samples
  if(process == "detection"){
    intercept = model_fit$sims.list$mu_alpha_null_l1
    coefs = model_fit$sims.list$alpha_coef_l1[,i,]
    colnames(coefs) = c("elev", "elev_sq", "day", "day_sq")
  } else if (process == "state"){
    if(length(dim(model_fit$sims.list$beta_null)) == 2){
      intercept = model_fit$sims.list$beta_null[,i]
    } else {
      intercept = rowMeans(model_fit$sims.list$beta_null[,i,])
    }
    coefs = model_fit$sims.list$beta_coef[,i,]
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
  plot_response_MSOM_std("MSOM_1", spec, "elev", "detection")
  Sys.sleep(0.1)
  plot_response_MSOM_std("MSOM_2", spec, "elev", "detection")
  Sys.sleep(0.1)
  plot_response_MSOM_std("MSOM_3", spec, "elev", "detection")
  Sys.sleep(0.1)
  plot_response_MSOM_hrc(spec, "elev", "detection")
  Sys.sleep(0.1)
}

for(spec in species_final$spec_id){
  plot_response_MSOM_std("MSOM_1", spec, "day", "detection")
  Sys.sleep(0.1)
  plot_response_MSOM_std("MSOM_2", spec, "day", "detection")
  Sys.sleep(0.1)
  plot_response_MSOM_std("MSOM_3", spec, "day", "detection")
  Sys.sleep(0.1)
  plot_response_MSOM_hrc(spec, "day", "detection")
  Sys.sleep(0.1)
}

for(spec in species_final$spec_id){
  plot_response_MSOM_std("MSOM_1", spec, "elev", "state")
  Sys.sleep(0.1)
  plot_response_MSOM_std("MSOM_2", spec, "elev", "state")
  Sys.sleep(0.1)
  plot_response_MSOM_std("MSOM_3", spec, "elev", "state")
  Sys.sleep(0.1)
  plot_response_MSOM_hrc(spec, "elev", "state")
  Sys.sleep(0.1)
}



