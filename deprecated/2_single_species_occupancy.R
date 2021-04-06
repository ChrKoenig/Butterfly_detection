library("tidyverse")
library("jagsUI")
library("foreach")
library("doParallel")

setwd("~/Butterfly_project")

rm(list=ls())
load("Data/observations_completed.RData")
load("Data/sample_sites.RData")
load("Data/species_final.RData")

# ------------------------------------------------------------------------------------- #
#                        SINGLE SPECIES OCCUPANCY MODEL
# --------------------------------------------------------------------------------------#
#### Define JAGS model ####
cat(file="Butterfly_detection/jags_models/SSOM.txt", "model{
  #          Priors            #
  # -------------------------- #
  ## 1. Detection model ##
  alpha_null ~ dnorm(0, 0.0001) T(-12,12)
  for(p in 1:n_pred_det) {
    alpha_coef[p] ~ dnorm(0, 0.001) T(-12,12)
  }

  ## 2. State model ##
  beta_null ~ dnorm(0, 0.0001) T(-12,12)
  for(p in 1:n_pred_occ){
    beta_coef[p] ~ dnorm(0, 0.0001) T(-12,12)
  }        
    
  # Likelihood
  for(i in 1:n_sites){ 
    z[i] ~ dbern(psi[i]) # True Occupancy
    logit(psi[i]) <- beta_null + inprod(beta_coef, x_state[i,]) # Probability of true occupancy
    for(j in 1:n_visits){
      y[i,j] ~ dbin(mu[i,j], 2) # Probability of observation in two visits
      mu[i,j] = z[i] * p[i,j]
      logit(p[i,j]) <- alpha_null + 
                       alpha_coef[1] * elev[i]  + 
                       alpha_coef[2] * elev_sq[i] +
                       alpha_coef[3] * day[i,j] + 
                       alpha_coef[4] * day_sq[i,j]
        
      # posterior predictive distribution
      y_pr[i,j] ~ dbin(mu[i,j], 2)
      
      # Chi-square test statistic and posterior predictive check
      chi2[i,j] <- pow((y[i,j] - 2*mu[i,j]),2) / (sqrt(2*mu[i,j]) + 0.0001) # observed
      chi2_pr[i,j] <- pow((y_pr[i,j] - 2*mu[i,j]),2) / (sqrt(2*mu[i,j]) + 0.0001) # expected
    }
  }
  # Derived measures
  fit = sum(chi2)
  fit_pr = sum(chi2_pr)
}")

# ------------------------------------------------------------------------------------- #
#### Fit model, loop over all species ####
cl = makeCluster(30)
registerDoParallel(cl)
specs = setdiff(species_final$spec_id, parse_number(list.files("Data/models_fit/SSOM/run_nov_24/")))

foreach(spec = specs, .packages = c("jagsUI", "tidyverse")) %dopar% {
  #### Prepare Data ####
  # Observations
  y = observations_completed %>% 
    filter(spec_id == as.character(spec)) %>% 
    dplyr::select(preabs1:preabs7, site_id, year) %>% 
    drop_na(preabs1:preabs7) %>% 
    droplevels() 
  
  # Predictors state model
  x_state = sample_sites@data[paste(y$site_id),] %>% 
    dplyr::select(ddeg0, bio_12, rad, asp, slp) %>%  
    mutate_all(scale) %>% # center and rescale
    mutate_all(list(sq = ~.*.)) %>%  # add quadratic terms 
    as.matrix()
  
  # Predictors detection model
  x_det = observations_completed %>% 
    filter(spec_id == spec) %>% 
    dplyr::select(observer, elevation, preabs1:preabs7, day1:day7) %>% 
    drop_na(preabs1:preabs7)
  
  day = matrix(apply(select(x_det, day1:day7), 2, as.numeric), ncol = 7)
  day[is.na(day)] = mean(day, na.rm = T)  
  day = (day - mean(day)) / sd(day)
  day_sq = day*day
  elev = as.vector(scale(as.numeric(x_det$elevation)))
  elev_sq = elev*elev
  
  # Final clean up
  y = y %>% select(-site_id, -year) %>% as.matrix()
  
  # Bundle data
  jags_data = list(y = y, n_sites = nrow(y), n_visits = ncol(y), n_pred_det = 4, n_pred_occ = ncol(x_state), 
                   x_state = x_state, day = day, day_sq = day_sq, elev = elev, elev_sq = elev_sq)
  
  # Initial values
  jags_inits = function(){
    list(z = rep(1, nrow(y)))
  }
  
  # Parameters to be monitored
  jags_params = c("alpha_null", "alpha_coef", "beta_null", "beta_coef", "fit", "fit_pr")
  
  # Fit model 
  jags_samples = jags(data = jags_data, 
                      inits = jags_inits, 
                      parameters.to.save = jags_params, 
                      model.file = "Butterfly_detection/jags_models/SSOM.txt",
                      n.chains = 4, n.adapt = 300, n.burnin = 4000, n.iter = 12000, n.thin = 80, 
                      parallel = F, DIC = T)
  save(jags_samples, file = paste0("Data/models_fit/SSOM", spec, "_MCMC.RData"))
  gc()
  return(spec)
}
stopCluster(cl)
