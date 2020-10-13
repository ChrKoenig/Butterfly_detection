library(tidyverse)
library(rjags)
library(foreach)
library(doParallel)

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())

load("Data/observations_completed.RData")
load("Data/sample_sites.RData")
load("Data/species_final.RData")

# ------------------------------------------------------------------------------------- #
#### Define JAGS model ####
cat(file="Butterfly_detection/jags_models/occupancy_single.txt", "model{
  # Priors for detection probability  
    for(o in 1:n_observer){
      alpha_0[o] ~ dnorm(mu_p, tau_p)
    }
    mu_p ~ dnorm(0, 0.00001)
    tau_p ~ dgamma(0.001,0.001) 
    for(n in 1:n_pred_alpha){
      alpha_1[n] ~ dnorm(0, 0.00001)
    }

  # Priors for occupancy probability
    for(t in 1:n_year){
      beta_0[t] ~ dnorm(mu_psi, tau_psi_year)
    }
    mu_psi ~ dnorm(0, 0.00001)
    tau_psi_year ~ dgamma(0.001,0.001) 
    for(m in 1:n_pred_beta){
      beta_1[m] ~ dnorm(0, 0.00001)
    }        
    
  # Likelihood
    for(i in 1:n_sites){ 
      z[i] ~ dbern(psi[i]) # True Occupancy
      logit(psi[i]) <- beta_0[year[i]] + inprod(beta_1, x_occ[i,]) # Probability of true occupancy
      for(j in 1:n_visits){
        y[i,j] ~ dbin(mu[i,j], 2) # Probability of observation in two visits
        mu[i,j] = z[i] * p[i,j]
        logit(p[i,j]) <- alpha_0[observer[i]] + # Random effect of obsever
                         alpha_1[1] * elev[i]  + 
                         alpha_1[2] * elev_sq[i] +
                         alpha_1[3] * day[i,j] + 
                         alpha_1[4] * day_sq[i,j]
      }
    }
  }")

# ------------------------------------------------------------------------------------- #
#### Fit model, loop over all species ####
cl = makeCluster(6)
registerDoParallel(cl)

foreach(spec = species_final$spec_id[27:nrow(species_final)], .packages = c("tidyverse", "rjags")) %dopar% {
  # cat(species_final$species[species_final$spec_id == paste(spec)], "\n") 
  
  #### Prepare Data ####
  # Observations
  y = observations_completed %>% 
    filter(spec_id == as.character(spec)) %>% 
    dplyr::select(preabs1:preabs7, site_id, year) %>% 
    drop_na(preabs1:preabs7) %>% 
    droplevels() 
  
  if(all(summarize(y, across(preabs1:preabs7, sum)) == 0)){
    cat("No observations. Skipping species.")
    next()
  }
  
  # Predictors state model
  x_occ = sample_sites@data[paste(y$site_id),] %>% 
    dplyr::select(ddeg0, bio_12, rad, asp, slp) %>%  
    mutate_all(scale) %>% # center and rescale
    mutate_all(list(sq = ~.*.)) %>%  # add quadratic terms 
    as.matrix()
  year = as.factor(y$year)

  # Predictors detection model
  x_det_raw = observations_completed %>% 
    filter(spec_id == spec) %>% 
    dplyr::select(observer, elevation, preabs1:preabs7, day1:day7) %>% 
    drop_na(preabs1:preabs7)
  
  day = matrix(apply(select(x_det_raw, day1:day7), 2, as.numeric), ncol = 7)
  day[is.na(day)] = mean(day, na.rm = T)  
  day = (day - mean(day)) / sd(day)
  day_sq = day*day
  elev = as.vector(scale(as.numeric(x_det_raw$elevation)))
  elev_sq = elev*elev
  observer = as.factor(x_det_raw$observer)
 
  # Final clean up
  y = y %>% select(-site_id, -year) %>% as.matrix()
  
  # Bundle data
  jags_data = list(y = y, n_sites = nrow(y), n_visits = ncol(y), n_pred_alpha = 4, n_pred_beta = ncol(x_occ), n_observer = nlevels(observer), 
                   n_year = nlevels(year), x_occ = x_occ, day = day, day_sq = day_sq, elev = elev, elev_sq = elev_sq, observer = observer, year = year)
  
  # Initial values
  jags_inits = function(){
    list(z = rep(1, nrow(y)), mu_p = runif(1), mu_psi = runif(1))
  }
  
  # Parameters to be monitored
  params = c("mu_p", "mu_psi", "alpha_1", "beta_1")
  
  # Run model 
  jags_model = jags.model(file = "Butterfly_detection/jags_models/occupancy_single.txt", 
                          data = jags_data, 
                          inits = jags_inits, 
                          n.chains = 2, n.adapt = 1000)
  update(jags_model, 2000)
  jags_samples = coda.samples(jags_model, params, n.iter = 5000, thin = 20)
  save(jags_samples, file = paste0("Data/model_fits/", spec, "_MCMC.RData"))
}

stopCluster(cl)
