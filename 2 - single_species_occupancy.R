library("tidyverse")
library("rjags")
library("coda")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())

load("Data/butterflies.RData")
load("Data/sample_sites.RData")
load("Data/traits_final.RData")

# ------------------------------------------------------------------------------------- #
#### Define JAGS model ####
cat(file="butterfly_detection/jags_models/occupancy_single.txt", "model{
  # Priors for detection probability  
    for(o in 1:n_observer){
      alpha_0[o] ~ dnorm(0, tau_p_obs)
    }
    tau_p_obs ~ dgamma(0.001,0.001) 
    for(n in 1:n_pred_alpha){
      alpha_1[n] ~ dnorm(0, 0.00001)
    }

  # Priors for occupancy probability´
    psi_mean ~ dunif(0,1) 
    beta_0 <- logit(psi_mean)
    for(m in 1:n_pred_beta){
      beta_1[m] ~ dnorm(0, 0.00001)
    }
    
  # Likelihood   TODO: RANDOM EFFECT OF FOR YEAR?
    for(i in 1:n_sites){ 
      z[i] ~ dbern(psi[i]) # True Occupancy
      logit(psi[i]) <- beta_0 + inprod(beta_1, x_occ[i,]) # Probability of true occupancy
        for(j in 1:n_visits){
        y[i,j] ~ dbin(mu[i,j], 2) # Probability of observation in two visits
        mu[i,j] = z[i] * p[i,j]
        logit(p[i,j]) <- alpha_0[observer[i]] + # Random effect of obsever
                         alpha_1[1] * elev[i]  + 
                         alpha_1[2] * elev_sq[i] +
                         alpha_1[3] * date[i,j] + 
                         alpha_1[4] * date_sq[i,j]
        }
    }
  }")

# ------------------------------------------------------------------------------------- #
#### Fit model ####
species_final = intersect(drop_na(traits_final)$species, butterflies$species)[1] # species with survey data and complete traits

for(spec in species_final){
  cat(spec)
  #### Prepare Data ####
  # Observations
  y = butterflies %>% 
    filter(species == spec) %>% 
    dplyr::select(preabs1:preabs7, site_id) %>% 
    drop_na(preabs1:preabs7) # TODO: Remove NA rows or replace with 0 and leave in?
  
  # Predictors state model
  x_occ = sample_sites@data[paste(y$site_id),] %>% 
    dplyr::select(bio_1, bio_4, bio_12, bio_15) %>%  # TODO: Log-Transform?
    mutate_all(scale) %>% # center and rescale
    mutate_all(list(sq = ~.*.)) %>%  # add quadratic terms 
    as.matrix()

  # Predictors detection model
  x_det_raw = butterflies %>% 
    filter(species == spec) %>% 
    select(observer, elevation, preabs1:preabs7, day1:day7) %>% 
    drop_na(preabs1:preabs7)
  
  date = matrix(apply(select(x_det_raw, day1:day7), 2, as.numeric), ncol = 7)
  date[is.na(date)] = mean(date, na.rm = T)  
  date = (date - mean(date)) / sd(date)
  date_sq = date*date
  elev = as.vector(scale(as.numeric(x_det_raw$elevation)))
  elev_sq = elev*elev
  observer = as.factor(x_det_raw$observer)
 
  y = y %>% select(-site_id) %>% as.matrix()
  
  # Bundle data
  jags_data = list(y = y, n_sites = nrow(y), n_visits = ncol(y), n_pred_alpha = 4, n_pred_beta = ncol(x_occ), n_observer = length(unique(observer)), 
                   x_occ = x_occ, date = date, date_sq = date_sq, elev = elev, elev_sq = elev_sq, observer = observer)
  
  # Initial values
  jags_inits = function(){
    list(z = rep(1, nrow(y)), psi_mean = runif(1))
  }
  
  # Parameters to be monitored
  params = c("alpha_0", "alpha_1", "beta_0", "beta_1")
  
  # Run model 
  jags_model = jags.model(file = "butterfly_detection/jags_models/occupancy_single.txt", 
                          data = jags_data, 
                          inits = jags_inits, 
                          n.chains = 2, n.adapt = 3000)
  jags_samples = coda.samples(jags_model, params, n.iter = 10000, thin = 50)
  save(jags_samples, file = paste0("Data/model_fits/", spec, " MCMC.RData"))
}

