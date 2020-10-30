library("tidyverse")
library("abind")
library("rjags")
library("coda")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")

rm(list=ls())
load("Data/observations_completed.RData")
load("Data/sample_sites.RData")
load("Data/species_final.RData")
load("Data/traits_final.RData")

# ------------------------------------------------------------------------------------- #
#### Define JAGS model ####
cat(file="Butterfly_detection/jags_models/occupancy_multi.txt", "model{
  #          Priors            #
  # -------------------------- #
  # Global trait effects
  alpha_brightness ~ dnorm(0, 0.0001)
  alpha_saturation ~ dnorm(0, 0.0001)
  
  # Global random observer effect
  for(o in 1:n_observer){
    alpha_null_obs[o] ~ dnorm(mu_alpha_null_obs, tau_alpha_null_obs)
  }
  mu_alpha_null_obs ~ dnorm(0, 0.0001)
  tau_alpha_null_obs ~ dgamma(0.01,0.01)
  
  # Species-level effects
  for (i in 1:n_spec){
    ## 1. Detection model intercept ##
    alpha_null_sp[i] ~ dnorm(mu_alpha_null_sp[i], tau_alpha_null_sp[i])
    mu_alpha_null_sp[i] ~ dnorm(0, 0.0001)
    tau_alpha_null_sp[i] ~ dgamma(0.01,0.01)
    
    ## 2. Detection model slope ##
    alpha_day[i] ~ dnorm(0, 0.0001)
    alpha_day_sq[i] ~ dnorm(0, 0.0001)
    alpha_elev[i] ~ dnorm(0, 0.0001)
    alpha_elev_sq[i] ~ dnorm(0, 0.0001)

    ## 3. State model intercept ##
    for(t in 1:n_year){
      beta_null[i,t] ~ dnorm(mu_beta_null[i], tau_beta_null[i])
    }
    mu_beta_null[i] ~ dnorm(0, 0.001)
    tau_beta_null[i] ~ dgamma(0.01,0.01)
    
    ## 4. State model slope ##
    for(p in 1:n_pred_occ) {
      beta_coef[i,p] ~ dnorm(mu_beta_coef[i], tau_beta_coef[i])
    }
    mu_beta_coef[i] ~ dnorm(0, 0.0001)
    tau_beta_coef[i] ~ dgamma(0.01,0.01)
  }
  
  #        Likelihood          #
  # -------------------------- #
  for(i in 1:n_spec){
    for(j in 1:n_sites){
      z[i,j] ~ dbern(psi[i,j])
      logit(psi[i,j]) = beta_null[i,year[j]] + inprod(beta_coef[i,], x_state[j,])
      for(k in 1:n_visits){
        y[i,j,k] ~ dbin(mu[i,j,k], 2)                           # Probability of observation in two visits (secondary samples)
        mu[i,j,k] = z[i,j] * p[i,j,k]
        logit(p[i,j,k]) = alpha_null_sp[i] +
                          alpha_null_obs[observer[i]] +
                          alpha_brightness * brightness[i] +
                          alpha_saturation * saturation[i] +
                          alpha_elev[i] * elev[j] +
                          alpha_elev_sq[i] * elev_sq[j] +
                          alpha_day[i] * day[j,k] +
                          alpha_day_sq[i] * day_sq[j,k]
                          
        # Chi-square test statistic and posterior predictive check
        #chi2[i,j,k] = pow((y[i,j,k] - 2*mu[i,j,k]),2) / (sqrt(2*mu[i,j,k]) + 0.00001) # observed
        #y_pr[i,j,k] ~ dbin(mu[i,j,k], 2)
        #chi2_pr[i,j,k] = pow((y_pr[i,j,k] - 2*mu[i,j,k]),2) / (sqrt(2*mu[i,j,k]) + 0.00001) # expected
      }
    }
  }
  
  #     Derived Quantities     #
  # -------------------------- #
}")

# ------------------------------------------------------------------------------------- #
#### Prepare Data ####
# Observations
y = abind(lapply(species_final$spec_id, FUN = function(spec){
  spec_site = observations_completed %>% 
    dplyr::filter(spec_id == spec) %>%
    arrange(site_id, year) %>%
    mutate(row_names = paste0(site_id, "_", year)) %>% 
    column_to_rownames("row_names") %>% 
    select(preabs1:preabs7) 
}), along = 3, new.names = species_final$spec_id)
y = aperm(y, c(3,1,2))
# Predictors state model
site_year = str_split_fixed(dimnames(y)[[2]], pattern = "_", n = 2)  
colnames(site_year) = c("site_id", "year")

x_state = site_year %>% 
  as_tibble() %>% 
  left_join(rownames_to_column(sample_sites@data, "site_id")) %>% 
  dplyr::select(ddeg0, bio_12, rad, asp, slp) %>%  
  mutate_all(scale) %>% # center and rescale
  mutate_all(list(sq = ~.*.)) %>%  # add quadratic terms 
  as.matrix()
year = as.factor(site_year[,"year"])

# Predictors detection model
x_det_site = observations_completed %>% 
  drop_na(preabs1:preabs7) %>% 
  dplyr::select(site_id, year, observer, elevation, day1:day7) %>% 
  arrange(site_id, year) %>% 
  distinct()

x_det_spec = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  select(saturation, brightness)
  
day = matrix(apply(select(x_det_site, day1:day7), 2, as.numeric), ncol = 7)
day[is.na(day)] = mean(day, na.rm = T)  
day = (day - mean(day)) / sd(day)
day_sq = day*day
elev = as.vector(scale(as.numeric(x_det_site$elevation)))
elev_sq = elev*elev
observer = as.factor(x_det_site$observer)

# Bundle data
jags_data = list(n_spec = dim(y)[1], n_sites = dim(y)[2], n_visits = dim(y)[3], 
                 n_pred_occ = ncol(x_state), n_observer = nlevels(observer), n_year = nlevels(year), 
                 y = y, x_state = x_state, year = year,
                 day = day, day_sq = day_sq, elev = elev, elev_sq = elev_sq, observer = observer,
                 saturation = x_det_spec$saturation, brightness = x_det_spec$brightness)

# Initial values
jags_inits = function(){
 list(z = array(1, dim(y)[c(1,2)]))
}

# Parameters to be monitored
params = c("alpha_brightness", "alpha_saturation", "alpha_null_sp", "mu_null_obs", "alpha_elev", "alpha_elev_sq", "alpha_day", "alpha_day_sq",
           "chi2", "chi2_pr")

# Fit model 
jags_model = jags.model(file = "Butterfly_detection/jags_models/occupancy_multi.txt", 
                        data = jags_data, 
                        inits = jags_inits, 
                        n.chains = 1, n.adapt = 100)
update(jags_model, 500) # Burn-In
jags_samples = coda.samples(jags_model, params, n.iter = 1000, thin = 5) # Sampling
save(jags_samples, file = paste0("Data/models_test/multi_MCMC_1000it.RData"))
