library("tidyverse")
library("abind")
library("jagsUI")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")

rm(list=ls())
load("Data/observations_completed.RData")
load("Data/sample_sites.RData")
load("Data/species_final.RData")
load("Data/traits_final.RData")

# ------------------------------------------------------------------------------------- #
#                        COMPLEX MULTI SPECIES OCCUPANCY MODEL
# Random effect of year on occupancy intercept
# Additional term in detection submodel representing observer-related variation
# Effects of traits on detection intercept
# Community hyperpriors on species-level parameters
# --------------------------------------------------------------------------------------#
#### Define JAGS model ####
cat(file="Butterfly_detection/jags_models/MSOM_3.txt", "model{
  #          Priors            #
  # -------------------------- #
  
  ## 1. Detection model ##
  # Intercept with random observer effect
  for(o in 1:n_observer){
    alpha_null[o] ~ dnorm(mu_alpha_null, tau_alpha_null)
  }
  mu_alpha_null ~ dnorm(0, 0.0001) T(-12,12)
  tau_alpha_null ~ dgamma(0.01,0.01)
  
  # Trait effects on intercept
  for(p in 1:n_pred_det_traits){
    alpha_coef_traits[p] ~ dnorm(0, 0.0001) T(-12,12)
  }
  
  # Hyperpriors for slopes
  for(p in 1:n_pred_det_env) {
    mu_alpha_coef_env[p] ~ dnorm(0, 0.0001) T(-12,12)
    tau_alpha_coef_env[p] ~ dgamma(0.01,0.01)
  }
  
  # Species-level slopes
  for (i in 1:n_spec){
    for(p in 1:n_pred_det_env) {
      alpha_coef_env[i,p] ~ dnorm(mu_alpha_coef_env[p], tau_alpha_coef_env[p])
    }
  }
  
  ## 2. State model ##
  # Hyperpriors for intercept and slopes
  mu_beta_null ~ dnorm(0, 0.0001) T(-12,12)
  tau_beta_null ~ dgamma(0.01,0.01) 
  for(p in 1:n_pred_occ) {
    mu_beta_coef[p] ~ dnorm(0, 0.0001) T(-12,12)
    tau_beta_coef[p] ~ dgamma(0.01,0.01)
  }
  
  for(i in 1:n_spec)
    for(t in 1:n_year){
      beta_null[i,t] ~ dnorm(mu_beta_null, tau_beta_null)
    }
    for(p in 1:n_pred_occ) {
      beta_coef[i,p] ~ dnorm(mu_beta_coef[p], tau_beta_coef[p])
    }
  }

  #        Likelihood          #
  # -------------------------- #
  for(i in 1:n_spec){
    for(j in 1:n_sites){
      logit(psi[i,j]) = beta_null[i,year[j]] + inprod(beta_coef[i,], x_state[j,])
      
      for(k in 1:n_visits){
        # Probability of observation in two visits (secondary samples)
        y[i,j,k] ~ dbin(mu[i,j,k], 2)
        mu[i,j,k] = z[i,j] * p[i,j,k]
        
        # occupancy state during visit k
        z[i,j,k] ~ dbern(psi[i,j])
        
        # detection probability
        logit(p[i,j,k]) = alpha_null[observer[i]] +
                          inprod(alpha_coef_traits, x_det_traits[i,]) +
                          alpha_coef_env[i,1] * x_det_elev[j] + 
                          alpha_coef_env[i,2] * pow(x_det_elev[j], 2) +
                          alpha_coef_env[i,3] * x_det_day[j,k] +
                          alpha_coef_env[i,4] * pow(x_det_day[j,k], 2)
        
        # posterior predictive distribution
        y_pr[i,j,k] ~ dbin(mu[i,j,k], 2)
        
        # Chi-square test statistic and posterior predictive check
        chi2[i,j,k] <- pow((y[i,j,k] - 2*mu[i,j,k]),2) / (sqrt(2*mu[i,j,k]) + 0.0001) # observed
        chi2_pr[i,j,k] <- pow((y_pr[i,j,k] - 2*mu[i,j,k]),2) / (sqrt(2*mu[i,j,k]) + 0.0001) # expected
      }
    }
  }
  
  #       Derived measures     #
  # -------------------------- #
  fit = sum(chi2)
  fit_pr = sum(chi2_pr)
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
    dplyr::select(preabs1:preabs7) 
}), along = 3, new.names = species_final$spec_id)
y = aperm(y, c(3,1,2)) # reorder to species, site, visit

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
year = as.integer(factor(site_year[,"year"]))

# Predictors detection model
x_det = observations_completed %>% 
  drop_na(preabs1:preabs7) %>% 
  dplyr::select(site_id, year, observer, elevation, day1:day7) %>% 
  arrange(site_id, year) %>% 
  distinct()
day = matrix(apply(dplyr::select(x_det, day1:day7), 2, as.numeric), ncol = 7)
day[is.na(day)] = mean(day, na.rm = T)  
x_det_day = (day - mean(day)) / sd(day)
x_det_elev = as.vector(scale(as.numeric(x_det$elevation)))
observer = as.integer(factor(x_det$observer))

traits_num = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  dplyr::select(Vol_min, WIn, HSI, FMo_Average, mean_sat_top, mean_sat_bottom, mean_lgt_top, mean_lgt_bottom)
color_top = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  mutate(main_color_top = fct_relevel(main_color_top, "none", "red", "orange", "yellow", "green", "blue")) %>% 
  model.matrix(~ main_color_top, data = .)
color_bottom = species_final %>% 
  left_join(traits_final, by = "species") %>%
  mutate(main_color_bottom = fct_relevel(main_color_bottom, "none", "red", "orange", "yellow", "green")) %>% 
  model.matrix(~ main_color_bottom, data = .)

x_det_traits = cbind(traits_num, color_top, color_bottom)

# Bundle data
jags_data = list(n_spec = dim(y)[1], n_sites = dim(y)[2], n_visits = dim(y)[3], 
                 y = y, x_state = x_state,
                 x_det_day = x_det_day, x_det_elev = x_det_elev, x_det_traits = x_det_traits,
                 year = year, n_year = n_distinct(year), observer = observer, n_observer = n_distinct(observer),
                 n_pred_occ = ncol(x_state), n_pred_det_env = 4, n_pred_det_traits = ncol(x_det_traits))

# Initial values
jags_inits = function(){
  list(z = array(1, dim(y))) # always start with z = 1
}

# Parameters to be monitored
jags_params = c("alpha_null", "alpha_coef_env", "alpha_coef_traits","beta_null", "beta_coef", "fit", "fit_pr")

# Fit model 
jags_samples = jags(data = jags_data, 
                    inits = jags_inits, 
                    parameters.to.save = jags_params, 
                    model.file = "Butterfly_detection/jags_models/MSOM_3.txt",
                    n.chains = 1, n.adapt = 10, n.burnin = 0, n.iter = 50, n.thin = 1, 
                    parallel = F, n.cores = 1, DIC = T)
saveRDS(jags_samples, file = "Data/models_fit/MSOM/MSOM_3_run1.RDS")