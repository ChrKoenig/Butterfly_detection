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
#                        ENHANCED MULTI SPECIES OCCUPANCY MODEL
# NO effects of year on occupancy
# NO effects of observer on detection
# Effects of traits on detection
# Community hyperpriors on species-level parameters
# --------------------------------------------------------------------------------------#
#### Define JAGS model ####
cat(file="Butterfly_detection/jags_models/MSOM_2.txt", "model{
  #          Priors            #
  # -------------------------- #
  # Priors for global trait effects
  for(t in 1:n_traits_num){
    alpha_traits_num[t] ~ dnorm(0,0.0001) T(-12,12)
  }
  for(c in 1:n_color_top){
    alpha_color_top[c] ~ dnorm(0,0.0001) T(-12,12)
  }
  for(c in 1:n_color_bottom){
    alpha_color_bottom[c] ~ dnorm(0,0.0001) T(-12,12)
  }
 
  # Hyperpriors for detection model
  mu_alpha_null ~ dnorm(0, 0.0001) T(-12,12)
  tau_alpha_null ~ dgamma(0.01,0.01) 
  for(p in 1:n_pred_det) {
    mu_alpha_coef[p] ~ dnorm(0, 0.0001) T(-12,12)
    tau_alpha_coef[p] ~ dgamma(0.01,0.01)
  }
  
  # Hyperpriors for state model
  mu_beta_null ~ dnorm(0, 0.0001) T(-12,12)
  tau_beta_null ~ dgamma(0.01,0.01) 
  for(p in 1:n_pred_occ) {
    mu_beta_coef[p] ~ dnorm(0, 0.0001) T(-12,12)
    tau_beta_coef[p] ~ dgamma(0.01,0.01)
  }
  
  # Priors for species-level coefficients
  for (i in 1:n_spec){
    ## 1. Detection model ##
    alpha_null[i] ~ dnorm(mu_alpha_null, tau_alpha_null)
    for(p in 1:n_pred_det) {
      alpha_coef[i,p] ~ dnorm(mu_alpha_coef[p], tau_alpha_coef[p])
    }

    ## 2. State model ##
    beta_null[i] ~ dnorm(mu_beta_null, tau_beta_null)
    for(p in 1:n_pred_occ) {
      beta_coef[i,p] ~ dnorm(mu_beta_coef[p], tau_beta_coef[p])
    }
  }

  #        Likelihood          #
  # -------------------------- #
  for(i in 1:n_spec){
    for(j in 1:n_sites){
      z[i,j] ~ dbern(psi[i,j])
      logit(psi[i,j]) = beta_null[i] + inprod(beta_coef[i,], x_state[j,])
      for(k in 1:n_visits){
        y[i,j,k] ~ dbin(mu[i,j,k], 2)             # Probability of observation in two visits (secondary samples)
        mu[i,j,k] = z[i,j] * p[i,j,k]
        logit(p[i,j,k]) = alpha_null[i] +
                          inprod(alpha_traits_num, traits_num[i,]) +
                          inprod(alpha_color_top, color_top[i,]) +
                          inprod(alpha_color_bottom, color_bottom[i,]) +
                          alpha_coef[i,1] * elev[j] + 
                          alpha_coef[i,2] * elev_sq[j] +
                          alpha_coef[i,3] * day[j,k] +
                          alpha_coef[i,4] * day_sq[j,k]
        
        # posterior predictive distribution
        y_pr[i,j,k] ~ dbin(mu[i,j,k], 2)
        
        # Chi-square test statistic and posterior predictive check
        chi2[i,j,k] <- pow((y[i,j,k] - 2*mu[i,j,k]),2) / (sqrt(2*mu[i,j,k]) + 0.0001) # observed
        chi2_pr[i,j,k] <- pow((y_pr[i,j,k] - 2*mu[i,j,k]),2) / (sqrt(2*mu[i,j,k]) + 0.0001) # expected
      }
    }
  }
  # Derived measures
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

# Predictors detection model
x_det = observations_completed %>% 
  drop_na(preabs1:preabs7) %>% 
  dplyr::select(site_id, year, observer, elevation, day1:day7) %>% 
  arrange(site_id, year) %>% 
  distinct()
day = matrix(apply(dplyr::select(x_det, day1:day7), 2, as.numeric), ncol = 7)
day[is.na(day)] = mean(day, na.rm = T)  
day = (day - mean(day)) / sd(day)
day_sq = day*day
elev = as.vector(scale(as.numeric(x_det$elevation)))
elev_sq = elev*elev

traits_num = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  dplyr::select(Vol_min, WIn, HSI, FMo_Average, mean_sat_top, mean_sat_bottom, mean_lgt_top, mean_lgt_bottom)
color_top = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  model.matrix(~ main_color_top, data = .)
color_bottom = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  model.matrix(~ main_color_bottom, data = .)

# Bundle data
jags_data = list(n_spec = dim(y)[1], n_sites = dim(y)[2], n_visits = dim(y)[3], 
                 n_pred_occ = ncol(x_state), n_pred_det = 4, y = y, x_state = x_state,
                 day = day, day_sq = day_sq, elev = elev, elev_sq = elev_sq,
                 traits_num = traits_num, n_traits_num = ncol(traits_num), 
                 color_top = color_top, n_color_top = ncol(color_top), color_bottom = color_bottom, n_color_bottom = ncol(color_bottom))

# Initial values
jags_inits = function(){
  list(z = array(1, dim(y)[c(1,2)])) # always start with z = 1
}

# Parameters to be monitored
jags_params = c("alpha_null", "alpha_coef_env", "alpha_traits_num", "alpha_color_top",  "alpha_color_bottom",
                "beta_null", "beta_coef", "fit", "fit_pr")

# Fit model 
jags_samples = jags(data = jags_data, 
                    inits = jags_inits, 
                    parameters.to.save = jags_params, 
                    model.file = "Butterfly_detection/jags_models/MSOM_2.txt",
                    n.chains = 4, n.adapt = 500, n.burnin = 5000, n.iter = 15500, n.thin = 100, 
                    parallel = T, n.cores = 4, DIC = T)
saveRDS(jags_samples, file = "Data/models_fit/MSOM/MSOM_2_run1.RDS")