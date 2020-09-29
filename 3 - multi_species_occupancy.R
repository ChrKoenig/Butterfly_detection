library("tidyverse")
library("abind")
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
#### Prepare Data ####
tmp_year = 2009 

# Final species selection (observed + traits available)
species_final = tibble(species = intersect(drop_na(traits_final)$species, butterflies$species)) %>% # species with survey data and complete traits
  inner_join(butterflies %>% filter(year == tmp_year) %>% select(spec_id, species) %>% distinct())

# Observations
y_df = butterflies %>% 
  filter(year == tmp_year & spec_id %in% species_final$spec_id) %>% 
  mutate_at(vars(contains("Day")), as.numeric) %>%
  rename_at(.vars = vars(contains("day")), .funs = funs(str_replace(. , "day", "day_"))) %>%
  rename_at(.vars = vars(contains("preabs")), .funs = funs(str_replace(. , "preabs", "preabs_"))) %>%
  pivot_longer(cols = day_1:preabs_7, names_to = c(".value", "visit"), names_sep = "_") %>% 
  mutate("site_id" = as.factor(site_id)) %>% 
  complete(site_id, visit, spec_id)

y = abind(  # i x j x k array of collapsed species observations (0/1/2)
  lapply(species_final$spec_id, function(x){
    tmp = y_df %>% filter(spec_id == x) %>% 
      pivot_wider(id_cols = site_id, names_from = visit, values_from = preabs) %>% 
      replace(is.na(.), 0) %>% 
      column_to_rownames("site_id")
  }), along = 3, new.names = as.character(species_final$spec_id) # set dimnames
) 

# Predictors state (occupancy) model
x_occ = sample_sites@data[dimnames(y)[[1]],] %>% # i x n_pred_vars matrix of site covariates
  dplyr::select(bio_1, bio_4, bio_12, bio_15) %>%  # TODO: Log-Transform?
  mutate_all(scale) %>% # center and rescale
  mutate_all(list(sq = ~.*.)) %>%  # add quadratic terms 
  as.matrix()

# Predictors detection model
day = y_df %>% # i x k matrix of visit dates
  select(site_id, visit, day) %>% 
  drop_na() %>% 
  distinct() %>%
  mutate(day = as.vector(scale(day))) %>% 
  pivot_wider(id_cols = site_id, names_from = visit, values_from = day) %>% 
  column_to_rownames("site_id") %>% # match column-order with observations
  select(dimnames(y)[[2]]) %>% 
  as.matrix() %>% 
  replace(is.na(.), 0) 
day_sq = day*day

elev = sample_sites@data[dimnames(y)[[1]],] %>% # i values for elevation
  dplyr::select(elevation) %>% 
  mutate_all(scale) %>% 
  pull(elevation) %>% 
  as.vector()
elev_sq = elev*elev # i values for elevation²

observer = y_df %>% # i x k matrix of observer identities
  select(site_id, yearXrepl, observer) %>% 
  drop_na() %>% 
  distinct() %>% 
  mutate(observer = as.factor(observer)) %>% 
  pivot_wider(id_cols = site_id, names_from = yearXrepl, values_from = observer) %>% 
  column_to_rownames("site_id") %>% # match column-order with observations
  select(dimnames(y)[[2]]) %>% 
  as.matrix()
observer = observer[dimnames(y)[[1]],] # match row-order with observations

traits = traits_final %>% 
  inner_join(species_final) %>% 
  select(spec_id, generations, body_area, RGB_top, RGB_bottom) %>% # TODO: binary-code categorical traits?
  column_to_rownames("spec_id") %>% 
  as.matrix()
traits =  traits[dimnames(y)[[3]],]


# Bundle data
jags_data = list(y = y, n_sites = dim(y)[1], n_visits = dim(y)[2], n_spec = dim(y)[2], 
                 n_pred_alpha = 4, n_pred_beta = ncol(x_occ), n_observer = length(unique(y_df$observer)), 
                 x_occ = x_occ, day = day, day_sq = day_sq, elev = elev, elev_sq = elev_sq, traits = traits, observer = observer)

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