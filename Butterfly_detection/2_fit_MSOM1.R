library(tidyverse)
library(abind)
library(runjags)
library(rjags)

setwd("~/Butterfly_project")

rm(list=ls())
load("Data/observations_completed.RData")
load("Data/sample_sites.RData")
load("Data/species_final.RData")

# ------------------------------------------------------------------------------------ #
#### Prepare Data ####
observations_completed = observations_completed %>% 
  dplyr::filter(year >= 2010, year <= 2014)

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

# ------------------------------------------------------------------------------------- #
#### Fit models ####
n_chains = 4
n_adapt = 1000 
n_burnin = 0
n_sample = 0
n_thin = 1
n_cores = 4

# MSOM 1
MSOM1_data = list(n_spec = dim(y)[1], n_sites = dim(y)[2], n_visits = dim(y)[3],
                  y = y, x_state = x_state,
                  x_det_day = x_det_day, x_det_elev = x_det_elev,
                  n_pred_occ = ncol(x_state), n_pred_det_env = 4)

MSOM1_inits = function(){list(z = array(1, dim(y)))} # always start with z = 1

MSOM1_params = c("alpha_null", "alpha_coef_env", "beta_null", "beta_coef", "deviance")

# The following code iteratively samples 1000 MCMC draws and saves the temporary result as separate files. Previous  samples are not deleted, so the 
# model objects become larger with every iteration. We did this because we needed some insight into the blackbox that is running JAGS models from R, which was
# characterized by frequent crashes, memory leaks or and other unexpected behaviour from the different wrapper packages. 
MSOM1_samples = run.jags(model = "Butterfly_detection/jags_models/MSOM_1.txt",
                         monitor = MSOM1_params, data = MSOM1_data, inits = MSOM1_inits,
                         n.chains = n_chains, burnin = n_burnin, sample = n_sample, adapt = n_adapt, thin = n_thin,
                         jags.refresh = 30, summarise = TRUE, method = "rjparallel")

saveRDS(MSOM1_samples, file = "//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM1/MSOM1_1.RDS")

for(i in 2:15){
  MSOM1 = readRDS(paste0("//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM1/MSOM1_", i-1,".RDS"))
  MSOM1_ext = extend.jags(MSOM1, sample = 1000, burnin = 0, adapt = 0)
  saveRDS(MSOM1_ext, file = paste0("//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM1/MSOM1_", i,".RDS"))
}