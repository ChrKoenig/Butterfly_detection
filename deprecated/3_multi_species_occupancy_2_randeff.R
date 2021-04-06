library("tidyverse")
library("abind")
library("jagsUI")

setwd("~/Butterfly_project")

rm(list=ls())
load("Data/observations_completed.RData")
load("Data/sample_sites.RData")
load("Data/species_final.RData")
load("Data/traits_final.RData")

# ------------------------------------------------------------------------------------- #
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
  model.matrix(~ 0 + main_color_top, data = .)
color_bottom = species_final %>% 
  left_join(traits_final, by = "species") %>%
  mutate(main_color_bottom = fct_relevel(main_color_bottom, "none", "red", "orange", "yellow", "green")) %>% 
  model.matrix(~ 0 + main_color_bottom, data = .)

x_det_traits = cbind(traits_num, color_top, color_bottom)

# ------------------------------------------------------------------------------------- #
#### Fit models ####
n_chains = 8
n_adapt = 300
n_burnin = 5000
n_iter = 7000
n_thin = 140
n_cores = 8

# MSOM 3
MSOM3_data = list(n_spec = dim(y)[1], n_sites = dim(y)[2], n_visits = dim(y)[3], 
                  y = y, x_state = x_state,
                  x_det_day = x_det_day, x_det_elev = x_det_elev, x_det_traits = x_det_traits,
                  year = year, n_year = n_distinct(year), observer = observer, n_observer = n_distinct(observer),
                  n_pred_occ = ncol(x_state), n_pred_det_env = 4, n_pred_det_traits = ncol(x_det_traits))

MSOM3_inits = function(){list(z = array(1, dim(y)))} # always start with z = 1

MSOM3_params = c("mu_alpha_null", "alpha_null", "alpha_coef_env", "alpha_coef_traits","beta_null", "beta_coef")

MSOM3_samples = jags(data = MSOM3_data, 
                     inits = MSOM3_inits, 
                     parameters.to.save = MSOM3_params, 
                     model.file = "Butterfly_detection/jags_models/MSOM_3.txt",
                     n.chains = n_chains, n.adapt = n_adapt, n.burnin = n_burnin, n.iter = n_iter, n.thin = n_thin, 
                     parallel = T, n.cores = n_cores, DIC = T)
saveRDS(MSOM3_samples, file = "Data/models_fit/MSOM/MSOM_3.RDS")