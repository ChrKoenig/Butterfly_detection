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

traits_num = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  dplyr::select(Vol_min, WIn, HSI, FMo_Average, mean_sat_top, mean_sat_bottom, mean_lgt_top, mean_lgt_bottom)
main_color = species_final %>% 
  left_join(traits_final, by = "species") %>% 
  mutate(main_color = ifelse(main_color_top == "none", main_color_bottom, main_color_top),
         main_color = fct_relevel(main_color, "none", "red", "orange", "yellow", "green", "blue")) %>% 
  model.matrix(~ main_color, data = .)
main_color = main_color[,-1] # remove intercept (i.e."none") column; will be estimated as separate param in JAGS 

x_det_traits = cbind(traits_num, main_color)
# qr(x_det_traits) # --> Full rank

# ------------------------------------------------------------------------------------- #
#### Fit models ####
n_chains = 4
n_adapt = 300
n_burnin = 4000
n_sample = 8000
n_thin = 80
n_cores = 4

# MSOM 2
MSOM2_data = list(n_spec = dim(y)[1], n_sites = dim(y)[2], n_visits = dim(y)[3], 
                  y = y, x_state = x_state,
                  x_det_day = x_det_day, x_det_elev = x_det_elev, x_det_traits = x_det_traits,
                  n_pred_occ = ncol(x_state), n_pred_det_env = 4, n_pred_det_traits = ncol(x_det_traits))

MSOM2_inits = function(){list(z = array(1, dim(y)))} # always start with z = 1

MSOM2_params = c("alpha_null", "alpha_coef_env", "alpha_coef_traits", "beta_null", "beta_coef")

MSOM2_samples = run.jags(model = "Butterfly_detection/jags_models/MSOM_2.txt", 
                         monitor = MSOM2_params, data = MSOM2_data, inits = MSOM2_inits, 
                         n.chains = n_chains, burnin = n_burnin, sample = n_sample, adapt = n_adapt, thin = n_thin, 
                         jags.refresh = 30, method = "rjparallel")

saveRDS(MSOM2_samples, file = "Data/models_fit/MSOM_2_rjags.RDS")