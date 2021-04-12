library("tidyverse")
library("ubms")
library("unmarked")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")

rm(list=ls())
load("Data/observations_completed.RData")
load("Data/sample_sites.RData")
load("Data/species_final.RData")

# ------------------------------------------------------------------------------------- #
spec = sample(species_final$spec_id, 1)

#### Prepare Data ####
# Observations
y = observations_completed %>% 
  filter(spec_id == as.character(spec)) %>% 
  dplyr::select(preabs1:preabs7, site_id, year) %>% 
  pivot_longer(cols = preabs1:preabs7, names_to = "visit", values_to = "preabs") %>% 
  mutate(preabs1 = ifelse(preabs == 1 | preabs == 2, 1, 0),
         preabs2 = ifelse(preabs == 2, 1, 0)) %>% 
  select(-visit, -preabs) %>% 
  droplevels() 

# Predictors state model
x_occ = sample_sites@data[paste(y$site_id),] %>% 
  dplyr::select(ddeg0, bio_12, rad, asp, slp) %>%  
  mutate_all(scale) %>% # center and rescale
  mutate_all(list(sq = ~.*.)) %>%  # add quadratic terms 
  mutate(year = as.factor(y$year))

# Predictors detection model
x_det_raw = observations_completed %>% 
  filter(spec_id == as.character(spec)) %>% 
  dplyr::select(observer, elevation, preabs1:preabs7, day1:day7) %>% 
  pivot_longer(cols = day1:day7, names_to = "visit", values_to = "day") 

day = matrix(scale(x_det_raw$day), nrow = nrow(y), ncol = 2)
day[is.na(day)] = mean(day, na.rm = T)  
day_sq = day*day
elev = matrix(scale(x_det_raw$elevation), nrow = nrow(y), ncol = 2)
elev_sq = elev*elev
observer = matrix(as.factor(x_det_raw$observer), nrow = nrow(y), ncol = 2)

x_det = list(day = data.frame(day), elev = data.frame(elev), day_sq = data.frame(day_sq), elev_sq = data.frame(elev_sq), observer = data.frame(observer))

# Final clean up
y = y %>% select(-site_id, -year) %>% as.matrix()

data_occu = unmarkedFrameOccu(y = y, siteCovs = x_occ, obsCovs = x_det)
model_fit = stan_occu(formula = ~day+day_sq+elev+elev_sq+(1|observer)
                                ~.-year + (1|year),  
                      data = data_occu, iter = 500)

model_fit
save(model_fit, file = "Data/model_fits/31111_MCMC_stan.RData")