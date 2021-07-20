library("tidyverse")
library("coda")       # mcmc/coda objects
library("runjags")    # JAGS wrapper
library("bayesplot")


setwd("~/ownCloud/Projects/Berlin/06_Butterfly_detection/")
rm(list=ls())

load("Data/species_final.RData")
load("Data/traits_final.RData")
load("Data/MSOM1_summary.RData")
load("Data/MSOM2_summary.RData")
load("Data/MSOM3_summary.RData")

# --------------------------------------------------------- #
#                    Compare model fit                   ####
# --------------------------------------------------------- #

# DIC analyses here

# --------------------------------------------------------- #
#                Merge model summaries                   ####
# --------------------------------------------------------- #
MSOM1_df = MSOM1_summary$statistics %>% 
  as_tibble(rownames = "parameter") %>%
  mutate(model = "MSOM1") %>% 
  relocate(model, parameter)

MSOM2_df = MSOM2_summary$statistics %>% 
  as_tibble(rownames = "parameter") %>%
  mutate(model = "MSOM2") %>% 
  relocate(model, parameter)

MSOM3_df = MSOM3_summary$statistics %>% 
  as_tibble(rownames = "parameter") %>%
  mutate(model = "MSOM3") %>% 
  relocate(model, parameter)

model_summary = bind_rows(MSOM1_df, MSOM2_df, MSOM3_df) %>% 
  mutate(row_id = str_extract(parameter, "\\d{1,}"))


# --------------------------------------------------------- #
#                Explore trait effects                   ####
# --------------------------------------------------------- #
antilogit = function(x){exp(x) / (1+exp(x))}    # Convert from logit scale back to probabilities

#### MSOM1 ####
p_det_1 = model_summary %>% 
  filter(model == "MSOM1", str_detect(parameter, "alpha_null")) %>% 
  pull(Mean) %>% 
  antilogit()

plot(traits_final$Vol_min, p_det_1)      # Voltinism
plot(traits_final$WIn, p_det_1)          # Wing index
plot(traits_final$HSI, p_det_1)          # Host specificity
plot(traits_final$FMo_Average, p_det_1)  # Flight period
plot(traits_final$mean_sat_top, p_det_1)            # Saturation top
plot(traits_final$mean_sat_bottom, p_det_1)         # Saturation bottom
plot(traits_final$mean_lgt_top, p_det_1)            # Lightness top
plot(traits_final$mean_lgt_bottom, p_det_1)         # Lightness bottom
boxplot(p_det_1 ~ traits_final$main_color_top)      # Color top
boxplot(p_det_1 ~ traits_final$main_color_bottom)   # Color bottom


#### MSOM2 ####
# Look at model coefficients directly
traits_design_matrix = readRDS("Data/traits_design_matrix.RDS")
trait_effects_MSOM2 = model_summary %>% 
  filter(model == "MSOM2", str_detect(parameter, "alpha_coef_trait")) %>% 
  mutate(trait_name = colnames(traits_design_matrix))
  
ggplot(trait_effects_MSOM2) +
  geom_point(aes(y = Mean, x = trait_name)) +
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD, x = trait_name), width = 0.2) +
  coord_flip() +
  theme_bw()

#### MSOM3 ####
p_det_3 = model_summary %>% 
  filter(model == "MSOM3", str_detect(parameter, "^alpha_null_l1")) %>% 
  pull(Mean) %>% 
  antilogit()

plot(traits_final$Vol_min, p_det_3)      # Voltinism
plot(traits_final$WIn, p_det_3)          # Wing index
plot(traits_final$HSI, p_det_3)          # Host specificity
plot(traits_final$FMo_Average, p_det_3)  # Flight period
plot(traits_final$mean_sat_top, p_det_3)            # Saturation top
plot(traits_final$mean_sat_bottom, p_det_3)         # Saturation bottom
plot(traits_final$mean_lgt_top, p_det_3)            # Lightness top
plot(traits_final$mean_lgt_bottom, p_det_3)         # Lightness bottom
boxplot(p_det_3 ~ traits_final$main_color_top)      # Color top
boxplot(p_det_3 ~ traits_final$main_color_bottom)   # Color bottom
