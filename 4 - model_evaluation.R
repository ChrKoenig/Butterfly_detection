library("tidyverse")
library("rjags")
library("coda")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list=ls())

# ------------------------------------------------------------------------------------- #
#### Convergence ####
load(sample(list.files("Data/model_fits/", full.names = T), 1)) # sample model for random species

effectiveSize(jags_samples)
plot(jags_samples)
gelman.plot(jags_samples)

#### Model fit ####

#### Response curves ####

