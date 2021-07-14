library("rjags") 

setwd("~/ownCloud/Projects/Berlin/06_Butterfly_detection/")
rm(list=ls())
load("Data/species_final.RData")

# ------------------------------------------------------- #
#              Sample DIC for fitted models            ####
# ------------------------------------------------------- #
MSOM1 = readRDS("~/Data/Butterfly_models_fit/MSOM1/MSOM1_15.RDS")
MSOM2 = readRDS("~/Data/Butterfly_models_fit/MSOM2/MSOM2_15.RDS")
MSOM3 = readRDS("~/Data/Butterfly_models_fit/MSOM3/MSOM3_15.RDS")

