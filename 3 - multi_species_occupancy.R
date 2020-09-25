library("tidyverse")
library("jagsUI")
library("sp")
library("unmarked")

rm(list=ls())
setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_occupancy/")

load("Data/butterflies.RData")
load("Data/sample_sites.RData")