##############################
# Static occupancy model with spatial autocovariates
# Author: mingjian
# Date: 22.10.2019
# Data: 1550_TF_Daten_blinded.csv
# Spatial model: CAR
##############################


setwd()
# library(plyr)
butterfly <- read.csv("butterfly.csv")


########## prepare data ##########
# normal recording
butterfly <- butterfly[butterfly$Aufnahmetyp=="Normalaufnahme_Z7", ]
# NA and ?
butterfly[butterfly=="?" | butterfly==""] <- NA
butterfly <- na.omit(butterfly)
# convert dates to numbers
butterfly[, 15:21] <- apply(butterfly[, 15:21], 2, as.numeric)


########## 2017 Aglais urticae ##########
butterfly_2017 <- butterfly[butterfly$Jahr==2017, ]
# count <- plyr::count(butterfly_2017, vars="NameTagfLat")
au_2017 <- butterfly_2017[butterfly_2017$NameTagfLat=="Aglais urticae", ]
nperiod=7

# sites
sites <- sort(unique(au_2017$aldStao_)) # code of sites
nsites <- length(sites)

# detection : xaug 1)sites 2)periods
xaug <- as.matrix(data.frame(list(
  D1=tapply(au_2017$PreAbs1, au_2017$aldStao_, function(x) x[1]),
  D2=tapply(au_2017$PreAbs2, au_2017$aldStao_, function(x) x[1]),
  D3=tapply(au_2017$PreAbs3, au_2017$aldStao_, function(x) x[1]),
  D4=tapply(au_2017$PreAbs4, au_2017$aldStao_, function(x) x[1]),
  D5=tapply(au_2017$PreAbs5, au_2017$aldStao_, function(x) x[1]),
  D6=tapply(au_2017$PreAbs6, au_2017$aldStao_, function(x) x[1]),
  D7=tapply(au_2017$PreAbs7, au_2017$aldStao_, function(x) x[1])
)))

xaug[xaug==2] <- 1
xaug[xaug==3] <- 2
xaug[xaug==-1] <- NA
xaug[xaug==6] <- 2

# elevation
elevation <- as.numeric(tapply(au_2017$Hoehe, au_2017$aldStao_, function(x) x[1]))
elevation1 <- (elevation-mean(elevation))/sd(elevation)
elevation2 <- (elevation1)^2

# dates
dates <- as.matrix(data.frame(list(
  D1=tapply(au_2017$Day1, au_2017$aldStao_, function(x) x[1]),
  D2=tapply(au_2017$Day2, au_2017$aldStao_, function(x) x[1]),
  D3=tapply(au_2017$Day3, au_2017$aldStao_, function(x) x[1]),
  D4=tapply(au_2017$Day4, au_2017$aldStao_, function(x) x[1]),
  D5=tapply(au_2017$Day5, au_2017$aldStao_, function(x) x[1]),
  D6=tapply(au_2017$Day6, au_2017$aldStao_, function(x) x[1]),
  D7=tapply(au_2017$Day7, au_2017$aldStao_, function(x) x[1])
)))
dates1 <- (dates-mean(dates))/sd(dates)
dates2 <- (dates1)^2

# Compute starting values for z matrix (with augmented species)
zst <- array(NA, c(nsites,nperiod)) # xaug
# numeric
for(i in 1:nsites){
  for (j in 1:nperiod) {
    zst[i,j]<- as.numeric(xaug[i,j]>0)
  }
}

# Initial values
txaug <- xaug
txaug[is.na(xaug)] <- 0
txaug[!is.na(xaug)] <- NA # all is na?


########## Occupancy model ##########

# bundle data
bundle_data <- list(xaug=xaug, nsites=nsites, nperiod=nperiod, elevation1=elevation1, elevation2=elevation2, dates1=dates1, dates2=dates2)
# xaug: detection
# txaug: initial values for xaug (all NA)
# zst: xaug
# nsites: 69
# nperiod: 7
# dates1: standerized date
# dates2: (dates1)^2
# elevation1: standerized elevation
# elevation2: (elevation1)^2

# initial values
initial_values <- function() {list(z=zst, xaug=txaug)}

# parameters
parameters <- c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta2", "beta3", "beta4")

# MCMC settings
nc=3; ni=1000; nb=500; nt=2

# jags
cat(file="model_butterfly.txt", "model{

# Priors
alpha0 <- logit(p_intercept)
p_intercept ~ dunif(0,1)
beta0 <- logit(psi_intercept)
psi_intercept ~ dunif(0,1)
beta1 ~ dnorm(0,0.001)
beta2 ~ dnorm(0,0.001)
beta3 ~ dnorm(0,0.001)
beta4 ~ dnorm(0,0.001)
alpha1 ~ dnorm(0,0.001)
alpha2 ~ dnorm(0,0.001)

# occupancy model
for(i in 1:nsites){ 
  for(j in 1:nperiod){
    logit(p[i,j]) <- alpha0 + alpha1*dates1[i,j] + alpha2*dates2[i,j]
    logit(psi[i,j]) <- beta0 + beta1*dates1[i,j] + beta2*dates2[i,j] + beta3*elevation1[i] + beta4*elevation2[i]
    z[i,j] ~ dbern(psi[i,j])
    mu2[i,j] <- z[i,j]*p[i,j]
    xaug[i,j] ~ dbin(mu2[i,j], 2)
  }
}
}")

# run
library("jagsUI")
out <- jags(bundle_data, initial_values, parameters, "model_butterfly.txt", nc, nt, ni, nb, parallel=T)
traceplot(out)
print(out, dig=3)


########## CAR: Occupancy model with spatial autocorrelation ##########

# Neighbourhood information
library(spdep)
neighbour <- dnearneigh(cbind(au_2017$coordx, au_2017$coordy), 0, 10) # boundaries can be adjusted
table(card(neighbour)) # Number of neighbours
# Convert the neighbourhood
wan <- nb2WB(neighbour)

# bundle data
bundle_data <- list(xaug=xaug, nsites=nsites, nperiod=nperiod, elevation1=elevation1, elevation2=elevation2, dates1=dates1, dates2=dates2, weights=wan$weights, adj=wan$adj, num=wan$num)

# initial values
initial_values <- function() {list(z=zst, xaug=txaug)}

# parameters
parameters <- c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta2", "beta3", "beta4", "tau")

# MCMC settings
nc=3; ni=1000; nb=500; nt=2

# bugs
cat(file="model_car.txt", "model {

# Priors
alpha0 <- logit(p_intercept)
p_intercept ~ dunif(0.00001,0.99999)
beta0 <- logit(psi_intercept)
psi_intercept ~ dunif(0.00001,0.99999)
beta1 ~ dnorm(0,0.001)
beta2 ~ dnorm(0,0.001)
beta3 ~ dnorm(0,0.001)
beta4 ~ dnorm(0,0.001)
alpha1 ~ dnorm(0,0.001)
alpha2 ~ dnorm(0,0.001)
 
# car
vrho ~ dnorm(0,0.2)
tau <- 1/vrho
rho[1:nsites] ~ car.normal(adj[], weights[], num[], tau)
   
# occupancy model
for(i in 1:nsites){ 
  for(j in 1:nperiod){
    logit(q[i,j]) <- alpha0 + alpha1*dates1[i,j] + alpha2*dates2[i,j]
    p[i,j] <- max(0.00001,min(0.99999,q[i,j]))
    logit(qsi[i,j]) <- rho[i] + beta0 + beta1*dates1[i,j] + beta2*dates2[i,j] + beta3*elevation1[i] + beta4*elevation2[i]
    psi[i,j] <- max(0.00001,min(0.99999,qsi[i,j]))
    z[i,j] ~ dbern(psi[i,j])
    mu2[i,j] <- z[i,j]*p[i,j]
    xaug[i,j] ~ dbin(mu2[i,j], 2)
  }
}
}")

# run
library(R2WinBUGS)
out <- bugs(data=bundle_data, parameters.to.save=parameters, model.file="model_car.txt", inits=initial_values, n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, DIC=F, bugs.dir="C:/Program Files/WinBUGS14", codaPkg=F, save.history=T, debug=T)

