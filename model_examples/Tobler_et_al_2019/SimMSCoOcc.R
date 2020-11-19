########################################################################################################
## Script to simulate correlated presence/absence data with imperfect detection (occupancy design)
## Supplemental material to: 
## Tobler et al. 2019 Joint species distribution models with species correlations 
## and imperfect detection, Ecology.
########################################################################################################

library(MCMCpack)
library(mclust)
library(Matrix)


#define parameters
n.env.var <- 2 #number of random environmental variables for occupancy
n.sites <- 500 #number of sites
n.species <- 10 #number of species
p.min <- 0.1 #minimum detection probability
p.max <- 0.7 #maximum detection probability
T=3  #number of survyes
cor.lv=5 #number of latent variables used to simulate correlation matrix (method #2)

#random environmental variables from Normal(0,1)
Xocc <-  matrix(runif(n.sites * n.env.var,-1,1), ncol=n.env.var)
Xocc <- cbind(1, Xocc)   #add intercept
head(Xocc)

#observation matrix, intercept only
Xobs <- matrix(rep(1,n.sites), ncol=1) #we could add observatio covaraites here as above
#Xobs <- cbind(1, Xobs) #add intercept, use if there are observation covariates
head(Xobs)


### Two ways to simulate a correlation matrix

###1: unstructured, random correlation matrix
R <- matrix(runif(n.species^2,-1,1),nrow=n.species,ncol=n.species)
R <- (R * lower.tri(R)) + t(R * lower.tri(R))
diag(R) <- 1
cm<-as.matrix(nearPD(R, posd.tol=1.e-04,corr=T)$mat)


###2: use latent variables to simulate structured correlation matrix 
trueloadings <- matrix(runif(n.species*cor.lv,-1,1), ncol = cor.lv) 
trueloadings[upper.tri(trueloadings, diag = FALSE)] <- 0 # upper diagonal 0
diag(trueloadings) <- abs(diag(trueloadings)) 
cm<-cov2cor(tcrossprod(trueloadings))


#occupancy coefficients for each species
bocc <- matrix(rnorm(n.species * (n.env.var + 1),0,0.8), ncol=n.env.var + 1)
bocc  #regression coefficients for occupancy


#occupancy probability on normal scale for each site and species
mu.psi<-Xocc %*% t(bocc)
#occupancy probability for each species and site without correlation
psi<-1-pnorm(0,mu.psi,1)
#correlated true presence/absence
mu.z<-t(apply(mu.psi,1,function(x)mvrnorm(n=1, mu=x, Sigma=cm)))
z<-ifelse(mu.z > 0, 1, 0) 



#simulate detection probabilities for all species
p<-runif(n.species,p.min, p.max)
#p<-plogis(Xobs %*% t(bobs)) #use this when adding detection covaraiates

#simulate detection history
y<-apply(t(t(z)*p),1:2,function(x)rbinom(1,T,x))


# number of covariates
Vocc <- ncol(Xocc)
Vobs <- ncol(Xobs)

# number of occasions for each site
K<-rep(T,n.sites)

#for compatibility with run code
n=n.species
J=n.sites
uspecies<-paste("sp",1:n)







