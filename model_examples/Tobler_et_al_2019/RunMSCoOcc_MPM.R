########################################################################################################
## Script to run a multivariate probit multi-species co-occurrence model with imperfect detection.
## Supplemental material to: 
## Tobler et al. 2019 Joint species distribution models with species correlations 
## and imperfect detection, Ecology.
########################################################################################################

library(MCMCpack)
library(mclust)
library(corrplot)

df=1
ID=diag(n)

### Mutltivarate probit mulsti-species co-occurence model
modelFile='MSCoOcc_MPM.txt'

#Specify the data
occ.data = list(n=n, J=J, k=K, y=y,
                Xocc=Xocc,Xobs=Xobs,Vocc=ncol(Xocc),Vobs=ncol(Xobs),ID=ID, df=n+df)

#Specify the parameters to be monitored
occ.params = c('z','u.b','v.b','mu.u.b','tau.u.b','mu.v.b','tau.v.b','tau.cor',
               'p.fit','p.fitnew','d2','dnew2')

#Specify the initial values
occ.inits = function() {
  library(MCMCpack)
  tau.cor <- rwish(n+df, ID)
  sigma<-solve(tau.cor)
  u.b.est<-t(sapply(seq_len(ncol(y)),
                    function(x) {unname(coef(glm(((y>0)*1)[, x] ~ Xocc[, -1],
                                                 family=binomial(link=probit))))}))
  u.b<-u.b.est*sqrt(diag(sigma))
  v.b<-matrix(rnorm(Vobs*n),c(n,Vobs))
  mu.psi<-matrix(rbinom((n)*J, size=1, prob=1),
                 nrow=J, ncol=(n))
  mu.u.b=apply(u.b,2,mean)
  mu.v.b=apply(v.b,2,mean)
  sigma.u.b=apply(u.b,2,sd)
  sigma.v.b=apply(v.b,2,sd) 
  e=matrix(6,nrow=J, ncol=n)
  list(u.b=u.b,v.b=v.b,mu.u.b=mu.u.b,mu.v.b=mu.v.b,sigma.u.b=sigma.u.b,sigma.v.b=sigma.v.b,mu.psi=mu.psi,e=e,tau.cor=tau.cor)
}


#run the model in JAGS with R2jags
library(R2jags)
fit <- jags.parallel(occ.data, occ.inits, occ.params, modelFile,
                     n.chains=3, n.iter=350000, n.burnin=300000, n.thin=50)
fit<-fit$BUGSoutput


##correlation and occupancy coefficients
###########################
Tau<-fit$sims.list$tau.cor
Beta.raw<-fit$sims.list$u.b


# Calcuate the covariance matrix from Tau
Sigma2 <- apply(Tau, 1, solve) #inverse of matrix
dim(Sigma2) <- rev(dim(Tau))
Sigma2 <- aperm(Sigma2, c(3, 2, 1))

# Calculate the correlation matrix from Sigma2
Rho <- apply(Sigma2, 1, cov2cor)
dim(Rho) <- rev(dim(Sigma2))
Rho <- aperm(Rho, c(3, 2, 1))

# Calculate corrected Beta from Beta.raw and Sigma2
BetaOcc <- apply(Beta.raw, 3, function(x) x / t(sqrt(apply(Sigma2, 1, diag))))
dim(BetaOcc) <- dim(Beta.raw)

#calculate mean correlation matrix
cmest<-apply(Rho,c(2,3),mean)

#compare simulated and estimated correlation matrix
cm
cmest
cor(c(cmest[lower.tri(cmest)]),c(cm[lower.tri(cm)]))
summary(lm(c(cmest[lower.tri(cmest)])~c(cm[lower.tri(cm)])))
(rmse <- sqrt(mean((c(cmest[lower.tri(cmest)])-c(cm[lower.tri(cm)]))^2)))
plot(cm,cmest)
abline(0,1)


#plot correlation matrix
corplot<-cmest
colnames(corplot)<-uspecies
rownames(corplot)<-uspecies
corrplot(corplot,method="circle",order="hclust",type="lower")
corrplot(corplot,method="color",order="original",type="lower",outline="black",tl.pos="d",tl.col="black")
corrplot(corplot,method="number",order="original",type="upper",col="black",addgrid.col="black", add=T,tl.pos="d",tl.col="black", cl.pos="n")

#plot large number of species (export 900x800)
corrplot(corplot,method="color",order="FPC",type="lower",outline="black",tl.pos="ld",tl.col="black",tl.cex=0.6,tl.srt = 45)



#occupancy coefficients
bocc1<-apply(BetaOcc,c(2,3),mean) #mean
bocc2<-apply(BetaOcc,c(2,3),quantile,0.025) #lower CI
bocc3<-apply(BetaOcc,c(2,3),quantile,0.975) #upper CI
bocc4<-((bocc2<0 & bocc3<0) | (bocc2>0 & bocc3>0)) #significant
bocc5<-bocc1*bocc4
row.names(bocc1)<-uspecies
colnames(bocc1)<-colnames(Xocc)
row.names(bocc5)<-uspecies
colnames(bocc5)<-colnames(Xocc)
write.table(round(bocc5,3),file="clipboard",sep="\t")
write.table(round(bocc1,3),file="clipboard",sep="\t")

##Occupancy probabilities for all species and sites
#occupancy probability on normal scale for each site and species
mu.psi<-Xocc %*% t(bocc1)
#occupancy probability for each site
psi<-1-pnorm(0,mu.psi,1)



#correlation between original and estimated betas
bocc1
bocc

cor(c(bocc1),c(bocc))
summary(lm(c(bocc1)~c(bocc)))
(rmse <- sqrt(mean((c(bocc1)-c(bocc))^2)))


#detection coefficients
BetaObs<-fit$sims.list$v.b

bobs1<-apply(BetaObs,c(2,3),mean)
bobs2<-apply(BetaObs,c(2,3),quantile,0.025)
bobs3<-apply(BetaObs,c(2,3),quantile,0.975)
bobs4<-((bobs2<0 & bobs3<0) | (bobs2>0 & bobs3>0))
bobs5<-bobs1*bobs4
row.names(bobs1)<-uspecies
colnames(bobs1)<-colnames(Xobs)
row.names(bobs5)<-uspecies
colnames(bobs5)<-colnames(Xobs)
write.table(round(bobs5,3),file="clipboard",sep="\t")
write.table(round(bobs1,3),file="clipboard",sep="\t")

##Detection probability
pest<-plogis(Xobs %*% t(bobs1))


