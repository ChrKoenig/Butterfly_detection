########################################################################################################
## Script to run a latent variable multi-species co-occurrence model with imperfect detection.
## Supplemental material to: 
## Tobler et al. 2019 Joint species distribution models with species correlations 
## and imperfect detection, Ecology.
########################################################################################################

library(MCMCpack)
library(mclust)
library(corrplot)

### Latent variable multi-species co-occurence model
modelFile='MSCoOcc_LVM.txt'

#Number of latent variables to use
nlv=5

#Specify the data
occ.data = list(n=n, J=J, k=K, y=y,
                Xocc=Xocc,Xobs=Xobs,Vocc=ncol(Xocc),Vobs=ncol(Xobs),nlv=nlv)

#Specify the parameters to be monitored
occ.params = c('z','u.b','v.b','mu.u.b','tau.u.b','mu.v.b','tau.v.b','LV','lv.coef')

#Specify the initial values
occ.inits = function() {
  lv.coef<-matrix(1,n,nlv)
  lv.coef[1:nlv,1:nlv]<-0
  for(l in 1:nlv-1){
    lv.coef[l,(l+1):nlv]<-NA
  }
  u.b<-t(sapply(seq_len(ncol(y)),
                function(x) {unname(coef(glm(((y>0)*1)[, x] ~ Xocc[, -1],
                                             family=binomial(link=probit))))}))
  list(
    #u.b=matrix(rnorm(Vocc*n),c(n,Vocc)),
    u.b=u.b,
    v.b=matrix(rnorm(Vobs*n),c(n,Vobs)),
    u=(y>0)-runif(1,0.1,0.8),
    #mu.a = matrix(rbinom((n)*J, size=1, prob=1),
    #              nrow=J, ncol=(n)),
    #lv.coef=matrix(runif(nlv*n,-1,1),n,nlv)*lv.coef,
    LV=matrix(rnorm(nlv*J),J,nlv)
    #z = matrix(rbinom((n)*J, size=1, prob=1),
    #           nrow=J, ncol=(n))
  )
}


occ.data = list(n=n, J=J, k=K, y=y,
                Xocc=Xocc,Xobs=as.matrix(Xobs),Vocc=ncol(Xocc),Vobs=ncol(Xobs),nlv=nlv)

#run the model in JAGS with R2jags
library(R2jags)
fit <- jags.parallel(occ.data, occ.inits, occ.params, modelFile,
            n.chains=3, n.iter=20000, n.burnin=10000, n.thin=10)
fit<-fit$BUGSoutput



#calculate the correlation matrix from the latent variables
lv.coef<-fit$sims.list$lv.coef
cmall<-array(NA,dim=c(dim(lv.coef)[1],n.species,n.species))
eps.res<-apply(lv.coef,c(1,2),function(x)1-sum(x^2))
for(i in 1:dim(lv.coef)[1]){ #for each mcmc sample
  cmall[i,,]<-cov2cor(tcrossprod(lv.coef[i,,]) + diag(apply(eps.res,2,mean)))
}
cmest<-apply(cmall,c(2,3),mean)

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

#plot large number of species
corrplot(corplot,method="color",order="FPC",type="lower",outline="black",tl.pos="ld",tl.col="black",tl.cex=0.6,tl.srt = 45)


#occupancy coefficients
BetaOcc<-fit$sims.list$u.b

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



