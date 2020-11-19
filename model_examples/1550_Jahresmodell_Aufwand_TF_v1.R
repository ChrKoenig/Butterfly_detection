#---------------------------------------------------------------------------------------------------
# SKRIPT ZUR ANALYSE DES TAGFALTERMODELS BESCHRIEBEN IN KERY ET AL: 2009
#
# Author: Ro
# Ref: 1550_Jahresmodell_Aufwand_TF_v1.R
# Datum: 01.12.2017
# 
# Bem.:	Das Model ist eine Erweiterung des Models in Kery et al 2009. Das Model in dieser 
#		Ausführung liefert wichtige Kenngrössen zu Tagfalteraufnahmen.
#
#---------------------------------------------------------------------------------------------------

# ERKLÄRUNG DER OUTPUT VARIABLEN
# Ntotal		Total der Arten in Bezugsraum
# N			Total Arten pro Site und Aufnahme
# Nsitetotal	Total der nachgewiesienen Arten pro Transekt
# alpha0out		Effekt der einzelnen Bearbeiter
# pout[]		Mittlere Entdeckungswahrscheinlichkeit einer Art

#---------------------------------------------------------------------------------------------------
# Voreinstellungen
#---------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))
setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_occupancy/")
library(R2jags)
library(readxl)
library(tidyverse)

filename<-"Output/1550 E3"

#---------------------------------------------------------------------------------------------------
# Modelleinstellungen
#---------------------------------------------------------------------------------------------------
# Jahre
jahre<-2017

# MCMC settings
nc = 3			# Number of chains
ni = 3000			# Number of iterations
nb = 1000		# Burnin length
nt = 2			# Thinning rate

# Augment Data with zeros
nzeroes = 30

# Set to zero psi-p correlation: add to data and eliminate from list of params to monitor
rho<-0			

# Anzahl Durchg?nge pro Saison
nperiod<-7

#---------------------------------------------------------------------------------------------------
# Some data manipulations
#---------------------------------------------------------------------------------------------------
# Daten einlesen
# Vorsicht: Nur Daten von von Tagfaltern mit vorhandener aIdErheb und keine Zusatzarten
alldat <- read.csv("Data/raw_data/1550_TF_Daten_blinded.csv",sep=",",as.is=T)
alldat<- alldat %>%
  filter(Aufnahmetyp == "Normalaufnahme_Z7") %>%
  filter(Flags.ZA == 0)
alldat$Aufnahmetyp <- NULL

# Nur relevante Daten ausw?hlen
alldat <- alldat[alldat$Jahr_ %in% jahre,]
alldat <- alldat[!is.na(alldat$aIdErhTagf),] #Kleinarten raus
alldat[alldat$BGR %in% c("Östliche Zentralalpen","Westliche Zentralalpen"),"BGR"] <- "Zentralalpen"
alldat$BGR[alldat$BGR=="Jura"]<-"JU"
alldat$BGR[alldat$BGR=="Mittelland"]<-"ML"
alldat$BGR[alldat$BGR=="Alpennordflanke"]<-"AN"
alldat$BGR[alldat$BGR=="Zentralalpen"]<-"ZA"
alldat$BGR[alldat$BGR=="Alpensüdflanke"]<-"AS"

# Fragezeichen in Datum loeschen, zu Zahlen umwandeln
alldat[alldat=="?" & !is.na(alldat)]<-NA
alldat[,paste("Day",1:7,sep="")]<-apply(alldat[,paste("Day",1:7,sep="")],2,as.numeric)

### Prepare matrix to save the estimates
a<-as.character()
namen<-c(unique(alldat$BGR),"HL","CH")
for (m in 1:length(namen)) {
  a[3*m-2]<-paste("n",namen[m])
  a[3*m-1]<-paste("p",namen[m]) 
  a[3*m] <- paste("ci",namen[m])
}
results<-matrix(NA,nrow=length(jahre)+1,ncol=length(a), 
                dimnames = list(c(jahre,"Mittel"),a))

#---------------------------------------------------------------------------------------------------
# BUGS model description
#---------------------------------------------------------------------------------------------------
sink("model.txt")
cat("
    model { 
    
    # Priors
    omega~dunif(0,1)
    
    mupsi~dnorm(0,.001)
    mup~dnorm(0,.001)
    
    tau.psi<-1/(sigma.psi*sigma.psi)
    tau.p<-1/(sigma.p*sigma.p)
    sigma.psi~dunif(0,10)
    sigma.p~dunif(0,10)
    
    rho~dunif(-1,1)
    var.eta<-tau.p/(1.-pow(rho,2))
    
    mupsi.season1~dnorm(0,.001)
    mupsi.season2~dnorm(0,.001)
    mupsi.elev1~dnorm(0,.001)
    mupsi.elev2~dnorm(0,.001)
    
    taupsi.season1 <-1/(sigmapsi.season1*sigmapsi.season1)
    taupsi.season2<-1/(sigmapsi.season2*sigmapsi.season2)
    taupsi.elev1 <-1/(sigmapsi.elev1*sigmapsi.elev1)
    taupsi.elev2<-1/(sigmapsi.elev2*sigmapsi.elev2)
    sigmapsi.season1~dunif(0,10)
    sigmapsi.season2~dunif(0,10)
    sigmapsi.elev1~dunif(0,10)
    sigmapsi.elev2~dunif(0,10)
    
    # mup.observer~dnorm(0,.001)
    mup.season1~dnorm(0,.001)
    mup.season2~dnorm(0,.001)
    taup.observer<-1/(sigmap.observer*sigmap.observer)
    taup.season1<-1/(sigmap.season1*sigmap.season1)
    taup.season2<-1/(sigmap.season2*sigmap.season2)
    sigmap.observer~dunif(0,10)
    sigmap.season1~dunif(0,10)
    sigmap.season2~dunif(0,10)
    
    # Priors for regression coefficients
    for(k in 1:(nspec+nzeroes)){  
      beta1[k]~dnorm(mupsi.season1,taupsi.season1) T(-12,12)
      beta2[k]~dnorm(mupsi.season2,taupsi.season2) T(-12,12)
      beta3[k]~dnorm(mupsi.elev1,taupsi.elev1) T(-12,12)
      beta4[k]~dnorm(mupsi.elev2,taupsi.elev2) T(-12,12)
      alpha1[k]~dnorm(mup.season1,taup.season1) T(-12,12)
      alpha2[k]~dnorm(mup.season2,taup.season2) T(-12,12)
    }
    for(o in 1:nobserver){
      alpha0[o] ~ dnorm(0,taup.observer) T(-12,12)  
    }
    
    # State model
    for(k in 1:(nspec+nzeroes)){
      w[k]~dbin(omega,1)
    
      lpsi[k]~dnorm(mupsi,tau.psi) T(-12,12)
      mu.lp[k]<- mup+(rho*sigma.p/sigma.psi)*(lpsi[k]-mupsi) 
      lp[k]  ~ dnorm(mu.lp[k],var.eta) T(-12,12)
    
      for(i in 1:nsite){  
        for(j in 1:nperiod){   
          # Regression of detection on date and observer
          logit(p[k,i,j])<- lp[k] + alpha0[observer[i]] + alpha1[k]*dates1[i,j]+alpha2[k]*dates2[i,j]
          # Regression of occurrence on date and elevation
          logit(psi[k,i,j])<- lpsi[k] + beta1[k]*dates1[i,j]+beta2[k]*dates2[i,j]+beta3[k]*elev1[i]+beta4[k]*elev2[i]
          mu.psi[k,i,j]<-psi[k,i,j]*w[k]
          z[k,i,j]~dbern(mu.psi[k,i,j])
        }
      }
    }
    
    # Binomial observation model for data 
    for(k in 1:(nspec+nzeroes)){  
      for (i in 1:nsite) { 
        for(j in 1:nperiod){  
          mu2[i,j,k]<- z[k,i,j]*p[k,i,j]  
          xaug[k,i,j]~dbin(mu2[i,j,k], 2)  
        }
      }
    }
    
    # Calculations to get number of species
    for(i in 1:nsite){
      for(j in 1:nperiod){
        N[i,j]<-sum(z[1:(nspec+nzeroes),i,j])  # Artenzahl pro Ort und Periode
    ",fill=TRUE)
sink()

#---------------------------------------------------------------------------------------------------
# Schleife um alle Jahre
#---------------------------------------------------------------------------------------------------
for (jahr in jahre){
  print(jahr)
  ### Daten des entprechneden Jahres vorbereiten
  # Nur Daten eines Jahres auswaehlen
  dat<-alldat[alldat$Jahr_==jahr,]
  sr<-table(dat$aldStao_) # Artenvielfalt inkl ZAs
  for (i in names(sr)) sr[i]->dat[dat$aldStao_==i,"sr"]

  ### Tabelle mit Angaben zu jedem Site vorbereiten
  siteinfo <- data.frame("SR"=as.numeric(sr),row.names=names(sr))
  tmp<-tapply(dat$BGR,dat$aldStao_,unique)
  siteinfo[names(tmp),"BGR"]<-tmp
  tmp<-tapply(dat$Hohe_Lagen,dat$aldStao_,unique)
  siteinfo[names(tmp),"Hohe_Lagen"]<-tmp
  tmp<-tapply(dat$Verdichtung,dat$aldStao_,unique)
  siteinfo[names(tmp),"Verdichtung"]<-tmp
  
  ### Angaben für Modell bereitstellen
  name_spec	<- sort(unique(dat$ArtCodeTagf))
  nspec  	<- length(name_spec)
  name_lat <- tapply(dat$NameTagfLat,dat$ArtCodeTagf,unique)
  name_site  <- sort(unique(dat$aldStao))
  nsite		<- length(name_site)
  
  # Array fuer Beobachtungsdaten
  xaug <- array(0, dim=c(nspec+nzeroes, nsite, nperiod), dimnames=list(c(name_spec,rep("zero",nzeroes)),name_site,1:7))
  for(k in 1:nspec){
    d <- dat[dat$ArtCodeTagf == name_spec[k],]
    for(i in 1:nsite){
      for(j in 1:nperiod){
        ausw <- d$aldStao_ == name_site[i]
        if(sum(ausw)==1) xaug[k,i,j] <- d[ausw, paste("PreAbs",j,sep="")]
        if(sum(ausw)>1) print(paste("Warnmeldung: Mehr als ein Nachweis von Art", name_spec[k], "in", name_site[i]))
      }
    }
  }
  
  xaug[xaug==2] <- 1
  xaug[xaug==3] <- 2
  xaug[xaug==-1] <- NA
  xaug[xaug==6] <- 2
  
  xaug[,tapply(dat$Day1, dat$aldStao_, function(x) is.na(x[1])),1] <- NA
  xaug[,tapply(dat$Day2, dat$aldStao_, function(x) is.na(x[1])),2] <- NA
  xaug[,tapply(dat$Day7, dat$aldStao_, function(x) is.na(x[1])),7] <- NA
  
  # Standardize elevation covariate and get squared values
  elev <- as.numeric(tapply(dat$Hoehe, dat$aldStao_, function(x) x[1]))
  mean.elev <- mean(elev)
  sd.elev <- sd(elev)
  elev1<-(elev-mean.elev)/sd.elev
  elev2<-elev1*elev1
  
  # Standardize date covariate and get squared values
  dates<-as.matrix(data.frame(list(
    D1=tapply(dat$Day1, dat$aldStao_, function(x) x[1]),
    D2=tapply(dat$Day2, dat$aldStao_, function(x) x[1]),
    D3=tapply(dat$Day3, dat$aldStao_, function(x) x[1]),
    D4=tapply(dat$Day4, dat$aldStao_, function(x) x[1]),
    D5=tapply(dat$Day5, dat$aldStao_, function(x) x[1]),
    D6=tapply(dat$Day6, dat$aldStao_, function(x) x[1]),
    D7=tapply(dat$Day7, dat$aldStao_, function(x) x[1])
  )))
  mean.date <- mean(dates[!is.na(dates)]) 
  sd.date <- sd(dates[!is.na(dates)])
  dates[is.na(dates)] <- mean.date
  dates1 <-(dates-mean.date)/sd.date
  dates2 <-dates1*dates1
  
  # Get observer identity
  dat$Bearbeiter<-as.factor(dat$Bearbeiter)
  observer.name <- as.character(tapply(dat$Bearbeiter, dat$aldStao_, function(x) as.character(x[1])))
  observer <- as.numeric(tapply(dat$Bearbeiter, dat$aldStao_, function(x) x[1]))
  nobserver<-length(unique(observer))
  
  # Compute starting values for z matrix (with augmented species)
  wst <- c(rep(1,nspec), rep(0,nzeroes))
  zst <- array(NA,c(nspec+nzeroes,nsite,nperiod))
  for(i in 1:nsite){
    for (j in 1:nperiod) {
      for(k in 1:(nspec+nzeroes)){
        zst[k,i,j]<- as.numeric(xaug[k,i,j]>0)
      }
    }
  }
  
  
  #---------------------------------------------------------------------------------------------------
  # Run model in WinBUGS
  #---------------------------------------------------------------------------------------------------
  # Bundle data
  jags.data <- list ("xaug","txaug","wst","zst","nsite","nperiod","nspec","nobserver","observer", "nzeroes", "dates1", "dates2", "elev1", "elev2","nc", "nb", "ni", "nt")
  
  ### Initial values
  txaug <- xaug
  txaug[is.na(xaug)] <- 0
  txaug[!is.na(xaug)] <- NA
  
  f.inits <- function(){
    list (w=wst,
          z = zst,
          xaug = txaug,
          mupsi=runif(1,-1,1),
          mup=runif(1,-1,1),
          mupsi.season1=runif(1,-1,1),
          mupsi.season2=runif(1,-1,1),
          mupsi.elev1=runif(1,-1,1),
          mupsi.elev2=runif(1,-1,1),
          mup.season1=runif(1,-1,1),
          mup.season2=runif(1,-1,1),
          sigma.psi=runif(1,0.1,3),
          sigma.p=runif(1,0.1,3),
          sigmapsi.season1=runif(1,0.1,3),
          sigmapsi.season2=runif(1,0.1,3),
          sigmapsi.elev1=runif(1,0.1,3),
          sigmapsi.elev2=runif(1,0.1,3),
          sigmap.season1=runif(1,0.1,3),
          sigmap.season2=runif(1,0.1,3),
          rho=0.5,
          beta1=rep(0, nspec+nzeroes),
          beta2=rep(0, nspec+nzeroes),
          beta3=rep(0, nspec+nzeroes),
          beta4=rep(0, nspec+nzeroes),
          alpha1=rep(0, nspec+nzeroes),
          alpha2=rep(0, nspec+nzeroes),
          lpsi=rep(0, nspec+nzeroes),
          lp=rep(0, nspec+nzeroes))
  }
  
  ### Parameters to estimate
  params <- c("Ntotal", "Nsitetotal", "lpsiout", "beta1out", "beta2out", "beta3out", "beta4out", "Ptotal")
  
  ### Laufen lassen
  start.time = Sys.time()
  out <- jags.parallel(jags.data,f.inits,params,"model.txt",nc,ni,nb,nt)
  end.time = Sys.time()
  elapsed.time = difftime(end.time, start.time, units='mins')
  cat(paste(paste('Posterior sample analysed in ', elapsed.time, sep=''), ' minutes\n', sep=''))
    
  #---------------------------------------------------------------------------------------------------
  # Make Result Files
  #---------------------------------------------------------------------------------------------------
  ### Allgemeine Resultate speichern
  pdf(paste(filename," Resultate Tagfaltermodel ",jahr,".pdf",sep=""), height=10, width=6, pointsize=10)
  par(mfrow=c(2,1))
  hist(out$BUGSoutput$sims.list$Ntotal, main="Total number of species", xlim=c(nspec-1, nspec+nzeroes),
       xlab="Anzahl Arten", breaks=40, freq=F)
  abline(v=nspec, col="red")
  hist(out$BUGSoutput$sims.list$Ptotal, main="Anteil nicht entdeckter Arten",
       xlab="P total", breaks=40, freq=F)
  dev.off()
    
  ### P-Werte je Site
  siteinfo$SR_SOM <- out$BUGSoutput$summary[paste("Nsitetotal[", 1:nsite, "]", sep=""),"50%"]
  siteinfo$p_werte <- siteinfo$SR/siteinfo$SR_SOM
  siteinfo$p_2 <- 1-((1-siteinfo$p_werte)*(1-siteinfo$p_werte))
  
  ### Resultate in Results
  # Je BGR
  tmp<-table(siteinfo$BGR)
  results[as.character(jahr),paste("n",names(tmp))]<-tmp
  tmp<-tapply(siteinfo$p_2,siteinfo$BGR,mean)
  results[as.character(jahr),paste("p",names(tmp))]<-round(tmp,2)
  tmp<-tapply(siteinfo$p_2,siteinfo$BGR,function(x) diff(t.test(x)$conf.int)/2)
  results[as.character(jahr),paste("ci",names(tmp))]<-round(tmp,3)
  
  # HL
  tmp<-siteinfo$p_2[siteinfo$Hohe_Lagen>90]
  results[as.character(jahr),"n HL"] <- length(tmp)
  results[as.character(jahr),"p HL"]<-round(mean(tmp),2)
  results[as.character(jahr),"ci HL"]<-round(diff(t.test(tmp)$conf.int)/2,3)
  
  # CH
  tmp<-siteinfo$p_2[siteinfo$Verdichtung=="nein"]
  results[as.character(jahr),"n CH"] <- length(tmp)
  results[as.character(jahr),"p CH"]<-round(mean(tmp),2)
  results[as.character(jahr),"ci CH"]<-round(diff(t.test(tmp)$conf.int)/2,3)
  
  ### Verteilung der Arten ?ber die Saison bei mittlerer Meereshoehe (ungefaehr 1200MueM)
  dat.pred2 <- seq(90,304, 2)
  mean.dat <- mean(dat.pred2)
  sd.dat <- sd(dat.pred2)
  dat1.pred <- (dat.pred2 - mean.dat) / sd.dat
  dat2.pred <- dat1.pred^2
  pdf(paste(filename," Vorkommen innerhalb Saison ",jahr,".pdf",sep=""), height=10, width=6, pointsize=8)
  par(mfrow=c(3,2))
  for(i in 1:nspec) {
    occu.pred <- plogis(out$BUGSoutput$mean$lpsiout[i] + out$BUGSoutput$mean$beta1out[i] *dat1.pred + out$BUGSoutput$mean$beta2out[i] *dat2.pred)
    plot(dat.pred2, occu.pred, xlab="Datum", ylab="", main=name_lat[as.character(name_spec[i])], axes=F, ty="l")
    axis(1,at=c(90,120,151,181,212,243,273,304), 
         labels=c("1. April","1. Mai","1. Juni","1. Juli","1. Aug","1. Sept","1. Okt", "1. Nov"))
  }
  dev.off()
    
  # Hoehenverbreitung der Arten in der Mitte der Saison
  elev.pred2 <- seq(200,3000, 10)
  elev1.pred <- (elev.pred2 - mean.elev) / sd.elev
  elev2.pred <- elev1.pred^2
  pdf(paste(filename," Vorkommen vs Meereshoehe ",jahr,".pdf",sep=""), height=10, width=6, pointsize=8)
  par(mfrow=c(3,2))
  for(i in 1:nspec) {
    occu.pred <- plogis(out$BUGSoutput$mean$lpsiout[i] + out$BUGSoutput$mean$beta3out[i] *elev1.pred + out$BUGSoutput$mean$beta4out[i] *elev2.pred)
    plot(elev.pred2, occu.pred, xlab="Hoehe ueber Meer (M)", ylab="", main=name_lat[as.character(name_spec[i])], axes=F, ty="l")
    axis(1)
  }
  dev.off()
  
  rm(out)
}

test<-as.data.frame(results, col.names = 1)
write_delim(test, paste(filename," Arterfassungsgrad.csv"), delim = ",")

#---------------------------------------------------------------------------------------------------
# Mean over all years
#---------------------------------------------------------------------------------------------------
library(xlsx)
res <- read_excel("1550 E3 Q-Kennzahlen e5_sk.xlsx", sheet = "Feldteam Arterfassung SiteOccup", skip = 11, n_max = 10)
#res <- read.xlsx("1360 E3 Q-Kennzahlen e1.xlsx", sheetName = "Feldteam Arterfassung SiteOccup", startRow = 12, endRow = 20)
names(res) <- c("Jahr", "JU", "JU.p", "JU.CI", "ML", "ML.p", "ML.CI", "AN", "AN.p", "AN.CI", "ZA", "ZA.p", "ZA.CI", 
                "AS", "AS.p", "AS.CI", "HL", "HL.p", "HL.CI", "CH", "CH.p", "CH.CI")

apply(res[, c("JU", "ML", "AN", "ZA", "AS", "HL", "CH")], 2, sum)
round(apply(res[, c("JU.p", "ML.p", "AN.p", "ZA.p", "AS.p", "HL.p", "CH.p")], 2, mean),2)
round(apply(res[, c("JU.CI", "ML.CI", "AN.CI", "ZA.CI", "AS.CI", "HL.CI", "CH.CI")], 2, function(x) mean(as.numeric(substr(as.character(x), start = 2, stop = 5)))),2)

yr <- res$Jahr
round(apply(res[, c("JU.p", "ML.p", "AN.p", "ZA.p", "AS.p", "HL.p", "CH.p")], 2, function(x) coef(lm(x ~ yr))[2]), 3)
round(apply(res[, c("JU.p", "ML.p", "AN.p", "ZA.p", "AS.p", "HL.p", "CH.p")], 2, function(x) {
  tt <- summary(lm(x ~ yr))$coefficients
  tt["yr", "Estimate"] - 2 * tt["yr", "Std. Error"]
  }), 3)
round(apply(res[, c("JU.p", "ML.p", "AN.p", "ZA.p", "AS.p", "HL.p", "CH.p")], 2, function(x) {
  tt <- summary(lm(x ~ yr))$coefficients
  tt["yr", "Estimate"] + 2 * tt["yr", "Std. Error"]
}), 3)



