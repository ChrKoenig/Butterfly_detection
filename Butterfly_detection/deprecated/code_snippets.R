# Fit unmarked model
unmarked_df = unmarkedFrameOccu(y = Y, siteCovs = data.frame(X_occ), obsCovs = list(day = data.frame(X_det[,1:7]), day_sq = data.frame(X_det[,8:14])))
unmarked_fit = occu(~ day + day_sq
                    ~ ., data = unmarked_df)
unmarked_smry = summary(unmarked_fit)

# observer random effect
taup.observer<-1/(sigmap.observer*sigmap.observer)
for(o in 1:nobserver){
  alpha0[o] ~ dnorm(0,taup.observer) T(-12,12)  
}
logit(p[k,i,j])<- lp[k] + alpha0[observer[i]] + alpha1[k]*dates1[i,j]+alpha2[k]*dates2[i,j]