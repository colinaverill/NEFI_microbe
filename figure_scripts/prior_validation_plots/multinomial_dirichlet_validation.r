#Plotting multi-nomial-dirichlet fits to the training data set - are our models any good "in-sample"?
#clear environment, source paths, pacakges and functions.
rm(list=ls())
source('paths.r')
library(runjags)
source('NEFI_functions/crib_fun.r')

#load data.
fit <- readRDS(ted_ITS.prior_dmulti.ddirch_fg_JAGSfit)
fit <- fit$all.preds


#CHECK JAGS diagnostic plots. Did chains converge?
#plot(fit$jags_model)
#check prsf scores. Pathogens could run a little longer, probably because low abundance.
fit$species_parameter_output

#check plots.
par(mfrow = c(2,2))
for(i in 2:ncol(fit$predicted)){
  plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i], pch = 16)
  abline(0,1,lwd = 2)
  abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
  mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  rsq <-round(summary(mod)$pseudo.r.squared, 2)
  mtext(colnames(fit$predicted)[i], side = 3)
  mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
  
}


#Check cosmo genera.
fit <- readRDS(ted_ITS.prior_dmulti.ddirch_cosmo_JAGSfit)
fit <- fit$all.preds

#check convergence plots and prsf scores. Did not all converge.
fit$species_parameter_output
#plot(fit$jags_model)

#check the plots.
par(mfrow = c(3,3))
for(i in 1:ncol(fit$predicted)){
  plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i], pch = 16)
  abline(0,1,lwd = 2)
  abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
  mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  rsq <-round(summary(mod)$pseudo.r.squared, 2)
  mtext(colnames(fit$predicted)[i], side = 3)
  mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
  
}
