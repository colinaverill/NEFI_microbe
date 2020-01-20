# check fits of 16S Bahram taxonomic and functional group models.

rm(list=ls())
source('paths.r')
library(runjags)
library(betareg)
library(coda)
library(ddpcr) # for quiet()
source('NEFI_functions/crib_fun.r')

#Load JAGS model.
allfits <- readRDS(bahram_16S_prior_dmulti.ddirch_all.group_JAGSfits)
allfits <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/JAGS_output/bahram_16S_prior_dmulti.ddirch_fg_JAGSfits")

# has everything converged? print any high prsf scores for all groups.
for (i in 1:length(allfits)) {
  print(names(allfits)[i])
  fit <- allfits[[i]]
  quiet(s <- summary(fit$jags_model))
  print(s[which(s[,11] > 1.1),])
}
fit$species_parameter_output

# save plots.
pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/prior_fit_dmulti_ddirch_all_groups.pdf")

#check the plots.
par(mfrow = c(2,2))
for (p in 1:length(allfits)) {
  fit <- allfits[[p]]

  for(i in 1:ncol(fit$predicted)) {
    if (colnames(fit$predicted)[i]=="other") next
    plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i], 
         pch = 16, xlab="predicted abundance", ylab="observed abundance")
    Axis(x="predicted", side=2)
    abline(0,1,lwd = 2)
    abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
    mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
    rsq <-round(summary(mod)$pseudo.r.squared, 3)
    mtext(colnames(fit$predicted)[i], side = 3)
    mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
  }
}
dev.off()
