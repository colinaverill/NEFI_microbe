#check fits, convergence plots and prsf scores for calibration models

rm(list=ls())
source('paths.r')
source('paths_fall2019.r')
library(runjags)
library(betareg)
library(coda)
library(ddpcr) # for quiet() function
source('NEFI_functions/crib_fun.r')

allfits <- readRDS(prior_delgado_ddirch_16S.path)

# has everything converged? print any high prsf scores.
for (i in 1:length(allfits)) {
  print(names(allfits)[[i]])
  fit <- allfits[[i]]
  summary(fit$jags_model)
  quiet(s <- summary(fit$jags_model))
  print(s[which(s[,11] > 1.1),])
}
fit$species_parameter_output

# save plots.
#pdf("/fs/data3/caverill/NEFI_data/16S/pecan_gen/figures/prior_fit_ddirch_all_groups.pdf")

#check the plots.
#allfits <- output.list
r2_all <- list()
par(mfrow = c(3,3))

for (p in 1:length(allfits)) {
fit <- allfits[[p]]
r2_lev <- list()
for(i in 1:ncol(fit$predicted)) {
  if (colnames(fit$predicted)[i]=="other") next()
  plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i],
       pch = 16, xlab="predicted abundance", ylab="observed abundance")
  Axis(x="predicted", side=2)
  abline(0,1,lwd = 2)
  abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
  
  
  
  
  #mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
  #rsq <-round(summary(mod)$pseudo.r.squared, 3)
  
  mod <- lm(fit$observed[,i] ~ fit$predicted[,i])
  rsq <- round(summary(mod)$r.squared, 2)
  #r2_out[[p]] <- rsq
  
  
  mtext(colnames(fit$predicted)[i], side = 3)
  mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
  r2_lev[[i]] <- rsq 
}
r2_all[[p]] <- r2_lev
}
dev.off()
names(r2_all) <- names(allfits)


# view summary by rank
lev.mu <- lapply(r2_all[1:5], unlist)
lev.mu <- lapply(lev.mu, mean)
tax.mu <- mean(unlist(r2_all[1:5]))
fg.mu <- mean(unlist(r2_all[6:18]))
n.mu <- mean(unlist(r2_all[7:13]))
c.mu <- mean(unlist(r2_all[c(6,14:16)]))
co.mu <- mean(unlist(r2_all[c(17:18)]))
lev.mu <- c(fg.mu, lev.mu)

names(lev.mu)[[1]] <- "functional"




lev.mu <- unlist(lev.mu)
xx <- seq(1,length(lev.mu))
ylim <- c(0,.7)
plot(lev.mu ~ xx, cex = 2.5, ylim = ylim, pch = 16, ylab = NA, xlab = NA, bty='n', xaxt = 'n', col="grey")
points(lev.mu ~ xx, cex = 2.5, pch = 16)
lines(xx, lev.mu, lty = 2, col="grey")
mtext(expression(paste("Calibration site-level R"^"2")), side = 2, line = 2.2, cex = o.cex)
axis(1, labels = F)
text(x=xx, y = -0.05, labels= names(lev.mu), srt=45, adj=1, xpd=TRUE, cex = 1.6)
text(x=xx, y = .17, labels= names(lev.mu), srt=45, adj=1, xpd=TRUE, cex = 1.8)



