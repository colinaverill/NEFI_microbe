#Show which groups fit very poorly at the site levelbut fit well in sample.
rm(list=ls())
source('paths.r')

#load data.----
#grab prior fits.
prior <- readRDS(ted_ITS_prior_all.groups_JAGSfits.path) #all phylo and functional groups.
#grab forecasts and observations of functional and phylogenetic groups.----
fcast <- readRDS(NEON_site_fcast_all_groups_1k_rare.path)
neon.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)

#get prior rsq values.----
prior.out <- list()
for(i in 1:length(prior)){
  lev <- prior[[i]]
  lev.out <- list()
  for(k in 1:ncol(lev$observed)){
    fit <- lm(lev$observed[,k] ~ lev$predicted[,k])
    rsq <- summary(fit)$r.squared
    lev.out[[k]] <- rsq
  }
  tax <- colnames(lev$predicted)
  lev.out <- unlist(lev.out)
  return <- data.frame(tax,lev.out)
  colnames(return) <- c('tax','rsq')
  return <- return[return$tax != 'other',]
  prior.out[[i]] <- return
}
prior.out <- do.call(rbind, prior.out)

#Get out of sample rsq values.----
fcast.out <- list()
for(i in 1:length(fcast)){
  pred <- fcast[[i]]$site.fit$mean
   obs <- neon.truth[[i]]$site.fit$mean
   lev.out <- list()
   for(k in 1:ncol(pred)){
     fit <- lm(obs[,k] ~ pred[,k])
     rsq <- summary(fit)$r.squared
     lev.out[[k]] <- rsq
   }
   tax <- colnames(pred)
   lev.out <- unlist(lev.out)
   return <- data.frame(tax,lev.out)
   colnames(return) <- c('tax','rsq')
   return <- return[return$tax != 'other',]
   fcast.out[[i]] <- return
}
fcast.out <- do.call(rbind, fcast.out)
out <- merge(prior.out, fcast.out, by = 'tax')

#Find taxa with highest in sample and lowest out of sample rsq.----
out <- out[order(out$rsq.x, decreasing = T),]

#great candidates are:
#1. Inocybe
#2. Hypocreales.
#3. Helotiales.

#grab predicted vs. observed for groups of interest----
namey <- c('Inocybe','Hypocreales','Helotiales')
level <- c('genus','order','order')
namey <- c('Inocybe','Thelephoraceae','Umbelopsis')
level <- c('genus','family','genus')

prior.plot <- list()
fcast.plot <- list()
for(i in 1:length(namey)){
  #prior pred vs. obs.
  pred <- prior[[which(names(prior) == level[i])]]$predicted
  pred <- pred[,colnames(pred) == namey[i]]
   obs <- prior[[which(names(prior) == level[i])]]$observed
   obs <-  obs[,colnames( obs) == namey[i]]
   prior.plot[[i]] <- data.frame(obs, pred)
  #fcst pred vs. obs.
  pred <- fcast[[which(names(fcast) == level[i])]]$site.fit$mean
   obs <- neon.truth[[which(names(fcast) == level[i])]]$site.fit$mean
   fcast.plot[[i]] <- data.frame(obs, pred)
}
names(prior.plot) <- namey
names(fcast.plot) <- namey

#Plot fits.----
png('test.png', height = 8, width = 6, units = 'in', res = 300)
par(mfrow = c(3,2),
    mar = c(1,1,1,1),
    oma = c(4,4,3,1))
trans = 0.4
for(i in 1:length(namey)){
  #grab observations, set x and y limits
   obs1 <- prior.plot[[i]][,1]
  pred1 <- prior.plot[[i]][,2]
   obs2 <- fcast.plot[[i]][,1]
  pred2 <- fcast.plot[[i]][,2]
   lab <- namey[i]
   limx <- c(min(pred1,pred2), max(pred1,pred2))
   limy <- c(min(obs1 ,obs2 ), max(obs1 , obs2))
   #prior plot panel.
   mod <- lm((obs1) ~ (pred1))
   rsq <- round(summary(mod)$r.squared, 2)
   plot((obs1) ~ (pred1), pch = 16, bty = 'l', ylab = NA, xlab = NA, xlim = limx, ylim = limy)
   abline(0,1,lwd = 2)
   abline(mod, lwd = 2, lty = 2, col= adjustcolor('purple', trans))
   mtext(lab, side = 3, adj = 0.05, line = -2.0)
   mtext(bquote(R^2 == .(rsq)), side = 3, adj = 0.05, line = -3.75)
  
  #fcast plot panel.
   mod <- lm((obs2) ~ (pred2))
   rsq <- round(summary(mod)$r.squared, 2)
   plot((obs2) ~ (pred2), pch = 16, bty = 'l', ylab = NA, xlab = NA, xlim = limx, ylim = limy)
   abline(0,1,lwd = 2)
   abline(mod, lwd = 2, lty = 2, col= adjustcolor('purple', trans))
   mtext(bquote(R^2 == .(rsq)), side = 3, adj = 0.05, line = -8)
 
}
#outer labels.
mtext('predicted', side = 1, cex = 1.5, outer = T, line = 2.5)
mtext('observed' , side = 2, cex = 1.5, outer = T, line = 2.0)
mtext('in-sample'    , side = 3, cex = 1.5, outer = T, line = -0.5, adj = 0.20)
mtext('out-of-sample', side = 3, cex = 1.5, outer = T, line = -0.5, adj = 0.875)

#end plot.
dev.off()