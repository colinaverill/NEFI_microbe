#Plot in-sample calibration fits for all 74 fungal taxonomic and functional groups.
rm(list=ls())
source('paths.r')

#Set output path motif for multiple plot outputs.----
#setting output path motif. Number of output paths depends on number of taxa we are plotting. 25 per panel.
#output.path_motif <- 'test'
output.path_motif <- ITS_prior_calibration_fits.path

#load data.----
fg <- readRDS(ted_ITS.prior_fg_JAGSfit_micronutrient)
pl <- readRDS(ted_ITS_prior_phylo.group_JAGSfits_micronutrient)


#Count number of taxa to be plotted, specify output paths.----
pl.count <- list()
n.panels <- 25 #needs to have an interger square root for this to work!
for(i in 1:length(pl)){pl.count[[i]] <- ncol(pl[[i]]$observed) - 1} #minus 1 because other group
pl.count <- do.call(sum, pl.count)
n.taxa <- pl.count + ncol(fg$all.preds$observed) - 1
n.plots <- round(n.taxa / n.panels)

#setup output paths.
output.paths <- paste0(output.path_motif,'_',1:n.plots,'.png')

#Generate a super matrices of predicted vs. observed.----
fg_pred <- fg$all.preds$predicted
fg_obs  <- fg$all.preds$observed
namey <- paste0('Function.Group_',colnames(fg_pred))
namey <- colnames(fg_pred)
colnames(fg_pred) <- namey
colnames(fg_obs ) <- namey
pl_pred <- list()
pl_obs  <- list()
for(i in 1:length(pl)){
  pl_pred[[i]] <- pl[[i]]$predicted
  pl_obs [[i]] <- pl[[i]]$observed
  namey <- paste0(names(pl)[i],'_',colnames(pl_pred[[i]]))
  namey <- Hmisc::capitalize(namey)
  colnames(pl_pred[[i]]) <- namey
  colnames(pl_obs [[i]]) <- namey
}
pred <- do.call(cbind,pl_pred)
 obs <- do.call(cbind,pl_obs )
pred <- cbind(fg_pred, pred)
 obs <- cbind(fg_obs , obs )
#Drop 'other' colums
pred <- pred[,-grep('other',colnames(pred))]
 obs <-  obs[,-grep('other',colnames( obs))]
#Convert to lists subsetted based on number of plots
pred.list <- list()
 obs.list <- list()
 for(i in 1:n.plots){
   range <- seq(1 + (i-1)*n.panels , i*n.panels, by = 1)
   if(i == n.plots){range <- range[range <= ncol(pred)]}
  pred.list[[i]] <- pred[,range]
   obs.list[[i]] <-  obs[,range]
}

#Start plot.----
for(k in 1:length(pred.list)){
  #grab plot subset.
  output.path <- output.paths[k]
  pred <- pred.list[[k]]
  obs  <-  obs.list[[k]]
  
  #Setup plot save path.
  
  png(filename=output.path,width=12,height=12,units='in',res=300)
  
  #global plot settings.
  par(mfrow=c(sqrt(n.panels), sqrt(n.panels)),
      mar = c(2,2,2,2),
      oma = c(5.6,5.8,1,1))
  bestfit.col <- 'purple'
  head.cex <- 1 #for plot titles and R2
  o.cex <- 2 #outer label cex.
  
  #loop over columns.
  for(i in 1:ncol(pred)){
    limy <- c(0, max(obs[,i]*1.1))
    plot(obs[,i] ~ pred[,i], pch = 16, cex = 0.8,ylab = NA, xlab = NA, ylim = limy)
    abline(0,1, lwd = 2)
    mod <- lm(obs[,i] ~ pred[,i])
    rsq <- round(summary(mod)$r.squared, 2)
    abline(mod, lty = 2, lwd = 2, col = bestfit.col)
    namey <- colnames(obs)[i]
    namey <- gsub('_',': ', namey)
    mtext(namey, side = 3, cex = head.cex, line = 0.2)
    lab <- bquote(R^2 == .(format(rsq, digits = 3)))
    mtext(lab, side = 3, line = -1.8, adj = 0.05, cex = head.cex)
  }
  mtext( 'Observed Relative Abundance', side = 2, cex = o.cex, line = 2.2, outer = T)
  mtext('Predicted Relative Abundance', side = 1, cex = o.cex, line = 3.0, outer = T)
  
  #End plot.
  dev.off()
}

#end script.----

