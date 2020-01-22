#Plotting ECM, Russula and Ascomycota.
#Building so its easy to update which representative groups we use.
#General method to combine phylo and functional data/casts, and then loop over to plot.
#Core level plots work. Next is plot-level.
#plot level close- seems to be plotting too many points? summary stats seem right.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- 'test.png'
output.path <- NEON_cps_rep.groups_forecast_figure.path

#groups of interest----
namey <- c('Ectomycorrhizal','Russula','Ascomycota')
level <- c('function_group','genus','phylum')

#grab forecasts and observations of functional and phylogenetic groups.----
pl.cast <- readRDS(NEON_site_fcast_all_groups_1k_rare.path)
pl.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)
names(pl.cast )[names(pl.cast ) == 'fg'] <- 'function_group'
names(pl.truth)[names(pl.truth) == 'fg'] <- 'function_group'

#grab the data of interest based on 'namey', set above.
pl.core_mu <- list()
pl.plot_mu <- list()
pl.site_mu <- list()
pl.plot_lo95 <- list()
pl.plot_hi95 <- list()
pl.site_lo95 <- list()
pl.site_hi95 <- list()
for(i in 1:length(pl.truth)){
  #rel.core <- (pl.core[[i]]$abundances + 1) / pl.core[[i]]$seq_total
  col.namey <- colnames(pl.truth[[i]]$core.fit)[colnames(pl.truth[[i]]$core.fit) %in% namey]
  pl.core_mu  [[i]] <- data.frame(pl.truth[[i]]$core.fit     [,colnames(pl.truth[[i]]$core.fit     ) %in% namey])
  pl.plot_mu  [[i]] <- data.frame(pl.truth[[i]]$plot.fit$mean[,colnames(pl.truth[[i]]$plot.fit$mean) %in% namey])
  pl.site_mu  [[i]] <- data.frame(pl.truth[[i]]$site.fit$mean[,colnames(pl.truth[[i]]$site.fit$mean) %in% namey])
  pl.plot_lo95[[i]] <- data.frame(pl.truth[[i]]$plot.fit$lo95[,colnames(pl.truth[[i]]$plot.fit$mean) %in% namey])
  pl.plot_hi95[[i]] <- data.frame(pl.truth[[i]]$plot.fit$hi95[,colnames(pl.truth[[i]]$plot.fit$mean) %in% namey])
  pl.site_lo95[[i]] <- data.frame(pl.truth[[i]]$site.fit$lo95[,colnames(pl.truth[[i]]$site.fit$mean) %in% namey])
  pl.site_hi95[[i]] <- data.frame(pl.truth[[i]]$site.fit$hi95[,colnames(pl.truth[[i]]$site.fit$mean) %in% namey])
  colnames(pl.core_mu  [[i]]) <- col.namey
  colnames(pl.plot_mu  [[i]]) <- col.namey
  colnames(pl.site_mu  [[i]]) <- col.namey
  colnames(pl.plot_lo95[[i]]) <- col.namey
  colnames(pl.plot_hi95[[i]]) <- col.namey
  colnames(pl.site_lo95[[i]]) <- col.namey
  colnames(pl.site_hi95[[i]]) <- col.namey
}
pl.core_mu <- do.call(cbind,pl.core_mu)
pl.plot_mu <- do.call(cbind,pl.plot_mu)
pl.site_mu <- do.call(cbind,pl.site_mu)
pl.plot_lo95 <- do.call(cbind,pl.plot_lo95)
pl.plot_hi95 <- do.call(cbind,pl.plot_hi95)
pl.site_lo95 <- do.call(cbind,pl.site_lo95)
pl.site_hi95 <- do.call(cbind,pl.site_hi95)
#name problem.
rownames(pl.core_mu) <- gsub('-GEN','',rownames(pl.core_mu))

#DEFINE OUTLIER SITES for groups- DSNY-ECM in this case.----
out_sites <- c('DSNY')
out_spp   <- c('Ectomycorrhizal')

#png save line.----
png(filename=output.path,width=12,height=12,units='in',res=300)

#global plot settings.----
par(mfrow = c(3,3),
    mai = c(0.3,0.3,0.3,0.3),
    oma = c(4,6,3,1))
trans <- 0.3
limy <- c(0,1)
core.cex <- 0.7
plot.cex <- 1.0
site.cex <- 1.5
outer.cex <- 2
glob.pch <- 16
names <- namey
out.color <- 'gray'
bf_col <- 'magenta1' #best-fit regression line color.

#loop over functional groups.----
for(i in 1:length(names)){
  #core.level.----
  #organize data.
  fcast <- pl.cast[[level[i]]]$core.fit
  obs <- pl.core_mu
  obs <- obs[rownames(obs) %in% rownames(fcast$mean),]
  for(k in 1:length(fcast)){
    fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(obs),]
  }
  obs <- obs[,colnames(obs) %in% colnames(fcast$mean)]
  obs <- as.matrix(obs)
  colnames(obs) <- names[names %in% colnames(fcast$mean)]
  rownames(obs) <- rownames(pl.core_mu)
  pos <- which(colnames(fcast$mean) == names[i]) #position that matches the fungal type we are plotting.
  mu <- fcast$mean[,pos][order(fcast$mean[,pos])]
  ci_0.975 <- fcast$ci_0.975[,pos][order(match(names(fcast$ci_0.975[,pos]),names(mu)))]
  ci_0.025 <- fcast$ci_0.025[,pos][order(match(names(fcast$ci_0.025[,pos]),names(mu)))]
  pi_0.975 <- fcast$pi_0.975[,pos][order(match(names(fcast$pi_0.975[,pos]),names(mu)))]
  pi_0.025 <- fcast$pi_0.025[,pos][order(match(names(fcast$pi_0.025[,pos]),names(mu)))]
  fungi_name <- colnames(fcast$mean)[pos]
  obs.pos <- which(colnames(obs) == names[i])
  obs.mu   <- obs[,obs.pos][order(match(names(obs[,c(fungi_name)]),names(mu)))]
  
  #make DSNY sites light gray for Ectos.
  obs.cols <- rep('black',nrow(obs))
  if(names[i] %in% out_spp){
    obs.cols <- ifelse(substring(names(obs.mu),1,4) %in% out_sites,out.color,'black')
  }
  
  #get y-limit.
  obs_limit <- max(obs.mu)
  if(max(pi_0.975) > as.numeric(obs_limit)){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  if(names[i] == 'Arbuscular'){y_max = 0.1}
  
  #plot
  plot(obs.mu ~ mu, cex = core.cex, pch=glob.pch, ylim=c(0,y_max), ylab=NA, xlab = NA, col = obs.cols)
  mod_fit <- lm(obs.mu ~ mu)
  if(names[i] %in% out_spp){
    siteID <- substring(names(obs.mu),1,4)
    to_keep <- ifelse(siteID %in% out_sites, F, T)
    mod_fit <- lm(obs.mu[to_keep] ~mu[to_keep])
  }
  rsq <- round(summary(mod_fit)$r.squared,2)
  mtext(paste0('R2=',rsq), side = 3, line = -2.7, adj = 0.03)
  mtext(fungi_name, side = 2, line = 2.5, cex = 1.5)
  #add confidence interval.
  range <- mu
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(as.numeric(obs.mu) < pi_0.975 & as.numeric(obs.mu) > pi_0.025) / length(obs.mu),2) * 100
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.3, adj = 0.05)
  abline(0,1,lwd=2)
  abline(mod_fit, lwd =2, lty = 2, col = bf_col)
  
  
  #plot.level.----
  #organize data.
  fcast <- pl.cast[[level[i]]]$plot.fit
  obs <- list(pl.plot_mu,pl.plot_lo95,pl.plot_hi95)
  names(obs) <- c('mu','lo95','hi95')
  for(k in 1:length(obs)){
    rownames(obs[[k]]) <- gsub('.','_',rownames(obs[[k]]), fixed = T)
    obs[[k]] <- obs[[k]][rownames(obs[[k]]) %in% rownames(fcast$mean),]
  }
  for(k in 1:length(fcast)){
    fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(obs$mu),]
  }
  pos <- which(colnames(fcast$mean) == names[i]) #position that matches the fungal type we are plotting.
  mu <- fcast$mean[,pos][order(fcast$mean[,pos])]
  ci_0.975 <- fcast$ci_0.975[,pos][order(match(names(fcast$ci_0.975[,pos]),names(mu)))]
  ci_0.025 <- fcast$ci_0.025[,pos][order(match(names(fcast$ci_0.025[,pos]),names(mu)))]
  pi_0.975 <- fcast$pi_0.975[,pos][order(match(names(fcast$pi_0.975[,pos]),names(mu)))]
  pi_0.025 <- fcast$pi_0.025[,pos][order(match(names(fcast$pi_0.025[,pos]),names(mu)))]
  fungi_name <- colnames(fcast$mean)[pos]
  #make sure row order matches.
  for(k in 1:length(obs)){
    obs[[k]] <- obs[[k]][order(match(rownames(obs[[k]]), names(mu))),]
  }
  obs.mu   <- obs$mu  [,fungi_name]
  obs.lo95 <- obs$lo95[,fungi_name]
  obs.hi95 <- obs$hi95[,fungi_name]
  names(obs.mu)   <- names(mu)

  #Make out_sites sites gray for out_spp.
  obs.cols <- rep('black',length(obs.mu))
  if(names[i] %in% out_spp){
    obs.cols <- ifelse(substring(names(obs.mu),1,4) %in% out_sites,out.color,'black')
  }
  
  #get y-limit.
  obs_limit <- max(obs.hi95)
  if(max(pi_0.975) > obs_limit){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  if(names[i] == 'Arbuscular'){y_max = 0.1}
  
  #plot
  plot(obs.mu ~ mu, cex = plot.cex, pch=glob.pch, ylim=c(0,y_max), ylab=NA, xlab = NA, col = obs.cols)
  arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.00, angle=90, code=3, col = obs.cols)
  mod_fit <- lm(obs.mu ~ mu)
  if(names[i] %in% out_spp){
    siteID <- substring(names(obs.mu),1,4)
    to_keep <- ifelse(siteID %in% out_sites, F, T)
    mod_fit <- lm(obs.mu[to_keep] ~mu[to_keep])
  }
  rsq <- round(summary(mod_fit)$r.squared,2)
  mtext(paste0('R2=',rsq), side = 3, line = -2.7, adj = 0.03)
  #1:1 line
  abline(0,1, lwd = 2)
  abline(mod_fit, lwd =2, lty = 2, col = bf_col)
  
  #add confidence interval.
  range <- mu
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.3, adj = 0.05)
  
  
  #site.level.----
  #organize data.
  fcast <- pl.cast[[level[i]]]$site.fit
  obs <- list(pl.site_mu,pl.site_lo95,pl.site_hi95)
  names(obs) <- c('mu','lo95','hi95')
  for(k in 1:length(obs)){
    rownames(obs[[k]]) <- gsub('.','_',rownames(obs[[k]]), fixed = T)
    obs[[k]] <- obs[[k]][rownames(obs[[k]]) %in% rownames(fcast$mean),]
  }
  for(k in 1:length(fcast)){
    fcast[[k]] <- fcast[[k]][rownames(fcast[[k]]) %in% rownames(obs$mu),]
  }
  pos <- which(colnames(fcast$mean) == names[i]) #position that matches the fungal type we are plotting.
  mu <- fcast$mean[,pos][order(fcast$mean[,pos])]
  ci_0.975 <- fcast$ci_0.975[,pos][order(match(names(fcast$ci_0.975[,pos]),names(mu)))]
  ci_0.025 <- fcast$ci_0.025[,pos][order(match(names(fcast$ci_0.025[,pos]),names(mu)))]
  pi_0.975 <- fcast$pi_0.975[,pos][order(match(names(fcast$pi_0.975[,pos]),names(mu)))]
  pi_0.025 <- fcast$pi_0.025[,pos][order(match(names(fcast$pi_0.025[,pos]),names(mu)))]
  fungi_name <- colnames(fcast$mean)[pos]
  #make sure row order matches.
  for(k in 1:length(obs)){
    obs[[k]] <- obs[[k]][order(match(rownames(obs[[k]]), names(mu))),]
  }
  obs.mu   <- obs$mu  [,fungi_name]
  obs.lo95 <- obs$lo95[,fungi_name]
  obs.hi95 <- obs$hi95[,fungi_name]
  names(obs.mu)   <- names(mu)
  
  #Make out_sites sites gray for out_spp.
  obs.cols <- rep('black',length(obs.mu))
  if(names[i] %in% out_spp){
    obs.cols <- ifelse(substring(names(obs.mu),1,4) %in% out_sites,out.color,'black')
  }
  
  #get y-limit.
  obs_limit <- max(obs.hi95)
  if(max(pi_0.975) > obs_limit){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  if(names[i] == 'Arbuscular'){y_max = 0.1}
  
  #plot
  plot(obs.mu ~ mu, cex = site.cex, pch=glob.pch, ylim=c(0,y_max), ylab=NA, xlab = NA, col = obs.cols)
  arrows(c(mu), obs.lo95, c(mu), obs.hi95, length=0.0, angle=90, code=3, col = obs.cols)
  mod_fit <- lm(obs.mu ~ mu)
  if(names[i] %in% out_spp){
    siteID <- substring(names(obs.mu),1,4)
    to_keep <- ifelse(siteID %in% out_sites, F, T)
    mod_fit <- lm(obs.mu[to_keep] ~mu[to_keep])
  }
  rsq <- round(summary(mod_fit)$r.squared,2)
  mtext(paste0('R2=',rsq), side = 3, line = -2.7, adj = 0.03)
  #1:1 line
  abline(0,1, lwd = 2)
  abline(mod_fit, lwd =2, lty = 2, col = bf_col)
  
  #add confidence interval.
  range <- mu
  polygon(c(range, rev(range)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(range, rev(range)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.3, adj = 0.05)
}

#outer labels.----
mtext('core-level', side = 3, line = -1.8, adj = 0.11, cex = outer.cex, outer = T)
mtext('plot-level', side = 3, line = -1.8, adj = 0.50, cex = outer.cex, outer = T)
mtext('site-level', side = 3, line = -1.8, adj = 0.88, cex = outer.cex, outer = T)
mtext( 'observed relative abundance', side = 2, line = 3, cex = outer.cex, outer = T)
mtext('predicted relative abundance', side = 1, line = 2, cex = outer.cex, outer = T)

#end plot.----
dev.off()

