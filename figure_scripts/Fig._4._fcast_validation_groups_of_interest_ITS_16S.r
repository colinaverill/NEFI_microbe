#Plotting 2 fungi and a bacteria - from Dirichlet fit.
#Built so its easy to update which representative groups we use.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
source('NEFI_functions/rsq_1.1.r')

#set output path.----
output.path <- 'figures/Fig._4._representative_groups.pdf'
#output.path <- 'figures/test.pdf'

#groups of interest----
namey <- c('Russula','Ascomycota','oligotroph')
level <- c('genus','phylum','oligotroph')

#grab forecasts and observations of functional and phylogenetic groups.----
fcast.ITS <- readRDS(NEON_site_fcast_all_groups_1k_rare.path)
truth.ITS <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)
fcast.16S <- readRDS(NEON_cps_fcast_ddirch_16S.path)
truth.16S <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)
names(fcast.ITS)[names(fcast.ITS) == 'fg'] <- 'function_group'
names(truth.ITS)[names(truth.ITS) == 'fg'] <- 'function_group'
#map of 16S sample codes.
map.16S <- readRDS(core_obs.path)


#Grab names of ITS and 16S taxa.----
names.ITS <- list()
names.16S <- list()
#ITS
for(i in 1:length(fcast.ITS)){
  nam <- colnames(fcast.ITS[[i]]$core.fit$mean)
  nam <- nam[!nam %in% c('other')]
  names.ITS[[i]] <- nam
}
names.ITS <- unlist(names.ITS)
#16S
for(i in 1:length(fcast.16S)){
  nam <- colnames(fcast.16S[[i]]$core.fit$mean)
  nam <- nam[!nam %in% c('other')]
  names.16S[[i]] <- nam
}
names.16S <- unlist(names.16S)

#grab the observed data of interest based on 'namey', set above.----
core_mu   <- list()
plot_mu   <- list()
site_mu   <- list()
plot_lo95 <- list()
plot_hi95 <- list()
site_lo95 <- list()
site_hi95 <- list()
for(i in 1:length(namey)){
  #check if its a bacteria or fungus, and level.
  if(namey[i] %in% names.ITS){
    truth <- truth.ITS
    fcast <- fcast.ITS
  }
  if(namey[i] %in% names.16S){
    truth <- truth.16S
    fcast <- fcast.16S
  }
  core_mu  [[i]] <- data.frame(truth[[level[i]]]$core.fit     [,colnames(truth[[level[i]]]$core.fit     ) %in% namey[i]])
  plot_mu  [[i]] <- data.frame(truth[[level[i]]]$plot.fit$mean[,colnames(truth[[level[i]]]$plot.fit$mean) %in% namey[i]])
  site_mu  [[i]] <- data.frame(truth[[level[i]]]$site.fit$mean[,colnames(truth[[level[i]]]$site.fit$mean) %in% namey[i]])
  plot_lo95[[i]] <- data.frame(truth[[level[i]]]$plot.fit$lo95[,colnames(truth[[level[i]]]$plot.fit$lo95) %in% namey[i]])
  plot_hi95[[i]] <- data.frame(truth[[level[i]]]$plot.fit$hi95[,colnames(truth[[level[i]]]$plot.fit$hi95) %in% namey[i]])
  site_lo95[[i]] <- data.frame(truth[[level[i]]]$site.fit$lo95[,colnames(truth[[level[i]]]$site.fit$lo95) %in% namey[i]])
  site_hi95[[i]] <- data.frame(truth[[level[i]]]$site.fit$hi95[,colnames(truth[[level[i]]]$site.fit$hi95) %in% namey[i]])
  colnames(core_mu  [[i]]) <- namey[i]
  colnames(plot_mu  [[i]]) <- namey[i]
  colnames(site_mu  [[i]]) <- namey[i]
  colnames(plot_lo95[[i]]) <- namey[i]
  colnames(plot_hi95[[i]]) <- namey[i]
  colnames(site_lo95[[i]]) <- namey[i]
  colnames(site_hi95[[i]]) <- namey[i]
}
#name problem.
for(i in 1:length(core_mu)){
  rownames(core_mu[[i]]) <- gsub('-GEN','',rownames(core_mu[[i]]))
}

#png save line.----
#png(filename=output.path,width=12,height=12,units='in',res=300)
pdf(file = output.path, width = 7.087, height = 8.5)

#global plot settings.----
par(mfrow = c(3,3),
    mai = c(0.3,0.3,0.3,0.3),
    oma = c(4,6,1,0.2))
trans <- 0.3
limy <- c(0,1)
core.cex <- 0.7
plot.cex <- 1.0
site.cex <- 1.5
outer.cex <- 1.25
rsq.lab.cex <- 0.65
glob.pch <- 16
names <- namey
out.color <- 'gray'
bf_col <- 'magenta1' #best-fit regression line color.

#loop over functional groups.----
for(i in 1:length(namey)){
  #define if working with 16S or ITS fcast.----
  if(namey[i] %in% names.ITS){fcast <- fcast.ITS}
  if(namey[i] %in% names.16S){fcast <- fcast.16S}
  
  #core.level.----
  #organize data.
  lev.cast <- fcast[[level[i]]]$core.fit
  obs <- core_mu[[i]]
  #If this is 16S, names are messed up.
  if(namey[i] %in% names.16S){
    map.name <- map.16S[,c('geneticSampleID','deprecatedVialID')]
    map.name <- map.name[map.name$deprecatedVialID %in% rownames(obs),]
    obs <- obs[rownames(obs) %in% map.name$deprecatedVialID,,drop = F]
    map.name <- map.name[order(match(map.name$deprecatedVialID, rownames(obs))),]
    map.name$geneticSampleID <- gsub('-GEN','',map.name$geneticSampleID)
    rownames(obs) <- map.name$geneticSampleID
  }
  for(k in 1:length(lev.cast)){
    lev.cast[[k]] <- lev.cast[[k]][rownames(lev.cast[[k]]) %in% rownames(obs),]
  }
  x        <- lev.cast$mean    [,namey[i]][order(lev.cast$mean[,namey[i]])]
  ci_0.975 <- lev.cast$ci_0.975[,namey[i]][order(match(names(lev.cast$ci_0.975[,namey[i]]),names(x)))]
  ci_0.025 <- lev.cast$ci_0.025[,namey[i]][order(match(names(lev.cast$ci_0.025[,namey[i]]),names(x)))]
  pi_0.975 <- lev.cast$pi_0.975[,namey[i]][order(match(names(lev.cast$pi_0.975[,namey[i]]),names(x)))]
  pi_0.025 <- lev.cast$pi_0.025[,namey[i]][order(match(names(lev.cast$pi_0.025[,namey[i]]),names(x)))]
  obs      <- obs[rownames(obs) %in% names(x),,drop=F]
  y        <- obs[order(match(rownames(obs),names(x))),]
  
  #get y-limit.
  obs_limit <- max(y)*1.1
  if(max(pi_0.975) > as.numeric(obs_limit)){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  if(i == 1 | i == 2){y_max <- 1.15}
  y_min <- min(pi_0.025 * 0.95)
  if(y_min < 0.05){y_min <- 0}
  limy <- c(y_min, y_max)

  #plot
  plot(y ~ x, cex = core.cex, pch=glob.pch, ylim=limy, ylab=NA, xlab = NA, col = 'black', bty = 'l')
  mod_fit <- lm(y ~ x)
  rsq <- round(summary(mod_fit)$r.squared,2)
  rsq.1 <- round(rsq_1.1(y, x), 2)
  if(rsq.1 < 0){rsq.1 <- 0}
  rsq.lab <- bquote(R^2 == .(rsq))
  rsq.1.lab <- bquote({R^2} [1:1] == .(rsq.1))
  mtext(rsq.lab  , side = 3, line = -2.6, adj = 0.03, cex = rsq.lab.cex)
  mtext(rsq.1.lab, side = 3, line = -4.1, adj = 0.03, cex = rsq.lab.cex)
  #uppercase first letter.
  laby <- paste0(toupper(substr(namey[i], 1, 1)), substr(namey[i], 2, nchar(namey[i])))
  if(laby == 'Oligotroph'){laby <- 'Oligotrophic Bacteria'}
  mtext(laby, side = 2, line = 2.5, cex = outer.cex)
  #add confidence interval.
  polygon(c(x, rev(x)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0) #predictive interval.
  polygon(c(x, rev(x)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0) #credible interval.
  #fraction within 95% predictive interval.
  in_it <- round(sum(as.numeric(y) < pi_0.975 & as.numeric(y) > pi_0.025) / length(y),2) * 100
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.2, adj = 0.05, cex = rsq.lab.cex)
  abline(0,1,lwd=2)
  abline(mod_fit, lwd =2, lty = 2, col = bf_col)
  
  
  #plot.level.----
  #organize data.
  lev.cast <- fcast[[level[i]]]$plot.fit
  obs <- list(plot_mu[[i]],plot_lo95[[i]],plot_hi95[[i]])
  names(obs) <- c('mu','lo95','hi95')
  for(k in 1:length(obs)){
    rownames(obs[[k]]) <- gsub('.','_',rownames(obs[[k]]), fixed = T)
    obs[[k]] <- obs[[k]][rownames(obs[[k]]) %in% rownames(lev.cast$mean), ,drop = F]
  }
  for(k in 1:length(lev.cast)){
    lev.cast[[k]] <- lev.cast[[k]][rownames(lev.cast[[k]]) %in% rownames(obs$mu),]
  }
  x        <- lev.cast$mean    [,namey[i]][order(lev.cast$mean[,namey[i]])]
  ci_0.975 <- lev.cast$ci_0.975[,namey[i]][order(match(names(lev.cast$ci_0.975[,namey[i]]),names(x)))]
  ci_0.025 <- lev.cast$ci_0.025[,namey[i]][order(match(names(lev.cast$ci_0.025[,namey[i]]),names(x)))]
  pi_0.975 <- lev.cast$pi_0.975[,namey[i]][order(match(names(lev.cast$pi_0.975[,namey[i]]),names(x)))]
  pi_0.025 <- lev.cast$pi_0.025[,namey[i]][order(match(names(lev.cast$pi_0.025[,namey[i]]),names(x)))]
  #make sure row order matches.
  for(k in 1:length(obs)){
    obs[[k]] <- obs[[k]][order(match(rownames(obs[[k]]), names(x))),,drop = F]
  }
  obs.mu   <- obs$mu  [,namey[i]]
  obs.lo95 <- obs$lo95[,namey[i]]
  obs.hi95 <- obs$hi95[,namey[i]]
  names(obs.mu)   <- names(x)
  
  #Make out_sites sites gray for out_spp.
  obs.cols <- rep('black',length(obs.mu))

  #get y-limit.
  obs_limit <- max(obs.hi95)
  if(max(pi_0.975) > obs_limit){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  if(i == 2){y_max <- 1.15}
  y_min <- min(pi_0.025 * 0.95)
  if(y_min < 0.05){y_min <- 0}
  limy <- c(y_min, y_max)
  
  
  #plot
  plot(obs.mu ~ x, cex = plot.cex, pch=glob.pch, ylim=limy, ylab=NA, xlab = NA, col = obs.cols, bty = 'l')
  arrows(c(x), obs.lo95, c(x), obs.hi95, length=0.00, angle=90, code=3, col = obs.cols)
  mod_fit <- lm(obs.mu ~ x)
  rsq <- round(summary(mod_fit)$r.squared,2)
  rsq.1 <- round(rsq_1.1(obs.mu, x), 2)
  if(rsq.1 < 0){rsq.1 <- 0}
  rsq.lab <- bquote(R^2 == .(rsq))
  rsq.1.lab <- bquote({R^2} [1:1] == .(rsq.1))
  mtext(rsq.lab  , side = 3, line = -2.6, adj = 0.03, cex = rsq.lab.cex)
  mtext(rsq.1.lab, side = 3, line = -4.1, adj = 0.03, cex = rsq.lab.cex)
  #1:1 line
  abline(0,1, lwd = 2)
  abline(mod_fit, lwd =2, lty = 2, col = bf_col)
  
  #add confidence interval.
  polygon(c(x, rev(x)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(x, rev(x)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.2, adj = 0.05, cex = rsq.lab.cex)
  
  
  #site.level.----
  #organize data.
  lev.cast <- fcast[[level[i]]]$site.fit
  obs <- list(site_mu[[i]],site_lo95[[i]],site_hi95[[i]])
  names(obs) <- c('mu','lo95','hi95')
  for(k in 1:length(obs)){
    rownames(obs[[k]]) <- gsub('.','_',rownames(obs[[k]]), fixed = T)
    obs[[k]] <- obs[[k]][rownames(obs[[k]]) %in% rownames(lev.cast$mean),, drop=F]
  }
  for(k in 1:length(lev.cast)){
    lev.cast[[k]] <- lev.cast[[k]][rownames(lev.cast[[k]]) %in% rownames(obs$mu),]
  }
  x        <- lev.cast$mean[,namey[i]][order(lev.cast$mean[,namey[i]])]
  ci_0.975 <- lev.cast$ci_0.975[,namey[i]][order(match(names(lev.cast$ci_0.975[,namey[i]]),names(x)))]
  ci_0.025 <- lev.cast$ci_0.025[,namey[i]][order(match(names(lev.cast$ci_0.025[,namey[i]]),names(x)))]
  pi_0.975 <- lev.cast$pi_0.975[,namey[i]][order(match(names(lev.cast$pi_0.975[,namey[i]]),names(x)))]
  pi_0.025 <- lev.cast$pi_0.025[,namey[i]][order(match(names(lev.cast$pi_0.025[,namey[i]]),names(x)))]
  #make sure row order matches.
  for(k in 1:length(obs)){
    obs[[k]] <- obs[[k]][order(match(rownames(obs[[k]]), names(x))),,drop = F]
  }
  obs.mu   <- obs$mu  [,namey[i]]
  obs.lo95 <- obs$lo95[,namey[i]]
  obs.hi95 <- obs$hi95[,namey[i]]

  #Make out_sites sites gray for out_spp.
  obs.cols <- rep('black',length(obs.mu))

  #get y-limit.
  obs_limit <- max(obs.hi95)
  if(max(pi_0.975) > obs_limit){obs_limit <- max(pi_0.975)}
  y_max <- as.numeric(obs_limit)*1.05
  if(y_max > 0.95){y_max <- 1}
  y_min <- min(pi_0.025 * 0.95)
  if(y_min < 0.05){y_min <- 0}
  limy <- c(y_min, y_max)
  
  #plot
  plot(obs.mu ~ x, cex = site.cex, pch=glob.pch, ylim=limy, ylab=NA, xlab = NA, col = obs.cols, bty = 'l')
  arrows(c(x), obs.lo95, c(x), obs.hi95, length=0.0, angle=90, code=3, col = obs.cols)
  mod_fit <- lm(obs.mu ~ x)
  rsq <- round(summary(mod_fit)$r.squared,2)
  rsq.1 <- round(rsq_1.1(obs.mu, x), 2)
  if(rsq.1 < 0){rsq.1 <- 0}
  rsq.lab <- bquote(R^2 == .(rsq))
  rsq.1.lab <- bquote({R^2} [1:1] == .(rsq.1))
  mtext(rsq.lab  , side = 3, line = -2.6, adj = 0.03, cex = rsq.lab.cex)
  mtext(rsq.1.lab, side = 3, line = -4.1, adj = 0.03, cex = rsq.lab.cex)
  #1:1 line
  abline(0,1, lwd = 2)
  abline(mod_fit, lwd =2, lty = 2, col = bf_col)
  
  #add confidence interval.
  polygon(c(x, rev(x)),c(pi_0.975, rev(pi_0.025)), col=adjustcolor('green', trans), lty=0)
  polygon(c(x, rev(x)),c(ci_0.975, rev(ci_0.025)), col=adjustcolor('blue' , trans), lty=0)
  #fraction within 95% predictive interval.
  in_it <- round(sum(obs.mu < pi_0.975 & obs.mu > pi_0.025) / length(obs.mu),2) *100
  state <- paste0(in_it,'% of obs. within interval.')
  mtext(state,side = 3, line = -1.2, adj = 0.05, cex = rsq.lab.cex)
}

#outer labels.----
mtext('core-level', side = 3, line = -1.8, adj = 0.11, cex = outer.cex, outer = T)
mtext('plot-level', side = 3, line = -1.8, adj = 0.50, cex = outer.cex, outer = T)
mtext('site-level', side = 3, line = -1.8, adj = 0.88, cex = outer.cex, outer = T)
mtext( 'observed relative abundance', side = 2, line = 3, cex = outer.cex, outer = T)
mtext('predicted relative abundance', side = 1, line = 2, cex = outer.cex, outer = T)

#end plot.----
dev.off()

