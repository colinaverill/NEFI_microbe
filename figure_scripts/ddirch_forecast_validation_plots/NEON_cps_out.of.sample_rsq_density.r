#Density plot of out of sample R2 values as a function of spatial scale (core-plot-site).
rm(list=ls())
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')

#set output path.----
output.path <- 'test.png'
output.path <- NEON_cps_out.of.sample_rsq_density_figure.path

#load in sample fits.----
fg.cal <- readRDS(ted_ITS.prior_fg_JAGSfit)
pl.cal <- readRDS(ted_ITS_prior_phylo.group_JAGSfits)
pl.cal$function_group <- fg.cal

#get in sample r2 values.----
cal.rsq <- list()
for(i in 1:length(pl.cal)){
  if(names(pl.cal)[i] != 'function_group'){
    y <- pl.cal[[i]]$observed
    x <- pl.cal[[i]]$predicted
  }
  if(names(pl.cal)[i] == 'function_group'){
    y <- pl.cal[[i]]$all.preds$observed
    x <- pl.cal[[i]]$all.preds$predicted
  }
  lev.rsq <- list()
  for(k in 1:ncol(y)){
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- colnames(y)[k]
    lev.rsq[[k]] <- rsq
  }
  lev.rsq <- unlist(lev.rsq)
  lev.rsq <- lev.rsq[-grep('other',names(lev.rsq))]
  cal.rsq[[i]] <- lev.rsq
}
cal.rsq <- unlist(cal.rsq)

#load forecasts predicted and observed.----
#functional groups
fg.cast <- readRDS(NEON_site_fcast_fg.path)
fg.core <- readRDS(NEON_ITS_fastq_taxa_fg.path)
#trick out to make work with phylo stuff.
fg.core$abundances$geneticSampleID <- NULL
fg.plot.truth <- readRDS(NEON_plot.level_fg_obs_fastq.path)
fg.site.truth <- readRDS(NEON_site.level_fg_obs_fastq.path)
#format fg stuff so we can pop it into the pl list.
fg.truth <- list(fg.plot.truth,fg.site.truth)
names(fg.truth) <- c('plot.fit','site.fit')

#phylogenetic groups.
pl.cast <- readRDS(NEON_site_fcast_all_phylo_levels.path)
pl.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq.path)
pl.core <- readRDS(NEON_ITS_fastq_all_cosmo_phylo_groups.path)
#merge in fg stuff
pl.core$function_group <- fg.core
pl.truth$function_group <- fg.truth
pl.cast$function_group <- fg.cast


#get core, plot site R2 values out of sample.----
all.core.rsq <- list()
all.plot.rsq <- list()
all.site.rsq <- list()
for(i in 1:length(pl.cast)){
  fcast <- pl.cast[[i]]
  core.rsq <- list()
  plot.rsq <- list()
  site.rsq <- list()
  #core.level----
  y <- (pl.core[[i]]$abundances + 1)/pl.core[[i]]$seq_total
  x <- fcast$core.fit$mean
  #make sure row and column orders match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(fcast$core.fit$mean)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    core.rsq[[k]] <- rsq
  }
  #plot.level----
  x <- fcast$plot.fit$mean
  y <- pl.truth[[i]]$plot.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(y)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    plot.rsq[[k]] <- rsq
  }
  #site.level----
  x <- fcast$site.fit$mean
  y <- pl.truth[[i]]$site.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(y)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    site.rsq[[k]] <- rsq
  }
  #wrap up for return.----
  all.core.rsq[[i]] <- unlist(core.rsq)
  all.plot.rsq[[i]] <- unlist(plot.rsq)
  all.site.rsq[[i]] <- unlist(site.rsq)
  
}
core.rsq <- unlist(all.core.rsq)
plot.rsq <- unlist(all.plot.rsq)
site.rsq <- unlist(all.site.rsq)
core.rsq <- core.rsq[-grep('other',names(core.rsq))]
plot.rsq <- plot.rsq[-grep('other',names(plot.rsq))]
site.rsq <- site.rsq[-grep('other',names(site.rsq))]

#Subset to observations that have a minimum calibration R2 value.----
pass <- cal.rsq[cal.rsq > .1]
core.rsq <- core.rsq[names(core.rsq) %in% names(pass)]
plot.rsq <- plot.rsq[names(plot.rsq) %in% names(pass)]
site.rsq <- site.rsq[names(site.rsq) %in% names(pass)]
core.d <- zero_truncated_density(core.rsq)
plot.d <- zero_truncated_density(plot.rsq)
site.d <- zero_truncated_density(site.rsq)

#png save line.----
png(filename=output.path,width=5,height=5,units='in',res=300)

#global plot settings.----
par(mfrow = c(1,1))
limx <- c(0,1)
limy <- c(0, 21)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
cols <- c('purple','cyan','yellow')
par(mfrow = c(1,1), mar = c(4.5,4,1,1))

#plot.----
plot(core.d,xlim = c(0, 0.8), ylim = limy, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(core.d, col = adjustcolor(cols[1],trans))
polygon(plot.d, col = adjustcolor(cols[2],trans))
polygon(site.d, col = adjustcolor(cols[3],trans))
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Validation R"^"2")), side = 1, line = 2.5, cex = o.cex)
legend(x = 0.6, y = 16, legend = c('core','plot','site'), col ='black', pt.bg=adjustcolor(cols,trans), bty = 'n', pch = 22, pt.cex = 1.5)

#end plot.----
dev.off()
