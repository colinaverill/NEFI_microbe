#Bray-Curtis community similarity ~ spac across phyogenetic scales.
rm(list=ls())
source('paths.r')
library(mgcv)

#Load data.----
d_all <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)
loc_all <- readRDS(dp1.10086.00_output.path)
#reorder and rename.
d_all <- d_all[c('fg','phylum','class','order','family','genus')]
names(d_all)[1] <- 'functional'


#Generate Bray-Curtis similarity matrices and spatial matrices.----
sim <- list()
spa <- list()
for(k in 1:length(d_all)){
  #Grab tax/functional level.----
  d <- d_all[[k]]
  d <- d$core.fit
  loc <- loc_all

  #match up data, subset to mineral soil, drop stuff from identical locations.----
  rownames(d) <- gsub('-GEN','',rownames(d))
  #drop organic horizons. Generates zero distances.
  loc <- loc[loc$horizon == 'M',]
  loc$lat_lon <- paste0(loc$adjDecimalLongitude, loc$adjDecimalLatitude)
  drop <- loc[duplicated(loc$lat_lon),]$sampleID #drop some duplicated locations.
  loc <- loc[!(loc$sampleID %in% drop),]
  loc <- loc[loc$sampleID %in% rownames(d),]
  d <-   d[rownames(d) %in% loc$sampleID,]
  loc <- loc[order(match(loc$sampleID, rownames(d))),]

  #Generate spatial distance matrix.----
  tax.loc <- loc[,c('adjDecimalLongitude','adjDecimalLatitude')]
  tax.dist <- geosphere::distm(tax.loc)
  tax.dist <- tax.dist[lower.tri(tax.dist)]
  
  #Generate Bray-Curtis similarity matrix.----
  bray <- 1 - vegan::vegdist(d)
  bray <- as.matrix(bray)
  bray <- bray[lower.tri(bray)]

  #return output.----
  sim[[k]] <- bray
  spa[[k]] <- tax.dist
}
names(sim) <- names(d_all)
names(spa) <- names(d_all)

#plot the results!----
par(mfrow = c(3,2))
for(i in 1:length(sim)){
  y <- log10(sim[[i]])
  x <- log10(spa[[i]])
  
  #fit lines within scales.
  dat <- data.frame(y, x)
  #dat <- dat[sample(nrow(dat), 1000),]
  core.fit <- lm(y ~ x, dat[dat$x < 2,])
  plot.fit <- lm(y ~ x, dat[dat$x > 2 & dat$x < 4.8,])
  site.fit <- lm(y ~ x, dat[dat$x > 4.8,])
  
  #plot similarity.
  group <- names(sim)[i]
  cat('plotting',group,'...\n')
  plot(y ~ x, cex = 0.1, ylab = NA, xlab = NA, data = dat)
  mtext(group, side = 3, adj = 0.05)
  
  #plot line across scales.
  abline(lm(y ~ x), lwd = 2, col = 'light gray')
  
  #gam 
  cat('Fitting gam...\n')
  x.gam <- seq(min(x), max(x), by = 0.1)
  newdat <- data.frame(x.gam)
  colnames(newdat) <- c('x')
  mod.gam <- gam(y ~ s(x, k=-1))
  gam.pred <- predict(mod.gam, newdata = newdat)
  lines(smooth.spline(gam.pred ~ x.gam), lwd = 2, col = 'green')
  
  #plot lines within scales.
  x.core <- seq(0,2,by = 0.1)
  y.core <- coef(core.fit)[1] + coef(core.fit)[2] * x.core
  lines(smooth.spline(x.core, y.core), lwd = 2, col = 'purple')
  x.plot <- seq(2,4.8, by = 0.2)
  y.plot <- coef(plot.fit)[1] + coef(plot.fit)[2]*x.plot
  lines(smooth.spline(x.plot, y.plot), lwd = 2, col = 'purple')
  x.site <- seq(4.8, 6.45, by = 0.1)
  y.site <- coef(site.fit)[1] + coef(site.fit)[2]*x.site
  lines(smooth.spline(x.site, y.site), lwd = 2, col = 'purple')
}

#outer labels.
mtext("log10(Bray-Curtis Similarity)", line = -2, side = 2, outer = T)
mtext("log10(meters)",side = 1, outer = T, line = -2)

