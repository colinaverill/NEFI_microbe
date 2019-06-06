#calculating spatial signal in taxa by phylo scale across NEON, using entire network.
rm(list=ls())
source('paths.r')
library(boot)

#output path.----
output.path <- 'Moran_I_figure.png'

#logit transform observed values and model residuals?----
do_logit <- F

#Load data.----
d_all <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)
loc_all <- readRDS(dp1.10086.00_output.path)
#spatial forecast.
fcast <- readRDS(NEON_site_fcast_all_groups_1k_rare.path)

#Calculate spatial statistics for all levels on forecast residuals.----
a_out <- list() #spatial signal of raw values.
b_out <- list() #spatial signal of model residuals.
for(k in 1:length(d_all)){
  #Grab tax/functional level.
  d <- d_all[[k]]
  d <- d$core.fit
  loc <- loc_all
  pred <- fcast[[k]]$core.fit$mean
  
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
  pred <- pred[rownames(pred) %in% rownames(d),]
  pred <- pred[order(match(rownames(pred), rownames(d))),]
  #make sure predictions and observation columns are in the same order.
  pred <- pred[,order(match(colnames(pred), colnames(d)))]
  
  #generate model residuals as a logit difference.
  resid <- (d) - (pred)
  if(do_logit == T){
    resid <- logit(d) - logit(pred)
    d <- logit(d)
  }
  
  #Generate spatial distance matrix.----
  tax.loc <- loc[,c('adjDecimalLongitude','adjDecimalLatitude')]
  tax.dist <- geosphere::distm(tax.loc)
  inv.tax.dist <- 1/tax.dist
  diag(inv.tax.dist) <- 0
  
  #calulate moran statistics.----
  moran_out <- list()
  resid_out <- list()
  for(i in 1:ncol(d)){
    #Get r.abundance distance matrix on logit scale.
    cat('Calculating ',colnames(d)[i],' moran statistics...\n')
    #make sure there is variance in taxa being calculated.
    if(length(unique(d[,i])) == 1){
      moran_out[[i]] <- NA
      next
    }
    #moran stat for raw residuals.
    moran.stat  <- ape::Moran.I((d[,i]), inv.tax.dist)
    moran.stat  <- c(moran.stat$observed, moran.stat$p.value)
    #moran stat for post-model residuals.
    resid.stat <- ape::Moran.I(resid[,i], inv.tax.dist)
    resid.stat <- c(resid.stat$observed, resid.stat$p.value)
    #moran.resid
    names(moran.stat) <- c('moran','p.val')
    names(resid.stat) <- c('moran','p.val')
    moran_out[[i]] <- moran.stat
    resid_out[[i]] <- resid.stat
  }
  
  #Drop taxa that couldn't fit, wrap up level output.----
  names(moran_out) <- colnames(d)
  names(resid_out) <- colnames(resid)
  moran_out <- moran_out[!(is.na(moran_out))]
  resid_out <- resid_out[!(is.na(resid_out))]
  moran_out <- data.frame(do.call(rbind, moran_out))
  resid_out <- data.frame(do.call(rbind, resid_out))
  
  #return output.----
  a_out[[k]] <- moran_out
  b_out[[k]] <- resid_out
  
}
names(a_out) <- names(d_all)[1:5]
names(b_out) <- names(d_all)[1:5]

#calculate average morans I for all groups.----
a_avg <- list()
b_avg <- list()
for(i in 1:length(a_out)){
  #no model.
  a.mu <- mean(a_out[[i]]$moran)
  a.sd <-   sd(a_out[[i]]$moran)
  a.se <- a.sd/sqrt(nrow(a_out[[i]]))
  to_return.a <- c(a.mu, a.sd, a.se)
  names(to_return.a) <- c('mu','sd','se')
  a_avg[[i]] <- to_return.a
  #residuals of model.
  b.mu <- mean(b_out[[i]]$moran)
  b.sd <-   sd(b_out[[i]]$moran)
  b.se <- b.sd/sqrt(nrow(b_out[[i]]))
  to_return.b <- c(b.mu, b.sd, b.se)
  names(to_return.b) <- c('mu','sd','se')
  b_avg[[i]] <- to_return.b
  
}
a_avg <- data.frame(do.call(rbind, a_avg))
b_avg <- data.frame(do.call(rbind, b_avg))
rownames(a_avg) <- names(fcast)
rownames(b_avg) <- names(fcast)

#assign x positions, functional groups first.
a_avg$x <- c(2:(nrow(a_avg)), 1)
b_avg$x <- c(2:(nrow(b_avg)), 1)

#setup figure output.----
png(filename=output.path,width=7,height=5,units='in',res=300)

#plot.-----
par(mfrow=c(1,2), mar = c(5,3,.5,1), oma = c(1,1.5,1,1))
#Raw spatial signal FUNGI.
limy <- c(0, max(a_avg$mu + a_avg$se))
plot(mu ~ x, data = a_avg, cex = 2, pch = 16, ylim = limy,
     ylab = NA, xlab = NA, xaxt = 'n', bty = 'l')
#error bars.
mu <- a_avg$mu
x <- a_avg$x
upr <- mu + a_avg$se
lwr <- mu - a_avg$se
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)
#x-axis.
x.lab <- rownames(a_avg)
x.lab[x.lab == 'fg'] <- 'functional'
axis(1, at=a_avg$x, labels= NA, cex = 1, srt = 45)
text(x= a_avg$x + .12, y = -0.06, labels= x.lab, srt=45, adj=1, xpd=TRUE, cex = 1)
mtext("Moran's I", side = 2, line = 2.5, cex = 1.2)
mtext('Fungi',side = 3, line = -2, cex = 1.7, adj  = 0.95)

#FUTURE BACTERIA PLOT. TURN CEX BACK ON WHEN READY FOR POINTS.
plot(mu ~ x, data = a_avg, cex = 0, pch = 16, ylim = limy,
     ylab = NA, xlab = NA, xaxt = 'n', bty = 'l')
#error bars.
#mu <- a_avg$mu
#x <- a_avg$x
#upr <- mu + a_avg$se
#lwr <- mu - a_avg$se
#arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)
#x-axis.
x.lab <- rownames(a_avg)
x.lab[x.lab == 'fg'] <- 'functional'
axis(1, at=a_avg$x, labels= NA, cex = 1, srt = 45)
text(x= a_avg$x + .12, y = -0.06, labels= x.lab, srt=45, adj=1, xpd=TRUE, cex = 1)
mtext('Bacteria',side = 3, line = -2, cex = 1.7, adj  = 0.95)


#end plot.----
dev.off()
