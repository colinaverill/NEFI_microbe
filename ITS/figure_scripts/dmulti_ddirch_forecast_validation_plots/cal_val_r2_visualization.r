#plott calibration and validation R2 distributions for Fungi and bacteria as 6-panel plot.
rm(list=ls())
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')

#set output path.----
output.path <- 'test.png'

#Load calibration data.----
cal <- readRDS(ted_ITS.prior_dmulti.ddirch_fg_JAGSfit) #all phylo and functional groups.
cal <- readRDS(ted_ITS_prior_all.groups_JAGSfits.path)
#re-order list, make function group first.
cal <- cal[c('fg','phylum','class','order','family','genus')]
names(cal)[1] <- 'functional'

#Get calibration r2 values.----
cal.r2 <- list()
lev.cal.out <- list()
for(i in 1:length(cal)){
  lev <- cal[[i]]
  obs <- lev$observed
  pred <- lev$predicted
  lev.r2 <- list()
  for(j in 1:ncol(obs)){lev.r2[[j]] <- summary(lm(obs[,j] ~ pred[,j]))$r.squared}
  lev.r2 <- unlist(lev.r2)
  names(lev.r2) <- colnames(obs)
  lev.r2 <- lev.r2[names(lev.r2) != 'other']
  cal.r2[[i]] <- lev.r2
  #names(pl.r2)[[i]] <- names(pl)[i]
}
names(cal.r2) <- names(cal)
#Subset to ones we can actually predict (r2 > 0.1).
lev.cal.r2 <- list()
for(i in 1:length(cal.r2)){
  lev.cal.r2[[i]] <- cal.r2[[i]][cal.r2[[i]] > 0.1]
}
cal.mu <- unlist(lapply(lev.cal.r2, mean  ))
cal.sd <- unlist(lapply(lev.cal.r2, sd    ))
cal.N  <- unlist(lapply(lev.cal.r2, length))
cal.se <- cal.sd / sqrt(cal.N)
names(cal.mu) <- names(cal.r2)
lev.cal.r2 <- unlist(lev.cal.r2)


#load forecasts predicted and observed.----
#val.cast <- readRDS(NEON_dmulti.ddirch_fcast_fg.path) 
val.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)
val.cast <- readRDS(NEON_site_fcast_all_groups_1k_rare.path)
#pl.core <- readRDS(NEON_ITS_fastq_all_cosmo_phylo_groups_1k_rare.path)

#Re-order.
val.cast  <- val.cast [c('fg','phylum','class','order','family','genus')]
val.truth <- val.truth[c('fg','phylum','class','order','family','genus')]
#pl.core <- pl.core[c('fg','phylum','class','order','family','genus')]

#names.
names(val.cast)[1] <- 'functional'

#get core, plot site R2 values out of sample.----
all.core.rsq <- list()
all.plot.rsq <- list()
all.site.rsq <- list()
for(i in 1:length(val.cast)){
  fcast <- val.cast[[i]]
  core.rsq <- list()
  plot.rsq <- list()
  site.rsq <- list()
  #core.level----
  #y <- (pl.core[[i]]$abundances + 1)/pl.core[[i]]$seq_total
  y <- val.truth[[i]]$core.fit
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
    if(fungi_name == 'Ectomycorrhizal'){
      sub.y <- y[-grep('DSNY',rownames(y)),]
      sub.x <- x[-grep('DSNY',rownames(x)),]
      rsq <- summary(lm(sub.y[,k] ~ sub.x[,k]))$r.squared
    }
    names(rsq) <- fungi_name
    if(fungi_name == 'other'){next}
    core.rsq[[k]] <- rsq
  }
  #plot.level----
  x <- fcast$plot.fit$mean
  y <- val.truth[[i]]$plot.fit$mean
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
    if(fungi_name == 'Ectomycorrhizal'){
      sub.y <- y[-grep('DSNY',rownames(y)),]
      sub.x <- x[-grep('DSNY',rownames(x)),]
      rsq <- summary(lm(sub.y[,k] ~ sub.x[,k]))$r.squared
    }
    names(rsq) <- fungi_name
    if(fungi_name == 'other'){next}
    plot.rsq[[k]] <- rsq
  }
  #site.level----
  x <- fcast$site.fit$mean
  y <- val.truth[[i]]$site.fit$mean
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
    if(fungi_name == 'Ectomycorrhizal'){
      sub.y <- y[-grep('DSNY',rownames(y)),]
      sub.x <- x[-grep('DSNY',rownames(x)),]
      rsq <- summary(lm(sub.y[,k] ~ sub.x[,k]))$r.squared
    }
    names(rsq) <- fungi_name
    if(fungi_name == 'other'){next}
    site.rsq[[k]] <- rsq
  }
  #wrap up for return.----
  all.core.rsq[[i]] <- unlist(core.rsq)
  all.plot.rsq[[i]] <- unlist(plot.rsq)
  all.site.rsq[[i]] <- unlist(site.rsq)
}

val.mu <- unlist(lapply(all.site.rsq, mean  ))
val.sd <- unlist(lapply(all.site.rsq, sd    ))
val.N  <- unlist(lapply(all.site.rsq, length))
val.se <- val.sd / sqrt(val.N)
names(val.mu) <- names(val.cast)
core.rsq <- unlist(all.core.rsq)
plot.rsq <- unlist(all.plot.rsq)
site.rsq <- unlist(all.site.rsq)
#core.rsq <- core.rsq[-grep('other',names(core.rsq))]
#plot.rsq <- plot.rsq[-grep('other',names(plot.rsq))]
#site.rsq <- site.rsq[-grep('other',names(site.rsq))]


#Subset to observations that have a minimum calibration R2 value.----
pass <- lev.cal.r2
core.rsq <- core.rsq[names(core.rsq) %in% names(pass)]
plot.rsq <- plot.rsq[names(plot.rsq) %in% names(pass)]
site.rsq <- site.rsq[names(site.rsq) %in% names(pass)]
core.d <- zero_truncated_density(core.rsq)
plot.d <- zero_truncated_density(plot.rsq)
site.d <- zero_truncated_density(site.rsq)

#png save line.----
png(filename=output.path,width=8,height=6,units='in',res=300)

#global plot settings.----
par(mfrow = c(2,3), mar = c(5,4,2,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.0 #outer label size.
cols <- c('purple','cyan','yellow')

#Calibration R2 denisty plot.-----
dat <- unlist(cal.r2)
dat <- zero_truncated_density(dat)
limy <- c(0, max(dat$y)*1.05)
limx <- c(0, max(dat$x))
#Density plot.
plot(dat,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(dat, col = adjustcolor(cols[3],trans))
mtext('Density', side = 2, line = 2.5, cex = o.cex)
mtext(expression(paste("Calibration R"^"2")), side = 1, line = 3, cex = o.cex)
mtext('(a)', side = 3, adj = 0.95, line = -2)

#Validation R2 denisty plot.-----
dat.a <- core.d
dat.b <- plot.d
dat.c <- site.d
limy <- c(0, max(c(dat.a$y, dat.b$y, dat.c$y))*1.05)
limx <- c(0, max(c(dat.a$x, dat.b$x, dat.c$x)))
#Density plot.
plot(dat.a,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(dat.a, col = adjustcolor(cols[1],trans))
polygon(dat.b, col = adjustcolor(cols[2],trans))
polygon(dat.c, col = adjustcolor(cols[3],trans))
mtext('Density', side = 2, line = 2.5, cex = o.cex)
mtext(expression(paste("Validation R"^"2")), side = 1, line = 3, cex = o.cex)
#legend
legend(x = 0.5, y = 12.5, legend = c('core','plot','site'), col ='black', pt.bg=adjustcolor(cols,trans), bty = 'n', pch = 22, pt.cex = 1.5, cex = 1.2)
mtext('(b)', side = 3, adj = 0.95, line = -2)


#Calibration and Validation rsq ~ function/phylo scale.----
x <- 1:length(cal.mu)
limy <- c(0,max(cal.mu)*1.1)
plot(cal.mu ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu, lty = 2)
lines(x, val.mu, lty = 2)
points(val.mu ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level R"^"2")), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = -.03, labels= names(cal.mu), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
legend(x = 1, y = 0.1, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(c)', side = 3, adj = 0.95, line = -2)

#Outer labels.----
mtext('Fungi'   , cex = 1.5, side = 3, outer = T, adj = 0.01, line = -1)
mtext('Bacteria', cex = 1.5, side = 3, outer = T, adj = 0.01, line = -23)

#end plot.----
dev.off()
