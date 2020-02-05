#Comparing R2 density between Tedersoo out-of-sample (OOS) forecast to NEON and NEON cross-validation (CV).
rm(list=ls())
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')
source('NEFI_functions/rsq_1.1.r')

#set figure out put path.----
output.path <- 'test.png'

#Load calibration and validation predictions and data for OOS and CV.----
#Out of sample Tedersoo calibration, and out of sample predictions for NEON.
#oos.val <- readRDS(NEON_dmulti.ddirch_fcast_fg.path)

#NEON out of sample observations across scale.
#oos.dat <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)

#Cross validation forecasts
cv.val.core <- readRDS(core.CV_NEON_fcast.path)
cv.val.plot <- readRDS(plot.CV_NEON_fcast.path)

#Cross validation predictions
#cv.cal.core <- readRDS(core.CV_NEON_dmulti.ddirch_all.path)
#cv.cal.plot <- readRDS(plot.CV_NEON_dmulti.ddirch_all.path)


#Cross validation data (calibration and validation)
cv.dat.core <- readRDS(core.CV_NEON_cal.val_data.path)
cv.dat.plot <- readRDS(plot.CV_NEON_cal.val_data.path)
cv.dat.core <- cv.dat.core$val$y.val
cv.dat.plot <- cv.dat.plot$val$y.val


#Get NEON Cross-validated stats at core and plot level.-----
#core level.
core.list <- list()
for(i in 1:length(cv.val.core)){
  x <- cv.val.core[[i]]$core.fit$mean
  y <- cv.dat.core[[i]]$rel.abundances
  #deal with trailing -GEN
  rownames(y) <- gsub('-GEN','', rownames(y))
  rownames(x) <- gsub('-GEN','', rownames(x)) 
  #match row and column names.
  y <- y[,order(match(colnames(y), colnames(x)))]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[rownames(y) %in% rownames(x),]
  y <- y[order(match(rownames(y), rownames(x))),]
  #Get R2 for each column.
  lev.stats <- list()
  for(j in 1:ncol(y)){
    #r square best fit.
    fit <- lm(y[,j] ~ x[,j])
    rsq <- summary(fit)$r.squared
    #r square 1:1.
    rss <- sum((x[,j] -      y[,j])  ^ 2)  ## residual sum of squares
    tss <- sum((y[,j] - mean(y[,j])) ^ 2)  ##    total sum of squares
    rsq1 <- 1 - rss/tss
    rsq1 <- rsq_1.1(y[,j], x[,j])          ## new rsq_1.1 function.
    if(rsq1 < 0){rsq1 <- 0}
    #RMSE.
    rmse <- sqrt(mean(fit$residuals^2))
    return <- c(rsq,rsq1,rmse)
    lev.stats[[j]] <- return
  }
  lev.stats <- do.call(rbind, lev.stats)
  colnames(lev.stats) <- c('rsq','rsq1','rmse')
  rownames(lev.stats) <- colnames(y)
  lev.stats <- lev.stats[-(rownames(lev.stats) %in% 'other'),]
  core.list[[i]] <- lev.stats
}
names(core.list) <- names(cv.val.core)

#plot level.
plot.list <- list()
for(i in 1:length(cv.val.plot)){
  x <- cv.val.plot[[i]]$plot.fit$mean
  y <- cv.dat.plot[[i]]$mean
  #deal with trailing -GEN
  rownames(y) <- gsub('-GEN','', rownames(y))
  rownames(x) <- gsub('-GEN','', rownames(x)) 
  #match row and column names.
  y <- y[,order(match(colnames(y), colnames(x)))]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[rownames(y) %in% rownames(x),]
  y <- y[order(match(rownames(y), rownames(x))),]
  #Get R2 for each column.
  lev.stats <- list()
  for(j in 1:ncol(y)){
    #r square best fit.
    fit <- lm(y[,j] ~ x[,j])
    rsq <- summary(fit)$r.squared
    #r square 1:1.
    rss <- sum((x[,j] -      y[,j])  ^ 2)  ## residual sum of squares
    tss <- sum((y[,j] - mean(y[,j])) ^ 2)  ## total sum of squares
    rsq1 <- 1 - rss/tss
    rsq1 <- rsq_1.1(y[,j], x[,j])          ## new rsq_1.1 function.
    if(rsq1 < 0){rsq1 <- 0}
    #RMSE.
    rmse <- sqrt(mean(fit$residuals^2))
    return <- c(rsq,rsq1,rmse)
    lev.stats[[j]] <- return
  }
  lev.stats <- do.call(rbind, lev.stats)
  colnames(lev.stats) <- c('rsq','rsq1','rmse')
  rownames(lev.stats) <- colnames(y)
  lev.stats <- lev.stats[-(rownames(lev.stats) %in% 'other'),]
  plot.list[[i]] <- lev.stats
}
names(plot.list) <- names(cv.val.plot)


#Get means by taxaonomic level.----
#core scale.
core.lev.mu <- list()
for(i in 1:length(core.list)){
  z <- data.frame(core.list[[i]])
  core.lev.mu[[i]] <- mean(z$rsq, na.rm = T)
}
core.lev.mu <- unlist(core.lev.mu)
names(core.lev.mu) <- names(core.list)
names(core.lev.mu)[length(core.lev.mu)] <- 'functional'
ref.order <- c(length(core.lev.mu), 1:(length(core.lev.mu) - 1))
core.lev.mu <- core.lev.mu[ref.order]

#plot scale.
plot.lev.mu <- list()
for(i in 1:length(plot.list)){
  z <- data.frame(plot.list[[i   ]])
  plot.lev.mu[[i]] <- mean(z$rsq, na.rm = T)
}
plot.lev.mu <- unlist(plot.lev.mu)
names(plot.lev.mu) <- names(plot.list)
names(plot.lev.mu)[length(plot.lev.mu)] <- 'functional'
ref.order <- c(length(plot.lev.mu), 1:(length(plot.lev.mu) - 1))
plot.lev.mu <- plot.lev.mu[ref.order]

#put together.
to.plot <- cbind(core.lev.mu, plot.lev.mu)

#png save line.-----
#png(filename=output.path,width=10,height=8,units='in',res=300)

#Begin plot: Global plot settings.-----
par(mfrow = c(1,2),
    mar = c(4,4,2,2))
limx <- c(0,1)
limy <- c(0, 25)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
y.cex <- 1.3

#ITS by spatial and taxonomic scale.----
cols <- brewer.pal(length(core.lev.mu) - 1,'Spectral')
cols <- c(cols,'green')
limy <- c(0, max(c(core.lev.mu, plot.lev.mu))*1.1)
x <- c(1:2)
limx <- c(min(x)*0.95, max(x)*1.05)
plot(to.plot[1,] ~ x, xlim = limx, ylim = limy, bty = 'l', col = 'black', bg = cols[1], 
     pch = 21, cex = 2.5, xlab= NA, ylab = NA, xaxt = 'n')
lines(x, to.plot[1,], lty = 1, col = cols[1], lwd = 1.5)
for(i in 2:nrow(to.plot)){
  points(to.plot[i,] ~ x, pch = 21, cex = 2.5, col = 'black', bg = cols[i])
  lines(x, to.plot[i,], lty = 1, lwd = 1.5,    col = cols[i])
}
axis(1, at = c(1,2,3), labels = F)
lab <- c('core','plot')
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= lab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
mtext(expression(paste("Cross-Validated R"^"2")), side = 2, line = 2.5, cex = y.cex)
leg.lab <- rownames(to.plot)
legend('topleft',leg.lab, pch = 21, col = 'black', pt.bg = cols, ncol = 1, bty= 'n')
mtext('(a)', side = 3, adj = 0.98, line = -2)


#end plot.-----
#dev.off()
