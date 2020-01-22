#Comparing R2 density between Tedersoo out-of-sample (OOS) forecast to NEON and NEON cross-validation (CV).
rm(list=ls())
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')

#set figure out put path.----
output.path <- 'test.png'

#Load calibration and validation predictions and data for OOS and CV.----
#Out of sample Tedersoo calibration, and out of sample predictions for NEON.
oos.cal <- readRDS(ted_ITS.prior_dmulti.ddirch_fg_JAGSfit)
oos.val <- readRDS(NEON_dmulti.ddirch_fcast_fg.path)

#NEON out of sample observations across scale.
oos.dat <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)

#Cross validation forecasts
cv.val.core <- readRDS(core.CV_NEON_fcast.path)
cv.val.plot <- readRDS(plot.CV_NEON_fcast.path)

#Cross validation predictions
cv.cal.core <- readRDS(core.CV_NEON_dmulti.ddirch_all.path)
cv.cal.plot <- readRDS(plot.CV_NEON_dmulti.ddirch_all.path)
 

#Cross validation data (calibration and validation)
 cv.dat.core <- readRDS(core.CV_NEON_cal.val_data.path)
 cv.dat.plot <- readRDS(plot.CV_NEON_cal.val_data.path)
 cv.dat.core <- cv.dat.core$val$y.val
 cv.dat.plot <- cv.dat.plot$val$y.val


#Generate Tedersoo out of sample fit statistics.----
oos.stats <- list()
for(i in 1:length(oos.val)){
  lev   <- oos.val[[i]]
  lev.y <- oos.dat[[i]]
  x <- list(lev$core.fit$mean, lev$plot.fit$mean, lev$site.fit$mean)
  y <- list(lev.y$core.fit, lev.y$plot.fit$mean, lev.y$site.fit$mean)
  cps.list <- list()
  for(k in 1:length(y)){
    yy <- y[[k]]
    xx <- x[[k]]
    #deal with trailing -GEN
    rownames(yy) <- gsub('-GEN','', rownames(yy))
    rownames(xx) <- gsub('-GEN','', rownames(xx)) 
    #match row and column names.
    yy <- yy[,order(match(colnames(yy), colnames(xx)))]
    xx <- xx[rownames(xx) %in% rownames(yy),]
    yy <- yy[order(match(rownames(yy), rownames(xx))),]
    #Get R2 for each column.
    lev.stats <- list()
    for(j in 1:ncol(yy)){
      #r square best fit.
      fit <- lm(yy[,j] ~ xx[,j])
      rsq <- summary(fit)$r.squared
      #r square 1:1.
      rss <- sum((xx[,j] -      yy[,j])  ^ 2)  ## residual sum of squares
      tss <- sum((yy[,j] - mean(yy[,j])) ^ 2)  ## total sum of squares
      rsq1 <- 1 - rss/tss
      if(rsq1 < 0){rsq1 <- 0}
      #RMSE.
      rmse <- sqrt(mean(fit$residuals^2))
      return <- c(rsq,rsq1,rmse)
      lev.stats[[j]] <- return
    }
    lev.stats <- do.call(rbind, lev.stats)
    colnames(lev.stats) <- c('rsq','rsq1','rmse')
    rownames(lev.stats) <- colnames(yy)
    cps.list[[k]] <- lev.stats
  }
  names(cps.list) <- c('core','plot','site')
  oos.stats[[i]] <- cps.list
}
names(oos.stats) <- names(oos.val)

#reformat tedersoo out of sample forecasts to match cross validation.----
oos.core.list <- list()
oos.plot.list <- list()
for(i in 1:length(oos.stats)){
  oos.core.list[[i]] <- oos.stats[[i]]$core
  oos.plot.list[[i]] <- oos.stats[[i]]$plot
}
names(oos.core.list) <- names(oos.stats)
names(oos.plot.list) <- names(oos.stats)

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
  y <- y[,order(match(colnames(y), colnames(xx)))]
  x <- x[rownames(x) %in% rownames(y),]
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
    if(rsq1 < 0){rsq1 <- 0}
    #RMSE.
    rmse <- sqrt(mean(fit$residuals^2))
    return <- c(rsq,rsq1,rmse)
    lev.stats[[j]] <- return
  }
  lev.stats <- do.call(rbind, lev.stats)
  colnames(lev.stats) <- c('rsq','rsq1','rmse')
  rownames(lev.stats) <- colnames(y)
  core.list[[i]] <- lev.stats
}
names(core.list) <- names(cv.dat.core)

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
    if(rsq1 < 0){rsq1 <- 0}
    #RMSE.
    rmse <- sqrt(mean(fit$residuals^2))
    return <- c(rsq,rsq1,rmse)
    lev.stats[[j]] <- return
  }
  lev.stats <- do.call(rbind, lev.stats)
  colnames(lev.stats) <- c('rsq','rsq1','rmse')
  rownames(lev.stats) <- colnames(y)
  plot.list[[i]] <- lev.stats
}
names(plot.list) <- names(cv.val.plot)


#Collapse to 6x2 data frames.----
cv.core.r2   <- list()
cv.plot.r2   <- list()
cv.core.r2.1 <- list()
cv.plot.r2.1 <- list()
cv.core.rmse <- list()
cv.plot.rmse <- list()
for(i in 1:length(core.list)){
  cv.core.r2  [[i]] <- core.list[[i]][,1]
  cv.plot.r2  [[i]] <- plot.list[[i]][,1]
  cv.core.r2.1[[i]] <- core.list[[i]][,2]
  cv.plot.r2.1[[i]] <- plot.list[[i]][,2]
  cv.core.rmse[[i]] <- core.list[[i]][,3]
  cv.plot.rmse[[i]] <- plot.list[[i]][,3]
}
cv.core.r2 <- unlist(cv.core.r2)
cv.plot.r2 <- unlist(cv.plot.r2)
cv.core.r2.1 <- unlist(cv.core.r2.1)
cv.plot.r2.1 <- unlist(cv.plot.r2.1)
cv.core.rmse <- unlist(cv.core.rmse)
cv.plot.rmse <- unlist(cv.plot.rmse)

#Same for Tedersoo forecast.
oos.core.r2   <- list()
oos.plot.r2   <- list()
oos.core.r2.1 <- list()
oos.plot.r2.1 <- list()
oos.core.rmse <- list()
oos.plot.rmse <- list()
for(i in 1:length(oos.core.list)){
  oos.core.r2  [[i]] <- oos.core.list[[i]][,1]
  oos.plot.r2  [[i]] <- oos.plot.list[[i]][,1]
  oos.core.r2.1[[i]] <- oos.core.list[[i]][,2]
  oos.plot.r2.1[[i]] <- oos.plot.list[[i]][,2]
  oos.core.rmse[[i]] <- oos.core.list[[i]][,3]
  oos.plot.rmse[[i]] <- oos.plot.list[[i]][,3]
}
oos.core.r2   <- unlist(oos.core.r2)
oos.plot.r2   <- unlist(oos.plot.r2)
oos.core.r2.1 <- unlist(oos.core.r2.1)
oos.plot.r2.1 <- unlist(oos.plot.r2.1)
oos.core.rmse <- unlist(oos.core.rmse)
oos.plot.rmse <- unlist(oos.plot.rmse)


#plot R2, R2-1 and RMSE stats by scale.-----
#png save line.
png(filename=output.path,width=10,height=8,units='in',res=300)

#Global plot settings.
par(mfrow = c(2,2),
    mar = c(4,4,2,2))
limx <- c(0,1)
limy <- c(0, 25)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
cols <- c('purple','cyan','yellow')
labs <- c('Global Forecast','NEON Cross-Validation')
labs <- c('Global Forecast (a)','NEON Cross-Validation (b)','Global Forecast (c)','NEON Cross-Validation (d)')

#Forecast OOS R2----
a <- oos.core.r2[!is.na(oos.core.r2)]
b <- oos.plot.r2[!is.na(oos.plot.r2)]
main.lab <- labs[1]
limx <- c(0,max(c(a, b))*1.1)
limy <- c(0, max(c(zero_truncated_density(a)$y,zero_truncated_density(b)$y))*1.02)
plot(cv.core.r2,ylim = limy, xlim = limx, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(zero_truncated_density(a), col = adjustcolor(cols[1],trans))
polygon(zero_truncated_density(b), col = adjustcolor(cols[2],trans))
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Validation R"^"2")), side = 1, line = 2.5, cex = o.cex)
mtext(main.lab,side = 3, line = -1.5, cex = o.cex)

##Plot legend
legend(x = 0.4, y = 25, legend = c('core','plot'), col ='black', pt.bg=adjustcolor(cols,trans), bty = 'n', pch = 22, pt.cex = 1.5, cex = 1.5)


#NEON CV R2----
a <- cv.core.r2[!is.na(cv.core.r2)]
b <- cv.plot.r2[!is.na(cv.plot.r2)]
limx <- c(0,max(c(a, b))*1.1)
limy <- c(0, max(c(zero_truncated_density(a)$y,zero_truncated_density(b)$y))*1.02)
main.lab <- labs[2]
plot(cv.core.r2,ylim = limy, xlim = limx, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(zero_truncated_density(a), col = adjustcolor(cols[1],trans))
polygon(zero_truncated_density(b), col = adjustcolor(cols[2],trans))
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Validation R"^"2")), side = 1, line = 2.5, cex = o.cex)
mtext(main.lab,side = 3, line = -1.5, cex = o.cex)

#Forecast OOS RMSE----
a <- oos.core.rmse
b <- oos.plot.rmse
limx <- c(0,max(c(a,b))*1.1)
limy <- c(0, max(c(zero_truncated_density(a)$y,zero_truncated_density(b)$y))*1.02)
main.lab <- labs[3]
plot(cv.core.rmse,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(zero_truncated_density(a), col = adjustcolor(cols[1],trans))
polygon(zero_truncated_density(b), col = adjustcolor(cols[2],trans))
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext('RMSE', side = 1, line = 2.5, cex = o.cex)
mtext(main.lab,side = 3, line = -1.5, cex = o.cex)


#NEON CV RMSE----
a <- cv.core.rmse
b <- cv.plot.rmse
limx <- c(0,max(c(a,b))*1.1)
limy <- c(0, max(c(zero_truncated_density(a)$y,zero_truncated_density(b)$y))*1.02)
main.lab <- labs[4]
plot(cv.core.rmse,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(zero_truncated_density(a), col = adjustcolor(cols[1],trans))
polygon(zero_truncated_density(b), col = adjustcolor(cols[2],trans))
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext('RMSE', side = 1, line = 2.5, cex = o.cex)
mtext(main.lab,side = 3, line = -1.5, cex = o.cex)

#end plot.----
dev.off()
