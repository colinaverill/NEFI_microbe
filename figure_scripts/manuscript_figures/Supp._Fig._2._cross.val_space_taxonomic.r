#Comparing R2 density between Tedersoo out-of-sample (OOS) forecast to NEON and NEON cross-validation (CV).
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')
source('NEFI_functions/rsq_1.1.r')

#set figure out put path.----
output.path <- 'figures/Supp._Fig._2._cross.validation_r2.png'

#Load ITS cross-validation data..----
#Cross validation forecasts
cv.val.core.ITS <- readRDS(core.CV_NEON_fcast.path)
cv.val.plot.ITS <- readRDS(plot.CV_NEON_fcast.path)

#Cross validation data (calibration and validation)
cv.dat.core.ITS <- readRDS(core.CV_NEON_cal.val_data.path)
cv.dat.plot.ITS <- readRDS(plot.CV_NEON_cal.val_data.path)
cv.dat.core.ITS <- cv.dat.core.ITS$val$y.val
cv.dat.plot.ITS <- cv.dat.plot.ITS$val$y.val


#Load 16S cross-validation data.-----
#Cross validation forecasts
cv.val.core.16S <- readRDS(core.CV_NEON_fcast_16S.path)
cv.val.plot.16S <- readRDS(plot.CV_NEON_fcast_16S.path)

#Cross validation data (calibration and validation)
cv.dat.core.16S <- readRDS(core.CV_NEON_cal.val_data_16S.path)
cv.dat.plot.16S <- readRDS(plot.CV_NEON_cal.val_data_16S.path)
cv.dat.core.16S <- cv.dat.core.16S$val$y.val
cv.dat.plot.16S <- cv.dat.plot.16S$val$y.val

# Create table to link deprecatedVialID and geneticSampleID
map <- readRDS(core_obs.path)[,c("geneticSampleID","deprecatedVialID")]
map$geneticSampleID <- gsub('-GEN','',map$geneticSampleID)
#drop observations with no deprecated Vial ID.
map <- map[!(map$deprecatedVialID == ''),]

#16S NEON Cross-validated stats at core and plot level.-----
#core level.
core.list.16S <- list()
for(i in 1:length(cv.val.core.16S)){
  x <- cv.val.core.16S[[i]]$core.fit$mean
  y <- cv.dat.core.16S[[i]]$rel.abundances
  #updated y rownames.
  ref <- map[map$deprecatedVialID %in% rownames(y),]
    y <-   y[rownames(y) %in% ref$deprecatedVialID,]
  ref <- ref[order(match(ref$deprecatedVialID, rownames(y))),]
  rownames(y) <- ref$geneticSampleID

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
  core.list.16S[[i]] <- lev.stats
}
names(core.list.16S) <- names(cv.val.core.16S)

#plot level.
plot.list.16S <- list()
for(i in 1:length(cv.val.plot.16S)){
  x <- cv.val.plot.16S[[i]]$plot.fit$mean
  y <- cv.dat.plot.16S[[i]]$mean
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
  plot.list.16S[[i]] <- lev.stats
}
names(plot.list.16S) <- names(cv.val.plot.16S)

#16S collapse functional groups.----
phylo <- c('phylum','class','order','family','genus')
ref <- core.list.16S[!(names(core.list.16S) %in% phylo)]
core.list.16S <- core.list.16S[names(core.list.16S) %in% phylo]
core.list.16S$functional <- do.call(rbind, ref)
ref <- plot.list.16S[!(names(plot.list.16S) %in% phylo)]
plot.list.16S <- plot.list.16S[names(plot.list.16S) %in% phylo]
plot.list.16S$functional <- do.call(rbind, ref)

#16S means by taxaonomic level.----
#core scale.
core.lev.mu.16S <- list()
for(i in 1:length(core.list.16S)){
  z <- data.frame(core.list.16S[[i]])
  core.lev.mu.16S[[i]] <- mean(z$rsq1, na.rm = T)
}
core.lev.mu.16S <- unlist(core.lev.mu.16S)
names(core.lev.mu.16S) <- names(core.list.16S)
ref.order <- c(length(core.lev.mu.16S), 1:(length(core.lev.mu.16S) - 1))
core.lev.mu.16S <- core.lev.mu.16S[ref.order]

#plot scale.
plot.lev.mu.16S <- list()
for(i in 1:length(plot.list.16S)){
  z <- data.frame(plot.list.16S[[i]])
  plot.lev.mu.16S[[i]] <- mean(z$rsq1, na.rm = T)
}
plot.lev.mu.16S <- unlist(plot.lev.mu.16S)
names(plot.lev.mu.16S) <- names(plot.list.16S)
ref.order <- c(length(plot.lev.mu.16S), 1:(length(plot.lev.mu.16S) - 1))
plot.lev.mu.16S <- plot.lev.mu.16S[ref.order]

#put together.
to.plot <- cbind(core.lev.mu, plot.lev.mu)
to.plot.16S <- cbind(core.lev.mu.16S, plot.lev.mu.16S)

#ITS NEON Cross-validated stats at core and plot level.-----
#core level.
core.list.ITS <- list()
for(i in 1:length(cv.val.core.ITS)){
  x <- cv.val.core.ITS[[i]]$core.fit$mean
  y <- cv.dat.core.ITS[[i]]$rel.abundances
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
  core.list.ITS[[i]] <- lev.stats
}
names(core.list.ITS) <- names(cv.val.core.ITS)

#plot level.
plot.list.ITS <- list()
for(i in 1:length(cv.val.plot.ITS)){
  x <- cv.val.plot.ITS[[i]]$plot.fit$mean
  y <- cv.dat.plot.ITS[[i]]$mean
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
  plot.list.ITS[[i]] <- lev.stats
}
names(plot.list.ITS) <- names(cv.val.plot.ITS)


#ITS means by taxaonomic level.----
#core scale.
core.lev.mu.ITS <- list()
for(i in 1:length(core.list.ITS)){
  z <- data.frame(core.list.ITS[[i]])
  core.lev.mu.ITS[[i]] <- mean(z$rsq, na.rm = T)
}
core.lev.mu.ITS <- unlist(core.lev.mu.ITS)
names(core.lev.mu.ITS) <- names(core.list.ITS)
names(core.lev.mu.ITS)[length(core.lev.mu.ITS)] <- 'functional'
ref.order <- c(length(core.lev.mu.ITS), 1:(length(core.lev.mu.ITS) - 1))
core.lev.mu.ITS <- core.lev.mu.ITS[ref.order]

#plot scale.
plot.lev.mu.ITS <- list()
for(i in 1:length(plot.list.ITS)){
  z <- data.frame(plot.list.ITS[[i]])
  plot.lev.mu.ITS[[i]] <- mean(z$rsq, na.rm = T)
}
plot.lev.mu.ITS <- unlist(plot.lev.mu.ITS)
names(plot.lev.mu.ITS) <- names(plot.list.ITS)
names(plot.lev.mu.ITS)[length(plot.lev.mu.ITS)] <- 'functional'
ref.order <- c(length(plot.lev.mu.ITS), 1:(length(plot.lev.mu.ITS) - 1))
plot.lev.mu.ITS <- plot.lev.mu.ITS[ref.order]

#put together.
to.plot.ITS <- cbind(core.lev.mu.ITS, plot.lev.mu.ITS)

#png save line.-----
png(filename=output.path,width=8,height=5,units='in',res=300)

#Begin plot: Global plot settings.-----
par(mfrow = c(1,2),
    mar = c(4,4,2,2))
limx <- c(0,1)
limy <- c(0, 25)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
y.cex <- 1.3
cols <- brewer.pal(nrow(to.plot.16S) - 1,'Spectral')
cols <- c(cols,'green')

#16S by spatial and taxonomic scale.----
limy <- c(0, max(to.plot.16S)*1.1)
x <- c(1:2)
limx <- c(min(x)*0.95, max(x)*1.05)
plot(to.plot.16S[1,] ~ x, xlim = limx, ylim = limy, bty = 'l', col = 'black', bg = cols[1], 
     pch = 21, cex = 2.5, xlab= NA, ylab = NA, xaxt = 'n')
lines(x, to.plot.16S[1,], lty = 1, col = cols[1], lwd = 1.5)
for(i in 2:nrow(to.plot.16S)){
  points(to.plot.16S[i,] ~ x, pch = 21, cex = 2.5, col = 'black', bg = cols[i])
  lines(x, to.plot.16S[i,], lty = 1, lwd = 1.5,    col = cols[i])
}
axis(1, at = c(1,2,3), labels = F)
lab <- c('core','plot')
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= lab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
rsq.1.lab <- bquote({'Cross-Validated R'^2} [1:1])
mtext(rsq.1.lab, side = 2, line = 2.5, cex = y.cex)
#mtext(expression(paste("Cross-Validated R"^"2")), side = 2, line = 2.5, cex = y.cex)
leg.lab <- rownames(to.plot.16S)
legend('bottomright',leg.lab, pch = 21, col = 'black', pt.bg = cols, ncol = 2, bty= 'n')
mtext('(a)', side = 3, adj = 0.98, line = -1)
mtext('Bacteria', side = 3, adj = 0.03, line = -1, cex = o.cex)

#ITS by spatial and taxonomic scale.----
#limy <- c(0, max(to.plot.ITS)*1.1)
x <- c(1:2)
limx <- c(min(x)*0.95, max(x)*1.05)
plot(to.plot.ITS[1,] ~ x, xlim = limx, ylim = limy, bty = 'l', col = 'black', bg = cols[1], 
     pch = 21, cex = 2.5, xlab= NA, ylab = NA, xaxt = 'n')
lines(x, to.plot.ITS[1,], lty = 1, col = cols[1], lwd = 1.5)
for(i in 2:nrow(to.plot.ITS)){
  points(to.plot.ITS[i,] ~ x, pch = 21, cex = 2.5, col = 'black', bg = cols[i])
  lines(x, to.plot.ITS[i,], lty = 1, lwd = 1.5,    col = cols[i])
}
axis(1, at = c(1,2,3), labels = F)
lab <- c('core','plot')
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= lab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
rsq.1.lab <- bquote({'Cross-Validated R'^2} [1:1])
mtext(rsq.1.lab, side = 2, line = 2.5, cex = y.cex)
#mtext(expression(paste("Cross-Validated R"^"2")), side = 2, line = 2.5, cex = y.cex)
leg.lab <- rownames(to.plot.ITS)
#legend('topleft',leg.lab, pch = 21, col = 'black', pt.bg = cols, ncol = 1, bty= 'n')
mtext('(b)', side = 3, adj = 0.98, line = -1)
mtext('Fungi', side = 3, adj = 0.03, line = -1, cex = o.cex)

#end plot.-----
dev.off()
