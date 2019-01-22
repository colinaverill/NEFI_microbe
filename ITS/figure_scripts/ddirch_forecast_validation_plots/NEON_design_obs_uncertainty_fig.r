#density plots of predictability given observation uncertainty and sampling effort.
rm(list=ls())
source('NEFI_functions/zero_truncated_density.r')
source('paths.r')

#set output path.
output.path <- 'test.png'
output.path <- NEON_design_observation_uncertainty.path

#load data.----
d <- readRDS(ITS_NEON_observation_uncertainty.path)

#get means and 95% intervals.
all <- list()
for( i in 1:length(d)){
  ck <- d[[i]]
  mean <- apply(ck[,1:3], 2, mean)
  lo95 <- apply(ck[,1:3], 2, quantile, probs = c(0.025))
  hi95 <- apply(ck[,1:3], 2, quantile, probs = c(0.975))
  output <- list(mean,lo95,hi95)
  names(output) <- c('mean','lo95','hi95')
  all[[i]] <- output
}
names(all) <- names(d)
mu <- c(all$core$mean, all$plot$mean, all$site$mean)
lo95 <- c(all$core$lo95, all$plot$lo95, all$site$lo95)
hi95 <- c(all$core$hi95, all$plot$hi95, all$site$hi95)

#png save line.----
png(filename=output.path,width=7,height=7,units='in',res=300)

#Global plot settings.----
par(mfrow = c(1,1), 
    mar = c(0.4,1,1,1),
    oma = c(0.1,.1,0.1,6))
space <- 0.5
x <- c(c(0.9, 1, 1.1),c(0.9,1,1.1)+space,c(0.9,1,1.1)+space*2)
limx <- c(min(x) - .1,max(x) + 0.1) 
limy <- c(0,1)
o.cex <- 1.3 #outer label size.
cols <- c('purple','cyan','yellow')
par(mfrow = c(1,1), mar = c(4.5,4,1,1))
plot(mu ~ x,ylim = limy, xlim = limx, bty = 'n', xaxt='n', xlab = NA, ylab = NA, cex =0)
abline(h = 0.95, col = 'gray', lty = 2, lwd = 2)
arrows(x, lo95, x, hi95, length=0.05, angle=90, code=3)
points(x, mu, pch = 21, cex = 1.3, col = 'black', bg = rep(cols,3))
#axis(1, labels = F)
#abline(h = 0) #x bounding line.
#x axis labels
mtext('core-scale', side = 1, adj = 0.1, cex = o.cex)
mtext('plot-scale', side = 1, adj = 0.5, cex = o.cex)
mtext('site-scale', side = 1, adj = 0.9, cex = o.cex)
mtext('R2', side = 2, line = 2, cex= o.cex)
#legend
legend(x = 2.2, y = 0.9, 
       legend = c('high abundance','low abundance','low variation'), 
       col ='black', pt.bg=cols, 
       bty = 'n', pch = 21, pt.cex = 1.5, xpd=NA)

#end plot.----
dev.off()