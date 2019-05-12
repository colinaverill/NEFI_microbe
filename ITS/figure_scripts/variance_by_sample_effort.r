#plotting effects of sampling effort NEON
rm(list=ls())
source('paths.r')
d <- readRDS(HARV_sampling_effort_analysis.path)

#get aggregated data.----
out <- list()
for(i in 1:length(d)){
  out[[i]] <- aggregate(mu ~ n.samp, data = d[[i]], FUN = mean)
}
names(out) <- names(d)
#drop ORNL, only 3 observations.
out$ORNL <- NULL

#assign to ecosystem types.
eco <- c('forest','forest','savannah','forest','forest','savannah','forest','grassland','forest','forest','grassland')
col <- c(rep('purple', length(eco)))
col <- ifelse(eco == 'savannah','orange',col)
col <- ifelse(eco == 'grassland','green',col)


#Get negative exponential fits + predicted values.----
fit <- list()
for(i in 1:length(out)){
  dat <- out[[i]]
  mod <- nls(mu ~ a*exp(b*n.samp) + c, data = dat, start = list(a = 32, b = -0.6, c = 3))
  range <- data.frame(seq(1, max(dat$n.samp), by = 0.1))
  colnames(range) <- c('n.samp')
  pred <- predict(mod, newdata = range)
  return <- list(mod, range$n.samp, pred)
  names(return) <- c('model','effort','pred')
  fit[[i]] <- return
}
names(fit) <- names(out)

#png save line.----
png('test.png', width = 5, height = 5, units = 'in', res = 300)

#plot them.----
#global plot settings.
par(mar = c(1,1,1,1), oma = c(3,3,1,1))
#get y-limit.
var.max <- list()
for(i in 1:length(out)){var.max[[i]] <- out[[i]]$mu}
var.max <- unlist(var.max)
limy <- c(0, max(var.max))
limx <- c(1, 30)

#plotting.
plot(mu ~ n.samp, data = dat, ylim = limy, xlim = limx, cex = 0)
#draw fits
for(i in 1:length(fit)){
  lines(smooth.spline(fit[[i]]$pred ~ fit[[i]]$effort), lty = 2, lwd = 1.5, col = col[i])
}
#draw legend.
legend(x = 17.5, y = 100,
       legend = c('forest','savannah','grassland'), lty = 2, lwd = 1.5, 
          col = c('purple','orange','green'),
       bty = 'n', y.intersp = 1)
#labels
mtext('Variance', side = 2, line = 2.2, cex = 1.2)
mtext('Number of soil cores', side = 1, line = 2.2, cex = 1.2)

#end plot.----
dev.off()
