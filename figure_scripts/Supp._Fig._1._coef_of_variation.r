#4 panel plots of calibration vs vaildation ~ functional/taxonomic scale, and moran's I across the NEON network.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
library(RColorBrewer)

#set output path.----
output.path <- 'figures/Supp._Fig._1._CV.png'

#Load calibration/validation rsq.1 data, Moran's I data.----
d.ITS <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)
d.16S <- readRDS(NEON_dmilti.ddirch_analysis_summary_16S.path)

#1. Workup ITS calibration/validation coef. of variation values across phylo/function scales.----
#calibration rsq.1 values by functional/taxonomic group.
cal.check.ITS <- list()
cal.mu.ITS    <- list()
for(i in 1:length(d.ITS$calibration$cal.stat)){
  check <- d.ITS$calibration$cal.stat[[i]]
  cal.check.ITS[[i]] <- check
  cal.mu.ITS   [[i]] <- mean(check$rmse.norm)
}
cal.mu.ITS <- unlist(cal.mu.ITS)
names(cal.mu.ITS) <- names(d.ITS$calibration$cal.stat)

#validation rsq.1 values by functional/taxonomic group.
val.check.ITS <- list()
val.mu.ITS    <- list()
for(i in 1:length(d.ITS$validation$val.stat$site.stat)){
  check <- d.ITS$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.ITS[[i]]$name,]
  val.check.ITS[[i]] <- check
  val.mu.ITS   [[i]] <- mean(check$rmse.norm)
}
val.mu.ITS <- unlist(val.mu.ITS)
names(val.mu.ITS) <- names(d.ITS$validation$val.stat$site.stat)

#2. Workup 16S calibration/validation coef. of variation values across phylo/function scales.----
#get calibration CV values by group.
cal.check.16S <- list()
cal.mu.16S <- list()
for(i in 1:length(d.16S$calibration$cal.stat)){
  check <- d.16S$calibration$cal.stat[[i]]
  cal.check.16S[[i]] <- check
  cal.mu.16S   [[i]] <- mean(check$rmse.norm)
}
cal.mu.16S <- unlist(cal.mu.16S)
names(cal.mu.16S) <- names(d.16S$calibration$cal.stat)

#get validation rsq.1 values by group.
val.check.16S <- list()
val.mu.16S    <- list()
for(i in 1:length(d.16S$validation$val.stat$site.stat)){
  check <- d.16S$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.16S[[i]]$name,]
  val.check.16S[[i]] <- check
  val.mu.16S   [[i]] <- mean(check$rmse.norm)
}
val.mu.16S <- unlist(val.mu.16S)
names(val.mu.16S) <- names(d.16S$validation$val.stat$site.stat)

#3. Workup ITS validation rsq.1 values across function/phylo by spatial scales.----
scale.list <- list()
for(i in 1:length(d.ITS$validation$val.stat$core.stat)){
  core <- d.ITS$validation$val.stat$core.stat[[i]]$rsq.1
  plot <- d.ITS$validation$val.stat$plot.stat[[i]]$rsq.1
  site <- d.ITS$validation$val.stat$site.stat[[i]]$rsq.1
  core <- ifelse(core < 0, 0, core)
  plot <- ifelse(plot < 0, 0, plot)
  site <- ifelse(site < 0, 0, site)
  to.return <- c(mean(core),mean(plot),mean(site))
  scale.list[[i]] <- to.return
}
scale.list <- (do.call(rbind, scale.list))
rownames(scale.list) <- names(d.ITS$validation$val.stat$core.stat)
colnames(scale.list) <- c('core','plot','site')



#png save line.----
png(filename = output.path, width = 8, height = 5, units = 'in', res = 300)

#Global plot settings.----
par(mfrow = c(1,2), mar = c(5,4,2,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.2 #outer label size.
y.cex <- 1.0 #y-axis label size.
x.cex <- 1.0 #x-axis label size.
cols <- c('purple','cyan','yellow')
#Calibration and Validation rsq ~ function/phylo scale 16S.----
x <- 1:length(cal.mu.16S)
limy <- c(-0.03,max(val.mu.16S)*1.1)
#Calibration.
plot(cal.mu.16S ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu.16S, lty = 2)
#Validation.
lines(x, val.mu.16S, lty = 2)
points(val.mu.16S ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext('Coefficient of Variation', side = 2, line = 2.7, cex = y.cex)
axis(1, labels = F)
text(x=x+0.03, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.16S), srt=45, adj=1, xpd=TRUE, cex = x.cex)
#legend
legend(x = 1, y = 3.0, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.1)
mtext('(a)', side = 1, adj = 0.98, line = -1.5)
mtext('Bacteria', side = 3, line = 0.5, cex = 1.5, adj = 0)

#Calibration and Validation rsq ~ function/phylo scale ITS.----
x <- 1:length(cal.mu.ITS)
limy <- c(-0.03,max(cal.mu.ITS)*1.1)
plot(cal.mu.ITS ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu.ITS, lty = 2)
lines(x, val.mu.ITS, lty = 2)
points(val.mu.ITS ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext('Coefficient of Variation', side = 2, line = 2.7, cex = y.cex)
axis(1, labels = F)
text(x=x+0.03, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.ITS), srt=45, adj=1, xpd=TRUE, cex = x.cex)
#legend
legend(x = 0.8, y = 0.06, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(b)', side = 1, adj = 0.98, line = -1.5)
mtext('Fungi', side = 3, line = 0.5, cex = 1.5, adj = 0)

#end plot.----
dev.off()
