#4 panel plots of calibration vs vaildation ~ functional/taxonomic scale, and moran's I across the NEON network.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
library(RColorBrewer)

#set output path.----
output.path <- 'Fig._1._phylo.moran_space.png'

#Load calibration/validation rsq.1 data, Moran's I data.----
d.ITS <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)
d.16S <- readRDS(NEON_dmilti.ddirch_analysis_summary_16S.path)
moran.ITS <- readRDS(NEON_ITS_morans_I_data.path)
moran.16S <- readRDS(NEON_16S_morans_I_data.path)

#1. Workup ITS calibration/validation rsq.1 values across phylo/function scales.----
#calibration rsq.1 values by functional/taxonomic group.
cal.check.ITS <- list()
cal.mu.ITS    <- list()
for(i in 1:length(d.ITS$calibration$cal.stat)){
  check <- d.ITS$calibration$cal.stat[[i]]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  cal.check.ITS[[i]] <- check
  cal.mu.ITS   [[i]] <- mean(check$rsq.1)
}
cal.mu.ITS <- unlist(cal.mu.ITS)
names(cal.mu.ITS) <- names(d.ITS$calibration$cal.stat)

#validation rsq.1 values by functional/taxonomic group.
val.check.ITS <- list()
val.mu.ITS    <- list()
for(i in 1:length(d.ITS$validation$val.stat$site.stat)){
  check <- d.ITS$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.ITS[[i]]$name,]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  val.check.ITS[[i]] <- check
  val.mu.ITS   [[i]] <- mean(check$rsq.1)
}
val.mu.ITS <- unlist(val.mu.ITS)
names(val.mu.ITS) <- names(d.ITS$validation$val.stat$site.stat)

#2. Workup 16S calibration/validation rsq.1 values across phylo/function scales.----
#get calibration rsq.1 values by group.
cal.check.16S <- list()
cal.mu.16S <- list()
for(i in 1:length(d.16S$calibration$cal.stat)){
  check <- d.16S$calibration$cal.stat[[i]]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  cal.check.16S[[i]] <- check
  cal.mu.16S   [[i]] <- mean(check$rsq.1)
}
cal.mu.16S <- unlist(cal.mu.16S)
names(cal.mu.16S) <- names(d.16S$calibration$cal.stat)

#get validation rsq.1 values by group.
val.check.16S <- list()
val.mu.16S    <- list()
for(i in 1:length(d.16S$validation$val.stat$site.stat)){
  check <- d.16S$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.16S[[i]]$name,]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  val.check.16S[[i]] <- check
  val.mu.16S   [[i]] <- mean(check$rsq.1)
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


#Global plot settings.----
png(filename = output.path, width = 8, height = 5, units = 'in', res = 300)
#par(mfrow = c(2,2), mar = c(5,4,2,1), oma = c(1,1,1,1))
m <- rbind(c(1,2,5), c(3,4,5))
layout(m)
layout.show(5)
par(mar = c(5,4,2,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.2 #outer label size.
y.cex <- 1.0 #y-axis label size.
x.cex <- 1.2 #x-axis label size.
cols <- c('purple','cyan','yellow')
#Calibration and Validation rsq ~ function/phylo scale 16S.----
x <- 1:length(cal.mu.16S)
limy <- c(-0.03,max(cal.mu.16S)*1.1)
#Calibration.
plot(cal.mu.16S ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu.16S, lty = 2)
#Validation.
lines(x, val.mu.16S, lty = 2)
points(val.mu.16S ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level R"^"2")), side = 2, line = 2.7, cex = y.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.16S), srt=45, adj=1, xpd=TRUE, cex = o.cex)
#legend
#legend(x = 1, y = 0.4, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(a)', side = 3, adj = 0.98, line = -2)
mtext('Bacteria', side = 3, line = 0.5, cex = 1.5, adj = 0)
#Bacterial Moran's I.----
limy <- c(0, max(moran.16S$mu + moran.16S$se))
plot(mu ~ x, data = moran.16S, cex = 2, pch = 16, ylim = limy,
     ylab = NA, xlab = NA, xaxt = 'n', bty = 'l', yaxs='i', las = 1)
#error bars.
mu <- moran.16S$mu
x <- moran.16S$x
upr <- mu + moran.16S$se
lwr <- mu - moran.16S$se
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)
#x-axis.
x.lab <- rownames(moran.16S)
x.lab[x.lab == 'fg'] <- 'functional'
axis(1, at=moran.16S$x, labels= NA, cex = 1, srt = 45)
text(x= moran.16S$x + .05, y = limy[1] - limy[2]*0.13, labels= x.lab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
mtext("Moran's I", side = 2, line = 2.5, cex = y.cex)
mtext('(b)', side = 3, adj = 0.98, line = -2)

#Calibration and Validation rsq ~ function/phylo scale ITS.----
x <- 1:length(cal.mu.ITS)
limy <- c(-0.03,max(cal.mu.ITS)*1.1)
plot(cal.mu.ITS ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu.ITS, lty = 2)
lines(x, val.mu.ITS, lty = 2)
points(val.mu.ITS ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level R"^"2")), side = 2, line = 2.7, cex = y.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.ITS), srt=45, adj=1, xpd=TRUE, cex = o.cex)
#legend
legend(x = 0.8, y = 0.06, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(c)', side = 3, adj = 0.98, line = -2)
mtext('Fungi', side = 3, line = 0.5, cex = 1.5, adj = 0)

#Fungal Moran's I.----
limy <- c(0, max(moran.ITS$mu + moran.ITS$se))
plot(mu ~ x, data = moran.ITS, cex = 2, pch = 16, ylim = limy,
     ylab = NA, xlab = NA, xaxt = 'n', bty = 'l', yaxs='i', las = 1)
#error bars.
mu <- moran.ITS$mu
x <- moran.ITS$x
upr <- mu + moran.ITS$se
lwr <- mu - moran.ITS$se
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)
#x-axis.
x.lab <- rownames(moran.ITS)
x.lab[x.lab == 'fg'] <- 'functional'
axis(1, at=moran.ITS$x, labels= NA, cex = 1, srt = 45)
#text(x= moran.ITS$x + .12, y = -0.06, labels= x.lab, srt=45, adj=1, xpd=TRUE, cex = 1)
text(x= moran.ITS$x + .05, y = limy[1] - limy[2]*0.13, labels= x.lab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
mtext("Moran's I", side = 2, line = 2.5, cex = y.cex)
mtext('(d)', side = 3, adj = 0.98, line = -2)

#Fungal validation rsq.1 across spatial scales.----
cols <- brewer.pal(nrow(scale.list) - 1,'Spectral')
cols <- c(cols,'green')
limy <- c(0, max(scale.list)*1.1)
x <- c(1:3)
limx <- c(min(x)*0.95, max(x)*1.05)
plot(scale.list[1,] ~ x, xlim = limx, ylim = limy, bty = 'l', col = 'black', bg = cols[1], 
     pch = 21, cex = 2.5, xlab= NA, ylab = NA, xaxt = 'n')
lines(x, scale.list[1,], lty = 1, col = cols[1], lwd = 1.5)
for(i in 2:nrow(scale.list)){
  points(scale.list[i,] ~ x, pch = 21, cex = 2.5, col = 'black', bg = cols[i])
  lines(x, scale.list[i,], lty = 1, lwd = 1.5,    col = cols[i])
}
axis(1, at = c(1,2,3), labels = F)
lab <- c('core','plot','site')
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= lab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
mtext(expression(paste("Fungal Validation R"^"2")), side = 2, line = 2.5, cex = y.cex)
leg.lab <- rownames(scale.list)
legend('topleft',leg.lab, pch = 21, col = 'black', pt.bg = cols, ncol = 1, bty= 'n')
mtext('(e)', side = 3, adj = 0.98, line = -2)

#end plot.----
dev.off()
