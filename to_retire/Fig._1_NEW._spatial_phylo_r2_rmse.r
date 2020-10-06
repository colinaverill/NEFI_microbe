#plott calibration and validation R2 distributions for Fungi and bacteria as 6-panel plot.
rm(list=ls())
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')

#set output path.----
output.path <- 'test.png'

#Load calibration data.----
d.ITS <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)
d.16S <- readRDS(NEON_dmilti.ddirch_analysis_summary_16S.path)

#Get ITS rsq.1 vectors within spatial scales and means across functional and taxonomic scales.----
#validation rsq.1 across spatial scales.
core.rsq.ITS <- d.ITS$validation$val.stat$core.stat
core.rsq.ITS <- data.frame(do.call(rbind, core.rsq.ITS))
core.rsq.ITS$rsq.1 <- ifelse(core.rsq.ITS$rsq.1 < 0, 0, core.rsq.ITS$rsq.1)
plot.rsq.ITS <- d.ITS$validation$val.stat$plot.stat
plot.rsq.ITS <- data.frame(do.call(rbind, plot.rsq.ITS))
plot.rsq.ITS$rsq.1 <- ifelse(plot.rsq.ITS$rsq.1 < 0, 0, plot.rsq.ITS$rsq.1)
site.rsq.ITS <- d.ITS$validation$val.stat$site.stat
site.rsq.ITS <- data.frame(do.call(rbind, site.rsq.ITS))
site.rsq.ITS$rsq.1 <- ifelse(site.rsq.ITS$rsq.1 < 0, 0, site.rsq.ITS$rsq.1)
scale.ITS <- list(core.rsq.ITS, plot.rsq.ITS, site.rsq.ITS)

#calibration rsq.1 values by functional/taxonomic group.
cal.check.ITS <- list()
cal.mu.ITS    <- list()
for(i in 1:length(d.ITS$calibration$cal.stat)){
  check <- d.ITS$calibration$cal.stat[[i]]
  #check <- check[check$rsq.1 > 0.0,]
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

#for 16S.----
#validation rsq.1 across spatial scales.
core.rsq.16S <- d.16S$validation$val.stat$core.stat
core.rsq.16S <- data.frame(do.call(rbind, core.rsq.16S))
core.rsq.16S$rsq.1 <- ifelse(core.rsq.16S$rsq.1 < 0, 0, core.rsq.16S$rsq.1)
plot.rsq.16S <- d.16S$validation$val.stat$plot.stat
plot.rsq.16S <- data.frame(do.call(rbind, plot.rsq.16S))
plot.rsq.16S$rsq.1 <- ifelse(plot.rsq.16S$rsq.1 < 0, 0, plot.rsq.16S$rsq.1)
site.rsq.16S <- d.16S$validation$val.stat$site.stat
site.rsq.16S <- data.frame(do.call(rbind, site.rsq.16S))
site.rsq.16S$rsq.1 <- ifelse(site.rsq.16S$rsq.1 < 0, 0, site.rsq.16S$rsq.1)
scale.16S <- list(core.rsq.16S, plot.rsq.16S, site.rsq.16S)
sum <- list()
for(i in 1:length(scale.16S)){
  z <- scale.16S[[i]]
  z <- z[z]
  sum[[i]] <- mean(scale.16S[[i]]$rsq.1)
}
unlist(sum)

#get calibration rsq.1 values by group.
cal.check.16S <- list()
   cal.mu.16S <- list()
for(i in 1:length(d.16S$calibration$cal.stat)){
  check <- d.16S$calibration$cal.stat[[i]]
  #check <- check[check$rsq.1 > 0.0,]
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

#png save line.----
#png(filename=output.path,width=9,height=6,units='in',res=300)
#global plot settings.----
par(mfrow = c(2,4), mar = c(5,4,2,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.0 #outer label size.
cols <- c('purple','cyan','yellow')

#Validation R2 denisty plot 16S.-----
dat.a <- core.d.2
dat.b <- plot.d.2
dat.c <- site.d.2
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
mtext('(a)', side = 3, adj = 0.95, line = -2)

#Validation RMSE denisty plot 16S.-----
dat.a <- core.rmse.d.2
dat.b <- plot.rmse.d.2
dat.c <- site.rmse.d.2
limy <- c(0, max(c(dat.a$y, dat.b$y, dat.c$y))*1.05)
limx <- c(0, max(c(dat.a$x, dat.b$x, dat.c$x)))
limx <- c(0, 10)
#Density plot.
plot(dat.a,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(dat.a, col = adjustcolor(cols[1],trans))
polygon(dat.b, col = adjustcolor(cols[2],trans))
polygon(dat.c, col = adjustcolor(cols[3],trans))
mtext('Density', side = 2, line = 2.5, cex = o.cex)
mtext(expression(paste("Validation RMSE"["norm"])), side = 1, line = 3, cex = o.cex)
#legend
legend(x = 0.5, y = 12.5, legend = c('core','plot','site'), col ='black', pt.bg=adjustcolor(cols,trans), bty = 'n', pch = 22, pt.cex = 1.5, cex = 1.2)
mtext('(b)', side = 3, adj = 0.95, line = -2)


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
mtext(expression(paste("Site-Level R"^"2")), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.16S), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
#legend(x = 1, y = 0.3, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(c)', side = 3, adj = 0.95, line = -2)

#Calibration and Validation RMSE ~ function/phylo scale 16S.----
x <- 1:length(cal.rmse.16S)
limy <- c(0,max(c(val.rmse.16S,cal.rmse.16S))*1.2)
plot(cal.rmse.16S ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
lines(x, cal.rmse.16S, lty = 2)
lines(x, val.rmse.16S, lty = 2)
points(val.rmse.16S ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level RMSE"["norm"])), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
#legend(x = 1, y = 3, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(d)', side = 3, adj = 0.95, line = -2)


#Validation R2 denisty plot ITS.-----
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
mtext('(e)', side = 3, adj = 0.95, line = -2)

#Validation RMSE denisty plot ITS.-----
dat.a <- core.rmse.d
dat.b <- plot.rmse.d
dat.c <- site.rmse.d
limy <- c(0, max(c(dat.a$y, dat.b$y, dat.c$y))*1.05)
limx <- c(0, max(c(dat.a$x, dat.b$x, dat.c$x)))
limx <- c(0, 10)
#Density plot.
plot(dat.a,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(dat.a, col = adjustcolor(cols[1],trans))
polygon(dat.b, col = adjustcolor(cols[2],trans))
polygon(dat.c, col = adjustcolor(cols[3],trans))
mtext('Density', side = 2, line = 2.5, cex = o.cex)
mtext(expression(paste("Validation RMSE"["norm"])), side = 1, line = 3, cex = o.cex)
#legend
legend(x = 0.5, y = 12.5, legend = c('core','plot','site'), col ='black', pt.bg=adjustcolor(cols,trans), bty = 'n', pch = 22, pt.cex = 1.5, cex = 1.2)
mtext('(f)', side = 3, adj = 0.95, line = -2)


#Calibration and Validation rsq ~ function/phylo scale ITS.----
x <- 1:length(cal.mu.ITS)
limy <- c(-0.03,max(cal.mu.ITS)*1.1)
plot(cal.mu.ITS ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu.ITS, lty = 2)
lines(x, val.mu.ITS, lty = 2)
points(val.mu.ITS ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level R"^"2")), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.ITS), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
#legend(x = 1, y = 0.1, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(g)', side = 3, adj = 0.95, line = -2)

#Calibration and Validation RMSE ~ function/phylo scale ITS.----
x <- 1:length(cal.rmse.ITS)
limy <- c(0,max(cal.rmse.ITS)*1.2)
plot(cal.rmse.ITS ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.rmse.ITS, lty = 2)
lines(x, val.rmse.ITS, lty = 2)
points(val.rmse.ITS ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level RMSE"["norm"])), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.ITS), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
#legend(x = 2, y = 0.6, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(h)', side = 3, adj = 0.95, line = -2)


#Outer labels.----
mtext('Bacteria', cex = 1.5, side = 3, outer = T, adj = 0.01, line = -1 )
mtext('Fungi'   , cex = 1.5, side = 3, outer = T, adj = 0.01, line = -23)

#end plot.----
#dev.off()
