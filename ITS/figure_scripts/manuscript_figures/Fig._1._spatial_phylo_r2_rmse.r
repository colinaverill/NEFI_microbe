#plott calibration and validation R2 distributions for Fungi and bacteria as 6-panel plot.
rm(list=ls())
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')

#set output path.----
output.path <- 'test.png'

#Load calibration data.----
d <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)

#Subset to observations that have a minimum calibration R2 value.----
core.rsq <- d$validation$val.stat.predictable$core.stat
core.rsq <- data.frame(do.call(rbind, core.rsq))
plot.rsq <- d$validation$val.stat.predictable$plot.stat
plot.rsq <- data.frame(do.call(rbind, plot.rsq))
site.rsq <- d$validation$val.stat.predictable$site.stat
site.rsq <- data.frame(do.call(rbind, site.rsq))
core.d <- zero_truncated_density(core.rsq$rsq)
plot.d <- zero_truncated_density(plot.rsq$rsq)
site.d <- zero_truncated_density(site.rsq$rsq)
core.rmse.d <- zero_truncated_density(core.rsq$rmse.norm)
plot.rmse.d <- zero_truncated_density(plot.rsq$rmse.norm)
site.rmse.d <- zero_truncated_density(site.rsq$rmse.norm)
#get calibration values by group.
cal.mu <- d$calibration$cal.stat.predictable.sum$rsq
val.mu <- d$validation$site.stat.predictable.sim$rsq
cal.rmse <- d$calibration$cal.stat.predictable.sum$rmse.norm
val.rmse <- d$validation$site.stat.predictable.sim$rmse.norm
names(cal.mu) <- rownames(d$calibration$cal.stat.predictable.sum)

#png save line.----
png(filename=output.path,width=9,height=6,units='in',res=300)
#global plot settings.----
par(mfrow = c(2,4), mar = c(5,4,2,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.0 #outer label size.
cols <- c('purple','cyan','yellow')


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
mtext('(a)', side = 3, adj = 0.95, line = -2)

#Validation RMSE denisty plot.-----
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
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
legend(x = 1, y = 0.1, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(c)', side = 3, adj = 0.95, line = -2)

#Calibration and Validation RMSE ~ function/phylo scale.----
x <- 1:length(cal.rmse)
limy <- c(0,max(cal.rmse)*1.2)
plot(cal.rmse ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.rmse, lty = 2)
lines(x, val.rmse, lty = 2)
points(val.rmse ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level RMSE"["norm"])), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
#legend(x = 2, y = 0.6, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(d)', side = 3, adj = 0.95, line = -2)


#Outer labels.----
mtext('Fungi'   , cex = 1.5, side = 3, outer = T, adj = 0.01, line = -1)
mtext('Bacteria', cex = 1.5, side = 3, outer = T, adj = 0.01, line = -23)

#end plot.----
dev.off()
