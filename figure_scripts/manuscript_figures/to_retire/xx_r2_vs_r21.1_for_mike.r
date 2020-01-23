#Comparing R2 and R2 1:1 for Mike.
rm(list=ls())
source('paths.r')
source('NEFI_functions/zero_truncated_density.r')

#set output path.----
output.path <- 'test.png'

#Load calibration data.----
d <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)

#get calibration/validation values by group.----
cal.mu.a <- d$calibration$cal.stat.predictable.sum$rsq
val.mu.a <- d$validation$site.stat.predictable.sim$rsq
cal.mu.b <- d$calibration$cal.stat.predictable.sum$rsq.1
val.mu.b <- d$validation$site.stat.predictable.sim$rsq.1
names(cal.mu.a) <- rownames(d$calibration$cal.stat.predictable.sum)
names(cal.mu.b) <- rownames(d$calibration$cal.stat.predictable.sum)

#png save line.----
png(filename=output.path,width=9,height=6,units='in',res=300)

#global plot settings.----
par(mfrow = c(1,2), mar = c(5,4,2,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.2 #outer label size.
cols <- c('purple','cyan','yellow')
 

#Calibration and Validation rsq ~ function/phylo scale.----
x <- 1:length(cal.mu.a)
limy <- c(0,max(cal.mu.a)*1.1)
plot(cal.mu.a ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu.a, lty = 2)
lines(x, val.mu.a, lty = 2)
points(val.mu.a ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level R"^"2")), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= names(cal.mu.a), srt=45, adj=1, xpd=TRUE, cex = 1.1)
#legend
legend(x = 1, y = 0.1, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2)
mtext('(a)', side = 3, adj = 0.95, line = -2)

#Calibration and Validation RMSE ~ function/phylo scale.----
x <- 1:length(cal.mu.b)
limy <- c(min(val.mu.b)*1.1,max(cal.mu.b)*8)
plot(cal.mu.b ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, cal.mu.b, lty = 2)
lines(x, val.mu.b, lty = 2)
points(val.mu.b ~ x, cex = 2.5, pch = 16, col = 'gray')
mtext(expression(paste("Site-Level R"^"2"," 1:1")), side = 2, line = 2.5, cex = o.cex)
axis(1, labels = F)
text(x=x+0.05, y = limy[1] - abs(limy[1])*0.1, labels= names(cal.mu.a), srt=45, adj=1, xpd=TRUE, cex = 1.1)
mtext('(b)', side = 3, adj = 0.95, line = -2)

#end plot.----
dev.off()
