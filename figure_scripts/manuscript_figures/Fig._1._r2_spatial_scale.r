#4 panel plots of calibration vs vaildation ~ functional/taxonomic scale, and moran's I across the NEON network.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
library(RColorBrewer)

#set output path.----
output.path <- 'figures/Fig._1._r2_spatial_scale.png'

#load data, work up means and standard errors.----
#requires the supplement file.
d <- read.csv('figures/Supplementary_Data_File_1.csv')


b.core <- mean(d[d$kingdom == 'Bacteria',]$rsq_1.1_core)
b.plot <- mean(d[d$kingdom == 'Bacteria',]$rsq_1.1_plot)
b.site <- mean(d[d$kingdom == 'Bacteria',]$rsq_1.1_site)
f.core <- mean(d[d$kingdom == 'Fungi'   ,]$rsq_1.1_core)
f.plot <- mean(d[d$kingdom == 'Fungi'   ,]$rsq_1.1_plot)
f.site <- mean(d[d$kingdom == 'Fungi'   ,]$rsq_1.1_site)

#standard errors.
b.core.se <- sd(d[d$kingdom == 'Bacteria',]$rsq_1.1_core) / sqrt(nrow(d[d$kingdom == 'Bacteria',]))
b.plot.se <- sd(d[d$kingdom == 'Bacteria',]$rsq_1.1_plot) / sqrt(nrow(d[d$kingdom == 'Bacteria',]))
b.site.se <- sd(d[d$kingdom == 'Bacteria',]$rsq_1.1_site) / sqrt(nrow(d[d$kingdom == 'Bacteria',]))
f.core.se <- sd(d[d$kingdom == 'Fungi'   ,]$rsq_1.1_core) / sqrt(nrow(d[d$kingdom == 'Bacteria',]))
f.plot.se <- sd(d[d$kingdom == 'Fungi'   ,]$rsq_1.1_plot) / sqrt(nrow(d[d$kingdom == 'Bacteria',]))
f.site.se <- sd(d[d$kingdom == 'Fungi'   ,]$rsq_1.1_site) / sqrt(nrow(d[d$kingdom == 'Bacteria',]))

#put together.
fun <- c(f.core, f.plot, f.site)
bac <- c(b.core, b.plot, b.site)
fun.se <- c(f.core.se, f.plot.se, f.site.se)
bac.se <- c(b.core.se, b.plot.se, b.site.se)

#png save line.----
png(filename = output.path, width = 8, height = 4, units = 'in', res = 300)

#Global plot settings.----
par(mfrow = c(1,2), mar = c(3,3,1,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.2 #outer label size.
y.cex <- 1.0 #y-axis label size.
x.cex <- 1.2 #x-axis label size.

#Fungal validation R2 across scales.----
x <- 1:length(fun)
xlab <- c('core','plot','site')
#error bars.
upr <- fun + fun.se
lwr <- fun - fun.se
limy <- c(0, max(upr)*1.1)
#points.
plot(fun ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
lines(x, fun, lty = 2)
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)

#labels
mtext(expression(paste("Validation R"^"2")[1:1]), side = 2, line = 2.7, cex = y.cex)
mtext('Fungi', side = 3, line = -1, adj = 0.03, cex = o.cex)
axis(1, at=c(1,2,3), labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= xlab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
#legend
mtext('(a)', side = 1, adj = 0.98, line = -2)

#Bacterial validation R2 across scales.----
x <- 1:length(bac)
xlab <- c('core','plot','site')
#error bars.
upr <- bac + bac.se
lwr <- bac - bac.se
limy <- c(0, max(upr)*1.1)
#Calibration.
plot(bac ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='l', xaxt = 'n', yaxs='i', las = 1, lwd = 0)
lines(x, bac, lty = 2)
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)

mtext(expression(paste("Validation R"^"2")[1:1]), side = 2, line = 2.7, cex = y.cex)
mtext('Bacteria', side = 3, line = -1, adj = 0.03, cex = o.cex)
axis(1, at=c(1,2,3), labels = F)
text(x=x+0.05, y = limy[1] - limy[2]*0.1, labels= xlab, srt=45, adj=1, xpd=TRUE, cex = o.cex)
#legend
mtext('(b)', side = 1, adj = 0.98, line = -2)


#end plot.----
dev.off()
