#Density plots of variable importance.
rm(list=ls())
source('paths.r')

#set output path.----
#output.path <- 'test.png'
output.path <- ITS_variable.imp_all_fungi.path

#load variable importance data.----
fg <- readRDS(NEON_fg_variable_importance_data.path)
pl <- readRDS(NEON_phylo.group_variable_importance_data.path)

#Get a matrix where columns are variables, rows are fungal groups, entries are variable importance scores.----
all <- list()
for(i in 1:length(pl)){all[[i]] <- pl[[i]]$mean}
all <- do.call(cbind, all)
all <- cbind(fg$mean,all)
all <- all[,-grep('other',colnames(all))] #drop 'other' category.
all <- t(all)
#all <- log(all)
#get density.
d <- list()
for(i in 1:ncol(all)){
  to.dens <- all[,i]
  #to.dens <- to.dens[to.dens < 2]
  to.dens <- log10(to.dens)
  d[[i]] <- density(to.dens)
  }

#get predictor names
namey <- c('%C','C:N','pH','NPP','MAT','MAP','Forest (0-1)','Conifers (0-1)','% ECM Trees')

#get common x and y limits.----
xmin <- list()
xmax <- list()
ymax <- list()
for(i in 1:length(d)){
  xmin[[i]] <- min(d[[i]]$x)
  xmax[[i]] <- max(d[[i]]$x)
  ymax[[i]] <- max(d[[i]]$y)
}
limx <- c(min(unlist(xmin)), max(unlist(xmax)))
limy <- c(0, max(unlist(ymax)))


#begin plot.----
#png setting.
png(filename=output.path,width=8,height=8,units='in',res=300)

#global plot settings.
trans <- 0.3 #shading transparency.
o.cex <- 1.3 #outer label size.
shade_col <- 'green'

par(mfrow = c(3,3),
    mar = c(5,1,1,1),
    oma = c(4,2.75,1,1))
for(i in 1:length(d)){
  dens <- d[[i]]
  plot(dens, xlim = limx, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0, axes=F)
  axis(1, pos = 0)
  polygon(dens, col = adjustcolor(shade_col,trans))
  mtext(namey[i], side = 1, line = 2.8, cex = o.cex)
}
#outer labels.
mtext(expression(paste('log' ['10'],' effect size')), side = 1, outer = T, cex = 1.8, line = 2)
mtext(expression(paste('Relative Density')), side = 2, outer = T, cex = 1.8, line = -1.5)

#end plot.
dev.off()

#Single panel, overlapping.
#cols <- RColorBrewer::brewer.pal(9, 'Pastel1')
#par(mfrow=c(1,1))
#plot(dens, xlim = c(-2.5,0.5), ylim = limy, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0, axes=F)
#hist(log10(all[,i]), xlim = c(-2.5,0.5), bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1)
#axis(1, pos = 0)
#axis(2, at=limx, labels=c("",""), lwd.ticks=0)
#axis(2, las = 1)
#for(i in 1:length(d)){
#  polygon(d[[i]], col = adjustcolor(cols[i],trans))
#}
#mtext('effect size', side = 1, line = 2.5, cex = o.cex)


