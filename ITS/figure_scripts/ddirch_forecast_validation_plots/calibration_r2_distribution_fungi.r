#Plot in-sample calibration R2 values as a density plot.
rm(list=ls())
source('paths.r')

#Set output path.----
output.path <- 'test.png'
#output.path <- ITS_prior_r2_distribution.density_figure.path

#load data.----
pl <- readRDS(ted_ITS_prior_all.groups_JAGSfits.path) #all phylo and functional groups.
#re-order list, make function group first.
pl <- pl[c('fg','phylum','class','order','family','genus')]

#phylogenetic groups.
 pl.r2 <- list()
lev.r2.out <- list()
for(i in 1:length(pl)){
   lev <- pl[[i]]
   obs <- lev$observed
  pred <- lev$predicted
  lev.r2 <- list()
  for(j in 1:ncol(obs)){lev.r2[[j]] <- summary(lm(obs[,j] ~ pred[,j]))$r.squared}
  lev.r2 <- unlist(lev.r2)
  names(lev.r2) <- colnames(obs)
  lev.r2 <- lev.r2[names(lev.r2) != 'other']
  pl.r2[[i]] <- lev.r2
  names(pl.r2)[[i]] <- names(pl)[i]
}
lev.mu <- unlist(lapply(pl.r2, mean))
lev.sd <- unlist(lapply(pl.r2, sd))
lev.N  <- unlist(lapply(pl.r2, length))
lev.se <- lev.sd / sqrt(lev.N)
all.r2 <- unlist(pl.r2)
names(lev.mu)[1] <- 'functional'


#Get truncated density so that we don't have densities less than zero.----
h <- density(all.r2)$bw  #get dansity bandwith
# Compute edge weights.
w <- 1 / pnorm(0, mean=all.r2, sd=h, lower.tail=FALSE)
#Generate truncated distribution.
h <- density(all.r2, bw=h, kernel="gaussian", weights=w / length(all.r2))
h$y[h$all.r2 < 0] <- 0
#Check: the integral ought to be close to 1:
if(sum(h$y * diff(h$x)[1]) > 1.1 | sum(h$y * diff(h$x)[1]) < 0.9){
  warning('Zero-truncated density sum is not close to 1. Consider doing something else...')
}


#Make density plot.----
#png save settings.
png(filename=output.path,width=8,height=5,units='in',res=300)

#global plot settings.
trans <- 0.3 #shading transparency.
o.cex <- 1.3 #outer label size.
par(mfrow = c(1,2), mar = c(5,4.2,1.5,1.5))

#R2 density plot.
plot(h,xlim = c(0, 0.8), ylim = c(0, round(max(h$y), 0)), bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1)
polygon(h, col = adjustcolor('purple',trans), border = NA)
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Calibration R"^"2")), side = 1, line = 2.5, cex = o.cex)

#calbiration R2 as a function of phylo/function scale.
x <- 1:length(lev.mu)
limy <- c(0,max(lev.mu + lev.se))
plot(lev.mu ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='n', xaxt = 'n')
lines(x, lev.mu, lty = 2)
mtext(expression(paste("Calibration R"^"2")), side = 2, line = 2.2, cex = o.cex)
axis(1, labels = F)
text(x=x, y = -.04, labels= names(lev.mu), srt=45, adj=1, xpd=TRUE)

#end plot.
dev.off()