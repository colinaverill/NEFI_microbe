#Plot in-sample calibration R2 values as a density plot.
rm(list=ls())
source('paths.r')

#Set output path.----
output.path <- 'test.png'
output.path <- ITS_prior_r2_distribution.density_figure.path

#load data.----
fg <- readRDS(ted_ITS.prior_fg_JAGSfit)
pl <- readRDS(ted_ITS_prior_phylo.group_JAGSfits)

#loop over lists and get R2 values.----
#functional groups.
obs <- fg$all.preds$observed
pred <- fg$all.preds$predicted
fg.r2 <- list()
for(i in 1:ncol(obs)){fg.r2[[i]] <- summary(lm(obs[,i] ~ pred[,i]))$r.squared}
fg.r2 <- unlist(fg.r2)
names(fg.r2) <- colnames(obs)
#drop 'other'.
fg.r2 <- fg.r2[names(fg.r2) != 'other']

#phylogenetic groups.
pl.r2 <- list()
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
pl.r2 <- unlist(pl.r2)

#Merge all r2 values.
all.r2 <- c(fg.r2, pl.r2)

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
png(filename=output.path,width=6,height=5,units='in',res=300)

#global plot settings.
trans <- 0.3 #shading transparency.
o.cex <- 1.3 #outer label size.
par(mfrow = c(1,1), mar = c(4.2,4.2,1.5,1.5))

#plot.
plot(h,xlim = c(0, 0.8), ylim = c(0, round(max(h$y), 0)), bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(h, col = adjustcolor('purple',trans), border = NA)
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Calibration R"^"2")), side = 1, line = 2.5, cex = o.cex)

#end plot.
dev.off()