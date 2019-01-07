#compare out of sample forecast to number of SVs within a group.
rm(list=ls())
source('paths.r')
library(data.table)

#load data.----
#load in sample diversity stats.
tax <- readRDS(tedersoo_ITS_common_phylo_groups_list.path)

#calibration fits.
cal <- readRDS(ted_ITS_prior_phylo.group_JAGSfits)

#out of sample forecast to all phylo levels
fcast <- readRDS(NEON_site_fcast_all_phylo_levels.path)

#validation data for forecasts.
d <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq.path)
d <- readRDS('/fs/data3/caverill/phylo_cps.rds')

#load diversity stats for out of sample data.
div <- readRDS(NEON_ITS_fastq_all_cosmo_phylo_groups.path)

#Get calibration R2 values.----
cal.r2_table <- list()
for(i in 1:length(cal)){
  cal_level <- cal[[i]]
  pred <- cal_level$predicted
   obs <- cal_level$observed
  cal.lev_r2 <- list()
  for(j in 1:ncol(pred)){
    mod <- lm( obs[,j] ~ pred[,j])
    cal.lev_r2[[j]] <- summary(mod)$r.squared
  }
  cal.lev_r2 <- data.frame(unlist(cal.lev_r2))
  cal.lev_r2 <- cbind(colnames(obs), cal.lev_r2)
  colnames(cal.lev_r2) <- c('group','cal_rsq')
  cal.r2_table[[i]] <- cal.lev_r2
}
names(cal.r2_table) <- names(cal)

#get R2 tables within each phylo level, merge in diversity statistics, calibration r2.----
output1 <- list() #full table
for(i in 1:length(d)){
pred <- fcast[[i]]$site.fit$mean
 obs <- d[[i]]$site.fit$mean
pred <- pred[,order(match(colnames(pred), colnames(obs)))]
cal.rsq <- cal.r2_table[[i]]
 
#make sure same order and identity of rownames
 obs <-  obs[rownames(obs ) %in% rownames(pred),]
pred <- pred[rownames(pred) %in% rownames( obs),]
rsq <- list()
rsq.1 <- list()
for(j in 1:ncol(obs)){
  #get best-fit rsq.
  mod <- lm(obs[,j] ~ pred[,j])
  rsq[[j]] <- summary(mod)$r.squared
  #get rsq relative to 1:1 line.
  rss <- sum((pred[,j] - obs[,j]) ^ 2)       ## residual sum of squares
  tss <- sum((obs[,j] - mean(obs[,j])) ^ 2)  ## total sum of squares
  rsq.1[[j]] <- 1 - rss/tss
}
rsq <- data.frame(unlist(rsq))
rsq.1 <- data.frame(unlist(rsq.1))
rsq <- cbind(colnames(obs), rsq, rsq.1)
colnames(rsq) <- c('group','rsq','rsq.1')

rsq$group <- as.character(rsq$group)
div.merge <- div[[i]]$group_frequencies
rsq <- merge(rsq, div.merge, by.x = 'group', by.y = 'groups', all.x = T)
rsq <- merge(rsq, cal.rsq, all.x = T)
rsq <- rsq[!(names(rsq) %in% 'other')]
rsq$phylo_level <- names(d)[i]
output1[[i]] <- rsq
}
names(output1) <- names(d)
test <- do.call(rbind, output1)

#if rsq.1 is less than 0, set it to zero.
test$rsq.1 <- ifelse(test$rsq.1 < 0, 0, test$rsq.1)

#plot results.----
#out of sample rsq vs. calibration.
plot(rsq.1 ~ N.SVs, data = test[test$cal_rsq > 0.1 & test$rsq.1 > 0,])
mod <- lm(rsq.1 ~ N.SVs, data = test[test$cal_rsq > 0.1 & test$rsq.1 > 0,])
abline(mod)
#rsq ~ diversity for groups with in sample rsq > 0.15.
test$div.even <- test$diversity * test$evenness
test$div.even2 <- test$diversity * boot::logit(test$evenness)
test$nsv.even <- test$N.SVs * test$evenness
plot(rsq ~ N.SVs, data = test[test$cal_rsq,], col = as.factor(test$phylo_level))

#aggregate phylo-level rsq, where calibration greater than a threshold.
pos <- data.frame(c(1:5), c('genus','family','order','class','phylum'))
colnames(pos) <- c('pos','phylo_level')
pos$phylo_level <- as.character(pos$phylo_level)
 mu <- aggregate(rsq.1 ~ phylo_level, data = test[test$cal_rsq > 0.1,], FUN = mean  )
 cr <- aggregate(cal_rsq ~ phylo_level, data = test[test$cal_rsq > 0.1,], FUN = mean)
 sd <- aggregate(rsq.1 ~ phylo_level, data = test[test$cal_rsq > 0.1,], FUN = sd    )
  N <- aggregate(rsq.1 ~ phylo_level, data = test[test$cal_rsq > 0.1,], FUN = length)
mu <- mu[order(match(mu$phylo_level,pos$phylo_level)),]
sd <- sd[order(match(sd$phylo_level,pos$phylo_level)),]
 N <-  N[order(match( N$phylo_level,pos$phylo_level)),]
cr <- cr[order(match(cr$phylo_level,pos$phylo_level)),]
to_plot <- cbind(pos, mu[,2],sd[,2],N[,2],cr[,2])
colnames(to_plot)[3:6] <- c('mu','sd','N','cal_rsq')
to_plot$std.err <- to_plot$sd / sqrt(to_plot$N)

#make a super simple plot for AGU.
out.path <- 'rsq_figure.png'
png(filename=out.path,width=8,height=4,units='in',res=300)
par(mfrow = c(1,2),
    mai = c(1,.5,.3,.1),
    oma = c(0.5,3,1,1))
#calibration rsq values.
plot(cal_rsq ~ pos, data = to_plot, cex = 3, pch = 16, xlab = NA, xaxt = 'n', 
     ylim = c(0, 0.4), xlim = c(0.5,5.5), bty = 'n', ylab = NA)
lines(to_plot$pos, to_plot$cal_rsq, lty = 2)   
text(x = c(1:5), y = -0.05, srt = 45, adj = 1, labels = to_plot$phylo_level, xpd = TRUE, cex = 1.25)  
#mtext('Calibration r-squared',line = 0, side = 3, cex = 2)

#out of sample rsq values.
plot(mu ~ pos, data = to_plot, cex = 3, pch = 16, xlab = NA, xaxt = 'n', 
     ylim = c(0, 0.4), xlim = c(0.5,5.5), bty = 'n', ylab = NA)
lines(to_plot$pos, to_plot$mu, lty = 2)   
text(x = c(1:5), y = -0.05, srt = 45, adj = 1, labels = to_plot$phylo_level, xpd = TRUE, cex = 1.25)  
#mtext('Out-of-sample r-squared',line = 0, side = 3, cex = 2)

dev.off()