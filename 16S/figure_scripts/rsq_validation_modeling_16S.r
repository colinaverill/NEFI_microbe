#compare out of sample forecast to number of SVs within a group, among other things.
rm(list=ls())
source('paths.r')
library(data.table)
library(MASS)
library(corrplot) 
library(caret)

#load in-/out-of-sample abundances and diversity stats.
cal.div <- readRDS(bahram_16S_common_phylo_fg_abun.path)
val.div <- readRDS(NEON_16S_phylo_fg_abundances.path)

# Load in calibration/validation models.
# cal <- readRDS(bahram_16S_prior_dmulti.ddirch_all.group_JAGSfits)
cal <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/prior_phylo_fg_JAGSfit_16S.rds"))
fg.cal <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/bahram_16S_prior_ddirch_fg_JAGSfits.rds"))
cal <- c(cal[1:5], fg.cal) # until the model is re-run, fg model must be loaded separately

# Load in forecast
val <- readRDS(NEON_cps_fcast_ddirch_16S.path)
val <- readRDS(paste0(pecan_gen_16S_dir, "/NEON_forecast_data/NEON_cps_fcast_ddirch_noLogMap.rds"))
# Load in validation data.
val.obs <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)


#Get calibration R2 values, merge in diversity stats.----
cal.r2_table <- list()
for(i in 1:length(cal)){
  cal_level <- cal[[i]]
  div <- cal.div[[i]]$group_frequencies
  pred <- cal_level$predicted
  obs <- cal_level$observed
  cal.lev_r2   <- list()
  cal.lev_r2.1 <- list()
  for(j in 1:ncol(pred)){
    mod <- lm(( obs[,j]) ~ (pred[,j]))
    cal.lev_r2[[j]] <- summary(mod)$r.squared
    #relative to 1:1 line.
    rss <- sum((pred[,j] - obs[,j]) ^ 2)       ## residual sum of squares
    tss <- sum((obs[,j] - mean(obs[,j])) ^ 2)  ## total sum of squares
    rsq.1 <- 1 - rss/tss
    rsq.1 <- ifelse(rsq.1 < 0, 0, rsq.1)
    cal.lev_r2.1[[j]] <- rsq.1
  }
  cal.lev_r2 <- data.frame(unlist(cal.lev_r2))
  cal.lev_r2 <- cbind(colnames(obs), cal.lev_r2)
  colnames(cal.lev_r2) <- c('group','cal_rsq')
  cal.lev_r2$cal_rsq.1 <- unlist(cal.lev_r2.1)
  #merge in diversity stats.
  colnames(div) <- c('group', 'cal_samp_freq','cal_N.SV','cal_diversity','cal_evenness')
  cal.lev_r2 <- merge(cal.lev_r2, div, all.x = T)
  cal.lev_r2 <- cal.lev_r2[-grep('other',cal.lev_r2$group),]
  cal.r2_table[[i]] <- cal.lev_r2
}
names(cal.r2_table) <- names(cal)

#get validation R2 values (site-level), merge in diversity stats.----
val.r2_table <- list() #full table
for(i in 1:length(val)){
  pred <-     val[[i]]$site.fit$mean
  obs <- val.obs[[i]]$site.fit$mean
  pred <- pred[,order(match(colnames(pred), colnames(obs)))]
  
  #make sure same order and identity of rownames
  obs <-  obs[rownames(obs ) %in% rownames(pred),]
  pred <- pred[rownames(pred) %in% rownames( obs),]
  rsq   <- list()
  rsq.1 <- list()
  abundance <- list()
  variance  <- list()
  range     <- list()
  for(j in 1:ncol(obs)){
    #get best-fit rsq.
    mod <- lm(obs[,j] ~ pred[,j])
    rsq[[j]] <- summary(mod)$r.squared
    #get rsq relative to 1:1 line.
   # rss <- sum((pred[,j] - obs[,j]) ^ 2)       ## residual sum of squares (-colin)
    rss <- sum(resid(mod)^2)      ## residual sum of squares
    tss <- sum((obs[,j] - mean(obs[,j])) ^ 2)  ## total sum of squares
    rsq_1 <- 1 - (rss/tss)
    rsq_1 <- ifelse(rsq_1 < 0, 0, rsq_1)
    rsq.1[[j]] <- rsq_1
    abundance[[j]] <- boot::inv.logit(mean(boot::logit(obs[,j])))
    variance[[j]] <- boot::inv.logit(  sd(boot::logit(obs[,j])))
    range[[j]] <- max(obs[,j]) - min(obs[,j])
  }
  rsq <- data.frame(unlist(rsq))
  rsq.1 <- data.frame(unlist(rsq.1))
  abundance <- data.frame(unlist(abundance))
  variance <- data.frame(unlist(variance))
  range <- data.frame(unlist(range))
  rsq <- cbind(colnames(obs), rsq, rsq.1,abundance, variance,range)
  colnames(rsq) <- c('group','val_rsq','val_rsq.1','val_abundance','val_variance','val_range')
  
  rsq$group <- as.character(rsq$group)
  div <- val.div[[i]]$group_frequencies
  colnames(div) <- c('group','val_samp_freq','val_N.SV','val_diversity','val_evenness')
  rsq <- rsq[-grep('other',rsq$group),]
  rsq <- merge(rsq, div, all.x = T)
  rsq$phylo_level <- names(val)[i]
  val.r2_table[[i]] <- rsq
}
names(val.r2_table) <- names(val)

#merge the tables.----
cal.stats <- do.call(rbind, cal.r2_table)
val.stats <- do.call(rbind, val.r2_table)
d <- merge(cal.stats, val.stats)

# rename functional groups
d[!d$phylo_level %in% c('phylum','class','order','family','genus'),]$phylo_level <- "fg"

#Model the R2 values.----

# Fungi results:
#vlaidation rsq of taxa with cal_rsq > 0.10.
#calibration and validation rsq are not related.
#validation rsq is related to validation variance, N.SVs (which are correlated). Best model just includes val_variance.
#diversity, N.SVs, variance all correlated, variance best predictor.
#also related to phylo level, but again, this mkaes val_variance model worse.
#true whether we model rsq or rsq.1 
#--------

#### now for bacteria ####
mod <- lm(val_rsq ~ cal_rsq + phylo_level,data = d[d$cal_rsq > 0.1,])
summary(mod)


mod <- lm(val_rsq ~ cal_rsq + phylo_level,data = d[d$cal_rsq > 0.33,])
summary(mod)
par(mfrow = c(2,2)) # no clear patterns, other than with variance
plot(val_rsq ~ cal_rsq        , data = d[d$cal_rsq > 0.33,])
plot(val_rsq ~ log10(val_N.SV), data = d[d$cal_rsq > 0.1,])
plot(val_rsq ~ val_diversity  ,data = d[d$cal_rsq > 0.1,])
plot(val_rsq ~ val_variance   , data = d[d$cal_rsq > 0.1,])


# subset for stepwise regression
h <- d[d$cal_rsq > 0.33,]
x <- h[,!colnames(h) %in% c("val_rsq", "val_rsq.1", "group")]
fit <- lm(h$val_rsq ~ .,data=x)
step <- stepAIC(fit, direction="both")
summary(step)
step$anova # display results

# best model includes 7 predictors (so, most of them) and has an adjusted R2 of .5296
# doesn't account for collinearity:
summary(lm(d$val_rsq ~ cal_rsq.1 + val_rsq.1 + val_abundance + val_variance + 
             val_range + val_samp_freq + val_N.SV + val_diversity + val_evenness, data=d)) 
# covariate correlation plot
dm = data.matrix(d[,!colnames(d) %in% c("val_rsq", "val_rsq1","phylo_level", "group")]) 
corrplot(cor(dm), order = 'hclust') 


# can we see why phylum forecasts split neatly into good and bad???
# no. n=11 is too small. val_sample_freq is significant, except it's the phyla present in ~100% of samples that have 0% predictability...
p <- d[d$phylo_level=="phylum",]
t <- p[,!colnames(p) %in% c("val_rsq", "phylo_level", "group")]
summary(lm(val_rsq ~ val_samp_freq, data = p)) #p=.46
fit <- lm(p$val_rsq~.,data=t)
step <- stepAIC(fit, direction="both")
summary(step)
vif(fit)
sqrt(vif(fit)) > 2
plot(val_rsq ~ val_samp_freq, data = p, ylab = 'validation R2', xlab = 'val sample freq R2') # too few points

# all tax groups except phyla 
np <- d[d$phylo_level != "phylum",]
summary(lm(val_rsq ~ val_samp_freq, data = np[np$cal_rsq > 0.1,])) 

# ok, do functional groups show any trends?
fg <- d[d$phylo_level=="fg",]
mod <- lm(val_rsq ~ cal_rsq,data = fg[fg$cal_rsq > 0.33,])
summary(mod)
fit <- lm(fg$val_rsq ~ . ,data=fg[,-c(4,8,9,13,17)])
step <- stepAIC(fit, direction="both")
summary(step)
step$anova # display results

#When you subset to calibration R2 > 0.33, calibration does predict validation, R2 = 0.49, unadjusted. <- not true for bacteria at all :(
par(mfrow=(c(1,2)))
mod1 <- lm(val_rsq ~ cal_rsq, data = d[d$cal_rsq > 0.33,])
plot(val_rsq ~ cal_rsq, data = d[d$cal_rsq > 0.33,], ylab = 'validation R2', xlab = 'calibration R2')
abline(mod1, lwd = 2)
mtext(paste0('R2 = ',round(summary(mod1)$r.squared,2)), side = 3, line = -2, adj = 0.05)

#You can explain 88% (R2 = 0.88, unadjusted) if you add phylo/functional group to the model. <- also not true for bacteria
mod2 <- lm(val_rsq ~ cal_rsq + phylo_level, data = d[d$cal_rsq > 0.33,])
y <- d[d$cal_rsq > 0.33,]$val_rsq
x <- fitted(mod2)
plot(y ~ x, ylab = 'observed validation R2', xlab = 'predicted validation R2')
abline(lm(y~x), lwd = 2)
mtext(paste0('R2 = ',round(summary(mod2)$r.squared,2)), side = 3, line = -2, adj = 0.05)



# You can explain .59 (.52 adj) if you add a lot to the model
mod2 <- lm(val_rsq ~ cal_rsq.1 + cal_evenness + val_abundance + val_variance + 
             val_range + val_samp_freq + phylo_level, data = d[d$cal_rsq > 0.33,])
y <- d[d$cal_rsq > 0.33,]$val_rsq
x <- fitted(mod2)
plot(y ~ x, ylab = NA, xlab = NA, bty =  'l', pch = 16, cex = 1.5, main = "val_rsq ~ cal_rsq.1 + cal_evenness \n+ val_abundance + val_variance + 
             val_range + val_samp_freq + phylo_level, \n for cal_rsq > 0.33")
mtext(mtext(expression(paste("Observed Validation R"^"2")) , side = 2, line = 2.5))
mtext(mtext(expression(paste("Predicted R"^"2")), side = 1, line = 3))
abline(lm(y~x), lwd = 2)
mtext(paste0('R2 = ',round(summary(mod2)$r.squared,2)), side = 3, line = -2, adj = 0.05, cex = 1.2)

par(mfrow=c(1,2))

# cal.rsq and phylo level, both subsetted and not
mod3 <- lm(val_rsq ~ cal_rsq + phylo_level, data = d)
y <- d$val_rsq
x <- fitted(mod3)
plot(y ~ x, ylab = NA, xlab = NA, bty =  'l', pch = 16, cex = 1.5, main="val_rsq ~ cal_rsq + phylo_level")
mtext(mtext(expression(paste("Observed Validation R"^"2")) , side = 2, line = 2.5))
mtext(mtext(expression(paste("Predicted R"^"2")), side = 1, line = 3))
abline(lm(y~x), lwd = 2)
mtext(paste0('R2 = ',round(summary(mod3)$r.squared,2)), side = 3, line = -2, adj = 0.05, cex = 1.2)
mtext('b.', side = 1, line = -1.5, adj = 0.95, cex = 1.2)


mod4 <- lm(val_rsq ~ cal_rsq + phylo_level, data = d[d$cal_rsq > 0.33,])
y <- d[d$cal_rsq > 0.33,]$val_rsq
x <- fitted(mod4)
plot(y ~ x, ylab = NA, xlab = NA, bty =  'l', pch = 16, cex = 1.5, main = "val_rsq ~ cal_rsq + phylo_level,\n for cal_rsq > 0.33")
mtext(mtext(expression(paste("Observed Validation R"^"2")) , side = 2, line = 2.5))
mtext(mtext(expression(paste("Predicted R"^"2")), side = 1, line = 3))
abline(lm(y~x), lwd = 2)
mtext(paste0('R2 = ',round(summary(mod4)$r.squared,2)), side = 3, line = -2, adj = 0.05, cex = 1.2)

dev.off()
