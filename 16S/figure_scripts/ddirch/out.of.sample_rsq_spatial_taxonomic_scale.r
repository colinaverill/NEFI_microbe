# View NEON out-of-sample validation density plots by spatial scale 

rm(list=ls())
source('paths.r')
#source('NEFI_functions/zero_truncated_density.r')
# source dmulti_ddirch_forecast
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/zero_truncated_density.r", ssl.verifypeer = FALSE)
eval(parse(text = script))


#set output path.----
output.path <- paste0(scc_gen_16S_dir, "figures/out.of.sample_R2_spatial_tax_scale_ddirch.png")

#Load calibration data.----
#pl <- readRDS(bahram_16S_prior_ddirch_all.group_JAGSfits) #all phylo and functional groups.
pl <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/prior_phylo_fg_JAGSfit_16S.rds"))
fg <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/bahram_16S_prior_ddirch_fg_JAGSfits.rds"))
pl <- c(pl[1:5], fg)

#re-order list, make function group first.
#pl <- pl[c('fg','phylum','class','order','family','genus')]

#Get calibration r2 values.----
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
  #names(pl.r2)[[i]] <- names(pl)[i]
}

fg <- unname(pl.r2[6:17])
fg.r2.all <- list(unlist(fg))
names(fg.r2.all)<- "fg"
pl <- pl.r2[1:5]
names(pl) <- c('phylum','class','order','family','genus')
all.r2 <- c(fg.r2.all, pl)
cal.rsq <- unlist(all.r2)


cal.lev.mu <- unlist(lapply(all.r2, mean))
cal.lev.sd <- unlist(lapply(all.r2, sd))
cal.lev.N  <- unlist(lapply(all.r2, length))
cal.lev.se <- cal.lev.sd / sqrt(cal.lev.N)

#load forecasts predicted and observed.----
#pl.cast <- readRDS(NEON_cps_fcast_dmulti.ddirch_16S.path)
pl.cast <- readRDS(NEON_cps_fcast_ddirch_16S.path)
#pl.cast <- readRDS(paste0(pecan_gen_16S_dir, "/NEON_forecast_data/NEON_cps_fcast_ddirch_old.hier_16S.rds"))
pl.truth <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)
pl.core <- readRDS(NEON_16S_phylo_fg_abundances.path)
# 
# #Re-order.
# pl.cast <- pl.cast[c('fg','phylum','class','order','family','genus')]
# pl.truth <- pl.truth[c('fg','phylum','class','order','family','genus')]
# pl.core <- pl.core[c('fg','phylum','class','order','family','genus')]

#names.
#names(pl.cast)[1] <- 'functional'

#get core, plot site R2 values out of sample.----
all.core.rsq <- list()
all.plot.rsq <- list()
all.site.rsq <- list()
for(i in 1:length(pl.cast)){
  fcast <- pl.cast[[i]]
  core.rsq <- list()
  plot.rsq <- list()
  site.rsq <- list()
  #core.level----
  y <- pl.truth[[i]]$core.fit
  #y <- (pl.core[[i]]$abundances + 1)/pl.core[[i]]$seq_total
  x <- fcast$core.fit$mean
  #make sure row and column orders match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  truth <- y
  map <- readRDS(core_obs_16S.path)
  truth <- as.data.frame(truth)
  truth$deprecatedVialID <- rownames(truth)
  truth1 <- merge(truth, map[,c("deprecatedVialID", "geneticSampleID")], by = "deprecatedVialID")
  truth1 <- truth1[!duplicated(truth1$geneticSampleID),]
  rownames(truth1) <- gsub('-GEN','',truth1$geneticSampleID)
  truth <- truth1
  y <- truth
  
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(x)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    core.rsq[[k]] <- rsq
  }
  #plot.level----
  x <- fcast$plot.fit$mean
  y <- pl.truth[[i]]$plot.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('\\.','_',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(y)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    plot.rsq[[k]] <- rsq
  }
  #site.level----
  x <- fcast$site.fit$mean
  y <- pl.truth[[i]]$site.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(y)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    names(rsq) <- fungi_name
    site.rsq[[k]] <- rsq
  }
  #wrap up for return.----
  all.core.rsq[[i]] <- unlist(core.rsq)
  all.plot.rsq[[i]] <- unlist(plot.rsq)
  all.site.rsq[[i]] <- unlist(site.rsq)
  
  all.core.rsq[[i]] <-  all.core.rsq[[i]][-grep('other',names(all.core.rsq[[i]]))]
  all.plot.rsq[[i]] <-  all.plot.rsq[[i]][-grep('other',names(all.plot.rsq[[i]]))]
  all.site.rsq[[i]] <-  all.site.rsq[[i]][-grep('other',names(all.site.rsq[[i]]))]
  
}
core.rsq <- unlist(all.core.rsq)
plot.rsq <- unlist(all.plot.rsq)
site.rsq <- unlist(all.site.rsq)

# fix the functional-group grouping for each level
core.rsq.fg <- list(unlist(all.core.rsq[6:17]))
names(core.rsq.fg)<- "functional"
core.rsq.pl <- all.core.rsq[1:5]
names(core.rsq.pl) <- c('phylum','class','order','family','genus')
core.rsq <- c(core.rsq.fg, core.rsq.pl)
core.rsq <- unlist(core.rsq)

plot.rsq.fg <- list(unlist(all.plot.rsq[6:17]))
names(plot.rsq.fg)<- "functional"
plot.rsq.pl <- all.plot.rsq[1:5]
names(plot.rsq.pl) <- c('phylum','class','order','family','genus')
plot.rsq <- c(plot.rsq.fg, plot.rsq.pl)
plot.rsq <- unlist(plot.rsq)

site.rsq.fg <- list(unlist(all.site.rsq[6:17]))
names(site.rsq.fg)<- "functional"
site.rsq.pl <- all.site.rsq[1:5]
names(site.rsq.pl) <- c('phylum','class','order','family','genus')
site.rsq.list <- c(site.rsq.fg, site.rsq.pl)
site.rsq <- unlist(site.rsq.list)

# shouldn't we be getting these values *after* subsetting by calibration R2?
lev.mu <- unlist(lapply(site.rsq.list, mean  ))
lev.sd <- unlist(lapply(site.rsq.list, sd    ))
lev.N  <- unlist(lapply(site.rsq.list, length))
lev.se <- lev.sd / sqrt(lev.N)
names(lev.mu) <- names(site.rsq.list)

#Subset to observations that have a minimum calibration R2 value.----
pass <- cal.rsq[cal.rsq > .1]
core.rsq <- core.rsq[names(core.rsq) %in% names(pass)]
plot.rsq <- plot.rsq[names(plot.rsq) %in% names(pass)]
site.rsq <- site.rsq[names(site.rsq) %in% names(pass)]


# get densities
core.d <- zero_truncated_density(core.rsq)
plot.d <- zero_truncated_density(plot.rsq)
site.d <- zero_truncated_density(site.rsq)





#png save line.----
png(filename=output.path,width=8,height=5,units='in',res=300)

#global plot settings.----
par(mfrow = c(1,2))
limx <- c(0,1)
limy <- c(0, 8)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
cols <- c('purple','cyan','yellow')
par(mfrow = c(1,2), mar = c(5,4.2,1.5,1.5))

#Density plot.----
plot(site.d,xlim = c(0, 0.89), ylim = limy, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(site.d, col = adjustcolor(cols[3],trans))
polygon(plot.d, col = adjustcolor(cols[2],trans))
polygon(core.d, col = adjustcolor(cols[1],trans))
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext(expression(paste("Validation R"^"2")), side = 1, line = 2.5, cex = o.cex)
legend(x = 0.6, y = 6, legend = c('core','plot','site'), col ='black', pt.bg=adjustcolor(cols,trans), bty = 'n', pch = 22, pt.cex = 1.5)


png(filename=paste0(pecan_gen_16S_dir, "figures/cal.val_by.tax.scale_ddirch_16S.png"),width=4,height=5,units='in',res=300)

#Calibration/Validation rsq ~ function/phylo scale.----
x <- 1:length(lev.mu)
limy <- c(0,max(cal.lev.mu + cal.lev.se))
plot(lev.mu ~ x, cex = 2.5, ylim = limy, pch = 16, ylab = NA, xlab = NA, bty='n', xaxt = 'n', col="grey")
points(cal.lev.mu ~ x, cex = 2.5, pch = 16)
#arrows(x, lev.mu - lev.se, x1 = x, y1 = lev.mu + lev.se, length=0.00, angle=90, code=3, col = 'black')
lines(x, lev.mu, lty = 2, col="grey")
lines(x, cal.lev.mu, lty = 2)
mtext(expression(paste("Site-Level R"^"2")), side = 2, line = 2.2, cex = o.cex)
axis(1, labels = F)
text(x=x, y = -.04, labels= names(lev.mu), srt=45, adj=1, xpd=TRUE)
legend("bottomleft", 
       legend = c("calibration", "validation"), 
       col = c("black","grey"), 
       pch = 19, 
       bty = "n",
       pt.cex = 2)

#end plot.----
dev.off()
