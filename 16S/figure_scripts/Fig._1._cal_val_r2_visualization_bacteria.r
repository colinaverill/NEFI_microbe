#plott calibration and validation R2 distributions for Fungi and bacteria as 6-panel plot.
rm(list=ls())
source('paths.r')
source('paths_fall2019.r')
#source('NEFI_functions/zero_truncated_density.r')
source('NEFI_functions/rsq_1.1.r')
library(Metrics)
library(RCurl)

script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/zero_truncated_density.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/rsq_1.1.r", ssl.verifypeer = FALSE)
eval(parse(text = script))


#set output path.----
output.path <- paste0(scc_gen_16S_dir, "figures/Fig1_cal_val_scale.png")

#Load calibration data.----
cal <- readRDS(prior_delgado_ddirch_16S.path)  #dirichlet-only.

#Get calibration r2 values.----
cal.r2 <- list()
lev.cal.out <- list()
for(i in 1:length(cal)){
  lev <- cal[[i]]
  obs <- lev$observed
  pred <- lev$predicted
  lev.r2 <- list()
  for(j in 1:ncol(obs)){lev.r2[[j]] <- summary(lm(obs[,j] ~ pred[,j]))$r.squared}
  lev.r2 <- unlist(lev.r2)
  names(lev.r2) <- colnames(obs)
  lev.r2 <- lev.r2[names(lev.r2) != 'other']
  cal.r2[[i]] <- lev.r2
  #names(pl.r2)[[i]] <- names(pl)[i]
}
#names(cal.r2) <- names(cal)


cal.r2.fg <- list(unlist(cal.r2[6:18]))
cal.r2 <- c(cal.r2.fg, cal.r2[1:5])
names(cal.r2) <- c('functional','phylum','class','order','family','genus')

#cal.rsq <- unlist(all.r2)

#Subset to ones we can actually predict (r2 > 0.1).
lev.cal.r2 <- list()
for(i in 1:length(cal.r2)){
  lev.cal.r2[[i]] <- cal.r2[[i]][cal.r2[[i]] > 0.1]
}
names(lev.cal.r2) <- names(cal.r2)
cal.mu <- unlist(lapply(lev.cal.r2, mean  ))
cal.sd <- unlist(lapply(lev.cal.r2, sd    ))
cal.N  <- unlist(lapply(lev.cal.r2, length))
cal.se <- cal.sd / sqrt(cal.N)
names(cal.mu) <- names(cal.r2)
lev.cal.r2 <- unlist(lev.cal.r2)


#load forecasts predicted and observed.----
#val.cast <- readRDS(NEON_dmulti.ddirch_fcast_fg.path) 
val.truth <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)
val.cast <- readRDS(NEON_cps_fcast_ddirch_16S.path)
map <- readRDS(core_obs_16S.path)

#get core, plot site R2 and RMSE values out of sample.----
all.core.rsq  <- list()
all.plot.rsq  <- list()
all.site.rsq  <- list()
all.core.rmse <- list()
all.plot.rmse <- list()
all.site.rmse <- list()
all.core.stat <- list()
all.plot.stat <- list()
all.site.stat <- list()
for(i in 1:length(val.cast)){
  fcast <- val.cast[[i]]
  core.rsq <- list()
  plot.rsq <- list()
  site.rsq <- list()
  core.rmse <- list()
  plot.rmse <- list()
  site.rmse <- list()
  core.stat <- list()
  plot.stat <- list()
  site.stat <- list()
  #core.level----
  y <- val.truth[[i]]$core.fit
  x <- fcast$core.fit$mean
  #fix rownames
  y.fix <- as.data.frame(y)
  y.fix$deprecatedVialID <- rownames(y.fix)
  y.fix <- merge(y.fix, map[,c("deprecatedVialID", "geneticSampleID")], by = "deprecatedVialID")
  y.fix <- y.fix[!duplicated(y.fix$geneticSampleID),]
  rownames(y.fix) <- gsub('-GEN','',y.fix$geneticSampleID )
  y.fix$geneticSampleID <- NULL
  y <- y.fix
  #make sure rownames and col
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #fit model, grab r2.
  for(k in 1:ncol(fcast$core.fit$mean)){
    fungi_name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    obs.rmse <- rmse(y[,k], x[,k])
    rsq.1 <- rsq_1.1(y[,k], x[,k])
    names(rsq) <- fungi_name
    if(fungi_name == 'other'){next}
    core.rsq [[k]] <- rsq
    core.rmse[[k]] <- obs.rmse
    stat.return <- c(fungi_name, rsq, rsq.1, obs.rmse)
    names(stat.return) <- c('name','rsq','rsq.1','rmse')
    core.stat[[k]] <- stat.return
  }
  #plot.level----
  x <- fcast$plot.fit$mean
  y <- val.truth[[i]]$plot.fit$mean
  
  #fix rownames
  rownames(y) <- gsub('\\.','_',rownames(y))
  #make sure row and column order match.
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
    rsq.1 <- rsq_1.1(y[,k], x[,k])
    obs.rmse <- rmse(y[,k], x[,k])
    names(obs.rmse) <- names(rsq)
    names(rsq) <- fungi_name
    if(fungi_name == 'other'){next}
    plot.rsq [[k]] <- rsq
    plot.rmse[[k]] <- obs.rmse
    stat.return <- c(fungi_name, rsq, rsq.1, obs.rmse)
    names(stat.return) <- c('name','rsq','rsq.1','rmse')
    plot.stat[[k]] <- stat.return
  }
  #site.level----
  x <- fcast$site.fit$mean
  y <- val.truth[[i]]$site.fit$mean
  #make sure row and column order match.
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
    rsq.1 <- rsq_1.1(y[,k], x[,k])
    obs.rmse <- rmse(y[,k],x[,k])
    names(obs.rmse) <- names(rsq)
    names(rsq) <- fungi_name
    if(fungi_name == 'other'){next}
    site.rsq [[k]] <- rsq
    site.rmse[[k]] <- obs.rmse
    stat.return <- c(fungi_name, rsq, rsq.1, obs.rmse)
    names(stat.return) <- c('name','rsq','rsq.1','rmse')
    site.stat[[k]] <- stat.return
  }
  #wrap up for return.----
  all.core.rsq [[i]] <- unlist(core.rsq)
  all.plot.rsq [[i]] <- unlist(plot.rsq)
  all.site.rsq [[i]] <- unlist(site.rsq)
  all.core.rmse[[i]] <- unlist(core.rmse)
  all.plot.rmse[[i]] <- unlist(plot.rmse)
  all.site.rmse[[i]] <- unlist(site.rmse)
  all.core.stat[[i]] <- do.call(rbind, core.stat)
  all.plot.stat[[i]] <- do.call(rbind, plot.stat)
  all.site.stat[[i]] <- do.call(rbind, site.stat)
}

# fix functional groupings for output
#core-level
core.rsq.fg <- list(unlist(all.core.rsq[6:18]))
all.core.rsq <- c(core.rsq.fg, all.core.rsq[1:5])
core.rmse.fg <- list(unlist(all.core.rmse[6:18]))
all.core.rmse <- c(core.rmse.fg, all.core.rmse[1:5])
names(all.core.rmse) <- c('functional','phylum','class','order','family','genus')
names(all.core.rsq) <- c('functional','phylum','class','order','family','genus')
#plot-level
plot.rsq.fg <- list(unlist(all.plot.rsq[6:18]))
all.plot.rsq <- c(plot.rsq.fg, all.plot.rsq[1:5])
plot.rmse.fg <- list(unlist(all.plot.rmse[6:18]))
all.plot.rmse <- c(plot.rmse.fg, all.plot.rmse[1:5])
names(all.plot.rsq) <- c('functional','phylum','class','order','family','genus')
names(all.plot.rmse) <- c('functional','phylum','class','order','family','genus')
#site-level
site.rsq.fg <- list(unlist(all.site.rsq[6:18]))
all.site.rsq <- c(site.rsq.fg, all.site.rsq[1:5])
site.rmse.fg <- list(unlist(all.site.rmse[6:18]))
all.site.rmse <- c(site.rmse.fg, all.site.rmse[1:5])
names(all.site.rsq) <- c('functional','phylum','class','order','family','genus')
names(all.site.rmse) <- c('functional','phylum','class','order','family','genus')

# means by taxonomic rank
val.mu <- unlist(lapply(all.site.rsq, mean  ))
val.sd <- unlist(lapply(all.site.rsq, sd    ))
val.N  <- unlist(lapply(all.site.rsq, length))
val.se <- val.sd / sqrt(val.N)

names(val.mu) <- names(all.site.rsq)
core.rsq <- unlist(all.core.rsq)
plot.rsq <- unlist(all.plot.rsq)
site.rsq <- unlist(all.site.rsq)
core.rmse <- unlist(all.core.rmse)
plot.rmse <- unlist(all.plot.rmse)
site.rmse <- unlist(all.site.rmse)
names(core.rmse) <- names(core.rsq)
names(plot.rmse) <- names(plot.rsq)
names(site.rmse) <- names(site.rsq)

#Subset to observations that have a minimum calibration R2 value.----
pass <- lev.cal.r2
core.rsq <- core.rsq[names(core.rsq) %in% names(pass)]
plot.rsq <- plot.rsq[names(plot.rsq) %in% names(pass)]
site.rsq <- site.rsq[names(site.rsq) %in% names(pass)]
core.rmse <- core.rmse[names(core.rmse) %in% names(pass)]
plot.rmse <- plot.rmse[names(plot.rmse) %in% names(pass)]
site.rmse <- site.rmse[names(site.rmse) %in% names(pass)]
core.d <- zero_truncated_density(core.rsq)
plot.d <- zero_truncated_density(plot.rsq)
site.d <- zero_truncated_density(site.rsq)
core.rmse.d <- zero_truncated_density(core.rmse)
plot.rmse.d <- zero_truncated_density(plot.rmse)
site.rmse.d <- zero_truncated_density(site.rmse)

#png save line.----
png(filename=output.path,width=8,height=6,units='in',res=300)

#global plot settings.----
#par(mfrow = c(2,3), mar = c(5,4,2,1), oma = c(1,1,1,1))
par(mfrow = c(1,3), mar = c(5,4,2,1), oma = c(1,1,1,1))
limx <- c(0,1)
limy <- c(0, 5.1)
trans <- 0.2 #shading transparency.
o.cex <- 1.0 #outer label size.
cols <- c('purple','cyan','yellow')

#Calibration R2 denisty plot (retired).-----
#dat <- unlist(cal.r2)
#dat <- zero_truncated_density(dat)
#limy <- c(0, max(dat$y)*1.05)
#limx <- c(0, max(dat$x))
##Density plot.
#plot(dat,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
#polygon(dat, col = adjustcolor(cols[3],trans))
#mtext('Density', side = 2, line = 2.5, cex = o.cex)
#mtext(expression(paste("Calibration R"^"2")), side = 1, line = 3, cex = o.cex)
#mtext('(a)', side = 3, adj = 0.95, line = -2)



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
mtext('(d)', side = 3, adj = 0.95, line = -2)

#Validation RMSE denisty plot.-----
dat.a <- core.rmse.d
dat.b <- plot.rmse.d
dat.c <- site.rmse.d
limy <- c(0, max(c(dat.a$y, dat.b$y, dat.c$y))*1.05)
limx <- c(0, max(c(dat.a$x, dat.b$x, dat.c$x)))
#Density plot.
plot(dat.a,xlim = limx, ylim = limy, bty = 'l', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
polygon(dat.a, col = adjustcolor(cols[1],trans))
polygon(dat.b, col = adjustcolor(cols[2],trans))
polygon(dat.c, col = adjustcolor(cols[3],trans))
mtext('Density', side = 2, line = 2.5, cex = o.cex)
mtext(expression(paste("Validation RMSE")), side = 1, line = 3, cex = o.cex)
#legend
legend(x = 0.5, y = 12.5, legend = c('core','plot','site'), col ='black', pt.bg=adjustcolor(cols,trans), bty = 'n', pch = 22, pt.cex = 1.5, cex = 1.2)
mtext('(e)', side = 3, adj = 0.95, line = -2)


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

text(x=x+0.05, y = -.03, labels= names(cal.mu), srt=45, adj=1, xpd=TRUE, cex = 1.25)
#legend
legend(x = 1, y = 0.1, legend = c('calibration','validation'), col =c('black','gray'), bty = 'n', pch = 16, pt.cex = 2.5, cex = 1.2, y.intersp = 2.6)
mtext('(f)', side = 3, adj = 0.95, line = -2)

#Outer labels.----
mtext('Bacteria'   , cex = 1.5, side = 3, outer = T, adj = 0.01, line = -1)
#mtext('Bacteria', cex = 1.5, side = 3, outer = T, adj = 0.01, line = -23)

#end plot.----
dev.off()
