#Summarize model calibration and validation fits.
#To be run after spatial calibration and forecast scripts have been run.
#Report calibration and validations rsq, rsq.1, RMSE, N.
rm(list=ls())
source('paths.r')
source('NEFI_functions/rsq_1.1.r')
library(Metrics)

#set output path.----
output.path <- NEON_dmilti.ddirch_analysis_summary.path

#Load calibration data.----
cal <- readRDS(ted_ITS_prior_all.groups_JAGSfits.path)  #dirichlet only on transformed relative abundances.

#re-order list, make function group first.
cal <- cal[c('fg','phylum','class','order','family','genus')]
names(cal)[1] <- 'functional'



#Get calibration statistics.----
cal.stat <- list()
lev.cal.out <- list()
for(i in 1:length(cal)){
  lev <- cal[[i]]
  obs <- lev$observed
  pred <- lev$predicted
  #normalize counts
  obs <- obs/rowSums(obs)
  pred <- pred / rowSums(pred)
  lev.stat <- list()
  for(j in 1:ncol(obs)){
    name <- colnames(obs)[j]
    if(name == 'other'){next}
    rsq   <- summary(lm(obs[,j] ~ pred[,j]))$r.squared
    rsq.1 <- rsq_1.1(obs[,j], pred[,j])
    rmse  <- rmse(obs[,j], pred[,j])
    abundance <- mean(obs[,j])
    rmse.norm <- rmse / abundance
    return <- c(name, rsq, rsq.1, rmse,abundance, rmse.norm)
    names(return) <- c('name','rsq','rsq.1','rmse','abundance','rmse.norm')
    lev.stat[[j]] <- return
    }
  lev.stat <- data.frame(do.call(rbind, lev.stat))
  lev.stat[,2:ncol(lev.stat)] <- apply(lev.stat[,2:ncol(lev.stat)], 2, as.character)
  lev.stat[,2:ncol(lev.stat)] <- apply(lev.stat[,2:ncol(lev.stat)], 2, as.numeric  )
  cal.stat[[i]] <- lev.stat
}
names(cal.stat) <- names(cal)

#Get 'predictable' calibration subset (rsq.1 > 0.1).----
cal.stat.predictable <- list()
for(i in 1:length(cal.stat)){
  cal.stat.predictable[[i]] <- cal.stat[[i]][cal.stat[[i]]$rsq > 0.1,]
}
names(cal.stat.predictable) <- names(cal.stat)

#Get calibration summaries for all and predictable taxa.----
#all taxa.
cal.stat.sum <- list()
for(i in 1:length(cal.stat)){
  rsq   <-   mean(cal.stat[[i]]$rsq  )
  rsq.1 <-   mean(cal.stat[[i]]$rsq.1)
  rmse  <-   mean(cal.stat[[i]]$rmse )
  rmse.norm <- mean(cal.stat[[i]]$rmse.norm)
  N     <- length(cal.stat[[i]]$rsq  )
  stat <- c(rsq, rsq.1, rmse, N, rmse.norm)
  cal.stat.sum[[i]] <- stat
}
cal.stat.sum <- data.frame(do.call(rbind, cal.stat.sum))
rownames(cal.stat.sum) <- names(cal.stat)
colnames(cal.stat.sum) <- c('rsq','rsq.1','rmse','N','rmse.norm')

#predictable taxa.
cal.stat.predictable.sum <- list()
for(i in 1:length(cal.stat)){
  rsq   <-   mean(cal.stat.predictable[[i]]$rsq  )
  rsq.1 <-   mean(cal.stat.predictable[[i]]$rsq.1)
  rmse  <-   mean(cal.stat.predictable[[i]]$rmse )
  N     <- length(cal.stat.predictable[[i]]$rsq  )
  rmse.norm <- mean(cal.stat.predictable[[i]]$rmse.norm)
  stat <- c(rsq, rsq.1, rmse, N, rmse.norm)
  cal.stat.predictable.sum[[i]] <- stat
}
cal.stat.predictable.sum <- data.frame(do.call(rbind, cal.stat.predictable.sum))
rownames(cal.stat.predictable.sum) <- names(cal.stat.predictable)
colnames(cal.stat.predictable.sum) <- c('rsq','rsq.1','rmse','N','rmse.norm')

#wrap all calibration data.----
calibration <- list(cal.stat, cal.stat.predictable, cal.stat.sum, cal.stat.predictable.sum)
names(calibration) <- c('cal.stat','cal.stat.predictable','cal.stat.sum','cal.stat.predictable.sum')

#Load validation data.----
val.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)
val.cast  <- readRDS(NEON_site_fcast_all_groups_1k_rare.path)
#Re-order.
val.cast  <- val.cast [c('fg','phylum','class','order','family','genus')]
val.truth <- val.truth[c('fg','phylum','class','order','family','genus')]
names(val.cast)[1] <- 'functional'


#get validation stats @ core, plot site scale.----
core.stat <- list()
plot.stat <- list()
site.stat <- list()
for(i in 1:length(val.cast)){
  fcast <- val.cast[[i]]
  core.lev.stat <- list()
  plot.lev.stat <- list()
  site.lev.stat <- list()
  #core.level----
  y <- val.truth[[i]]$core.fit
  x <- fcast$core.fit$mean
  #make sure row and column orders match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #grab statistics for all taxa within level.
  #core level.----
  for(k in 1:ncol(x)){
    name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    rsq.1 <- rsq_1.1(y[,k], x[,k])
    obs.rmse <- rmse(y[,k], x[,k])
    abundance.y <- mean(y[,k])
    abundance.x <- mean(x[,x])
    rmse.norm <- obs.rmse / abundance.y
    if(name == 'other'){next}
    stat.return <- c(name, rsq, rsq.1, obs.rmse,rmse.norm,abundance.y,abundance.x)
    names(stat.return) <- c('name','rsq','rsq.1','rmse','rmse.norm','abundance.y','abundance.x')
    core.lev.stat[[k]] <- stat.return
  }
  #plot.level----
  x <- fcast$plot.fit$mean
  y <- val.truth[[i]]$plot.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #grab statistics for all taxa within level.
  for(k in 1:ncol(x)){
    name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    rsq.1 <- rsq_1.1(y[,k], x[,k])
    obs.rmse <- rmse(y[,k], x[,k])
    abundance.y <- mean(y[,k])
    abundance.x <- mean(x[,k])
    rmse.norm <- obs.rmse / abundance.y
    if(name == 'other'){next}
    stat.return <- c(name, rsq, rsq.1, obs.rmse,rmse.norm,abundance.y,abundance.x)
    names(stat.return) <- c('name','rsq','rsq.1','rmse','rmse.norm','abundance.y','abundance.x')
    plot.lev.stat[[k]] <- stat.return
  }
  #site.level----
  x <- fcast$site.fit$mean
  y <- val.truth[[i]]$site.fit$mean
  #make sure row and column order match.
  rownames(y) <- gsub('-GEN','',rownames(y))
  y <- y[rownames(y) %in% rownames(x),]
  x <- x[rownames(x) %in% rownames(y),]
  y <- y[,colnames(y) %in% colnames(x)]
  x <- x[,colnames(x) %in% colnames(y)]
  x <- x[order(match(rownames(x),rownames(y))),]
  x <- x[,order(match(colnames(x),colnames(y)))]
  #grab statistics for all taxa within level.
  for(k in 1:ncol(x)){
    name <- colnames(x)[k]
    rsq <- summary(lm(y[,k] ~ x[,k]))$r.squared
    rsq.1 <- rsq_1.1(y[,k], x[,k])
    obs.rmse <- rmse(y[,k], x[,k])
    abundance.y <- mean(y[,k])
    abundance.x <- mean(x[,k])
    rmse.norm <- obs.rmse / abundance.y
    if(name == 'other'){next}
    stat.return <- c(name, rsq, rsq.1, obs.rmse,rmse.norm,abundance.y,abundance.x)
    names(stat.return) <- c('name','rsq','rsq.1','rmse','rmse.norm','abundance.y','abundance.x')
    site.lev.stat[[k]] <- stat.return
  }
  #wrap up for return.----
  core.stat[[i]] <- data.frame(do.call(rbind, core.lev.stat))
  plot.stat[[i]] <- data.frame(do.call(rbind, plot.lev.stat))
  site.stat[[i]] <- data.frame(do.call(rbind, site.lev.stat))
  core.stat[[i]][,2:ncol(core.stat[[i]])] <- apply(core.stat[[i]][,2:ncol(core.stat[[i]])], 2, as.character)
  core.stat[[i]][,2:ncol(core.stat[[i]])] <- apply(core.stat[[i]][,2:ncol(core.stat[[i]])], 2, as.numeric  )
  plot.stat[[i]][,2:ncol(plot.stat[[i]])] <- apply(plot.stat[[i]][,2:ncol(plot.stat[[i]])], 2, as.character)
  plot.stat[[i]][,2:ncol(plot.stat[[i]])] <- apply(plot.stat[[i]][,2:ncol(plot.stat[[i]])], 2, as.numeric  )
  site.stat[[i]][,2:ncol(site.stat[[i]])] <- apply(site.stat[[i]][,2:ncol(site.stat[[i]])], 2, as.character)
  site.stat[[i]][,2:ncol(site.stat[[i]])] <- apply(site.stat[[i]][,2:ncol(site.stat[[i]])], 2, as.numeric  )
}
#name the lists w/ their function/taxonomic levels.
names(core.stat) <- names(val.cast)
names(plot.stat) <- names(val.cast)
names(site.stat) <- names(val.cast)

#Get 'predictable' validation subset.----
core.stat.predictable <- list()
plot.stat.predictable <- list()
site.stat.predictable <- list()
for(i in 1:length(core.stat)){
  core <- core.stat[[i]]
  plot <- plot.stat[[i]]
  site <- site.stat[[i]]
  ref  <- cal.stat.predictable[[i]]
  core.stat.predictable[[i]] <- core[core$name %in% ref$name,]
  plot.stat.predictable[[i]] <- plot[plot$name %in% ref$name,]
  site.stat.predictable[[i]] <- site[site$name %in% ref$name,]
}
names(core.stat.predictable) <- names(core.stat)
names(plot.stat.predictable) <- names(plot.stat)
names(site.stat.predictable) <- names(site.stat)

#Get validation summaries @# site level for all and predictable taxa.----
#all validation summary - site level.
site.stat.sum <- list()
for(i in 1:length(site.stat)){
  rsq   <- mean(site.stat[[i]]$rsq  )
  rsq.1 <- mean(site.stat[[i]]$rsq.1)
  rmse  <- mean(site.stat[[i]]$rmse )
  N     <- nrow(site.stat[[i]])
  rmse.norm <- mean(site.stat[[i]]$rmse.norm)
  stat  <- c(rsq, rsq.1, rmse, N, rmse.norm)
  names(stat) <- c('rsq','rsq.1','rmse','N','rmse.norm')
  site.stat.sum[[i]] <- stat
}
site.stat.sum <- data.frame(do.call(rbind, site.stat.sum))
site.stat.sum <- cbind(names(site.stat), site.stat.sum)
colnames(site.stat.sum)[1] <- 'level'

#predictable validation summary - site level.
site.stat.predictable.sum <- list()
for(i in 1:length(site.stat)){
  rsq   <- mean(site.stat.predictable[[i]]$rsq  )
  rsq.1 <- mean(site.stat.predictable[[i]]$rsq.1)
  rmse  <- mean(site.stat.predictable[[i]]$rmse )
  N     <- nrow(site.stat.predictable[[i]])
  rmse.norm <- mean(site.stat.predictable[[i]]$rmse.norm)
  stat  <- c(rsq, rsq.1, rmse, N, rmse.norm)
  names(stat) <- c('rsq','rsq.1','rmse','N','rmse.norm')
  site.stat.predictable.sum[[i]] <- stat
}
site.stat.predictable.sum <- data.frame(do.call(rbind, site.stat.predictable.sum))
site.stat.predictable.sum <- cbind(names(site.stat), site.stat.predictable.sum)
colnames(site.stat.predictable.sum)[1] <- 'level'

#wrap all validation data.----
val.stat <- list(core.stat, plot.stat, site.stat)
names(val.stat) <- c('core.stat','plot.stat','site.stat')
val.stat.predictable <- list(core.stat.predictable, plot.stat.predictable, site.stat.predictable)
names(val.stat.predictable) <- c('core.stat','plot.stat','site.stat')
validation <- list(val.stat,val.stat.predictable,site.stat.sum, site.stat.predictable.sum)
names(validation) <- c('val.stat','val.stat.predictable','site.stat.sum','site.stat.predictable.sim')

#save output.----
output <- list(calibration, validation)
names(output) <- c('calibration','validation')
saveRDS(output, output.path)

#end script.
