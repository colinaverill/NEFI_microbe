#Comparing tedersoo vs. neonCV fits to core and plot scale data.
rm(list=ls())
source('paths.r')

#load forecasts, reformat for analyses.----
neo.core.cast <- readRDS(core.CV_NEON_fcast.path)
neo.plot.cast <- readRDS(plot.CV_NEON_fcast.path)
ted.cast <- readRDS(NEON_dmulti.ddirch_fcast_fg.path)

#break apart ted.cast to match NEON organization.
ted.core.cast <- list()
ted.plot.cast <- list()
for(i in 1:length(ted.cast)){
 ted.core.cast[[i]] <- list(ted.cast[[i]]$core.fit, ted.cast[[i]]$core.preds, ted.cast[[i]]$core.sd)
 ted.plot.cast[[i]] <- list(ted.cast[[i]]$plot.fit, ted.cast[[i]]$plot.preds, ted.cast[[i]]$plot.sd)
 names(ted.core.cast[[i]]) <- c('core.fit','core.preds','core.sd')
 names(ted.plot.cast[[i]]) <- c('plot.fit','plot.preds','plot.sd')
}
names(ted.core.cast) <- names(ted.cast)
names(ted.plot.cast) <- names(ted.cast)

#Load observed core and plot scale data, reformat for analysis.----
core.dat <- readRDS(core.CV_NEON_cal.val_data.path)
plot.dat <- readRDS(plot.CV_NEON_cal.val_data.path)

#Massage into nested by phylo/function level format.
neo.core.dat <- list()
neo.plot.dat <- list()
ted.core.dat <- list()
ted.plot.dat <- list()
for(i in 1:length(core.dat$cal$y.cal)){
  cal.core <- core.dat$cal$y.cal[[i]]$rel.abundances
  val.core <- core.dat$val$y.val[[i]]$rel.abundances
  cal.plot <- plot.dat$cal$y.cal[[i]]$mean
  val.plot <- plot.dat$val$y.val[[i]]$mean
  rownames(cal.core) <- gsub('-GEN','',rownames(cal.core))
  rownames(val.core) <- gsub('-GEN','',rownames(val.core))
  ted.core.dat[[i]] <- rbind(cal.core, val.core)
  ted.plot.dat[[i]] <- rbind(cal.plot, val.plot)
  neo.core.dat[[i]] <- val.core
  neo.plot.dat[[i]] <- val.plot
}
names(ted.core.dat) <- names(ted.cast)
names(ted.plot.dat) <- names(ted.cast)
names(neo.core.dat) <- names(ted.cast)
names(neo.plot.dat) <- names(ted.cast)

#Make model and data lists.
cast <- list(neo.core.cast, neo.plot.cast, ted.core.cast, ted.plot.cast)
data <- list(neo.core.dat , neo.plot.dat , ted.core.dat , ted.plot.dat )
names(cast) <- c('neo.core.cast','neo.plot.cast','ted.core.cast','ted.plot.cast')

#Get observed ~ predicted R2-best and R2-1:1 for each set of models.----
all.rsq <- list()
for(i in 1:length(cast)){
  this.cast <- cast[[i]]
  this.data <- data[[i]]
  this.rsq  <- list()
  for(j in 1:length(this.cast)){
    y <- this.data[[j]]
    x <- this.cast[[j]][[1]]$mean
    #make sure rownames are in same order
    x <- x[rownames(x) %in% rownames(y),]
    x <- x[order(match(rownames(x), rownames(y))),]
    #calculate R2 values
    rsq <- list()
    for(k in 1:ncol(x)){
      mod <- lm(y[,k] ~ x[,k])
      #rsq of best fit line.
      rsq.b <- summary(mod)$r.squared
      #resq of 1:1 line.
      rss <- sum((x[,k] -      y[,k])  ^ 2)  ## residual sum of squares
      tss <- sum((y[,k] - mean(y[,k])) ^ 2)  ## total sum of squares
      rsq.1 <- 1 - rss/tss
      rsq.1 <- ifelse(rsq.1 < 0, 0, rsq.1)
      return <- c(rsq.b, rsq.1)
      rsq[[k]] <- return
    }
    rsq <- do.call(rbind, rsq)
    rownames(rsq) <- colnames(x)
    colnames(rsq) <- c('rsq.b','rsq.1')
    this.rsq[[j]] <- rsq
  }
  names(this.rsq) <- names(this.cast)
  all.rsq[[i]] <- this.rsq
}
names(all.rsq) <- names(cast)

#Visualize differences in R2 by function/phylo group, scale (core vs. plot) and forecast (tedersoo vs. NEON.cv).-----
to_plot.rsq <- list()
for(i in 1:length(all.rsq)){
  this.rsq <- all.rsq[[i]]
  rsq.sum <- list()
  for(k in 1:length(this.rsq)){
   rsq.sum[[k]] <- mean(this.rsq[[k]][,1], na.rm = T)
  }
  rsq.sum <- unlist(rsq.sum)
  names(rsq.sum) <- names(this.rsq)
  to_plot.rsq[[i]] <- rsq.sum
}
names(to_plot.rsq) <- names(all.rsq)


