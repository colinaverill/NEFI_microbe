#Forecast using multinomial-dirichlet fit - NEON plot-level 50/50 cross validation.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/ddirch_forecast_noLogMap.r')

#set output path.----
output.path <- plot.CV_NEON_fcast_16S.path

#load model and NEON site predictors..----
all.mod <- readRDS(plot.CV_NEON_ddirch_16S_JAGSfit)
dat <- readRDS(plot.CV_NEON_cal.val_data_16S.path) #calibration/validation dat core-level NEON.

#Define x_mu and x_sd values.----
plot.preds <- dat$val$x_mu.val
plot.sd    <- dat$val$x_sd.val

#run forecast over all phylo/functional levels.----
all.output <- list()
for(i in 1:length(all.mod)){
  mod <- all.mod[[i]]
  plot.fit <- ddirch_forecast_noLogMap(mod, cov_mu = plot.preds, cov_sd = plot.sd, names = plot.preds$plotID)
  #store output as a list and save.----
  output <- list(plot.fit,plot.preds,plot.sd)
  names(output) <- c('plot.fit','plot.preds','plot.sd')
  all.output[[i]] <- output
  cat(names(all.mod)[i],'forecast complete.',i,'of',length(all.mod),'forecasts completed.\n')
}
names(all.output) <- names(all.mod)

#Save output.----
saveRDS(all.output, output.path)
