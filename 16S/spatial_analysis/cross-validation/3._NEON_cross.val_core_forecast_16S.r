#Forecast bacteria using multinomial-dirichlet fit - NEON core-level 50/50 cross validation.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/ddirch_forecast_noLogMap.r')

#set output path.----
output.path <- core.CV_NEON_fcast_16S.path

#load model and NEON site predictors..----
all.mod <- readRDS(core.CV_NEON_ddirch_16S_JAGSfit)
dat <- readRDS(core.CV_NEON_cal.val_data_16S.path) #calibration/validation dat core-level NEON.

#Define x_mu and x_sd values.----
core.preds <- dat$val$x_mu.val
core.sd    <- dat$val$x_sd.val

#run forecast over all phylo/functional levels.----
all.output <- list()
for(i in 1:length(all.mod)) {
  mod <- all.mod[[i]]
  core.fit <- ddirch_forecast_noLogMap(mod, cov_mu = core.preds, cov_sd = core.sd, names = core.preds$sampleID)
  #store output as a list and save.----
  output <- list(core.fit,core.preds,core.sd)
  names(output) <- c('core.fit','core.preds','core.sd')
  all.output[[i]] <- output
  cat(names(all.mod)[i],'forecast complete.',i,'of',length(all.mod),'forecasts completed.\n')
}
names(all.output) <- names(all.mod)

#Save output.----
saveRDS(all.output, output.path)
