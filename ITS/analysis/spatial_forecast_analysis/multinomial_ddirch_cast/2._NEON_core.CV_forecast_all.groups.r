#Forecast using multinomial-dirichlet fit - NEON core-level 50/50 cross validation.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/dmulti_ddirch_forecast.r')

#set output path.----
output.path <- core.CV_NEON_fcast.path

#load model and NEON site predictors..----
all.mod <- readRDS(core.CV_NEON_dmulti.ddirch_all.path)
dat <- readRDS(core.CV_NEON_cal.val_data.path) #calibration/validation dat core-level NEON.

#Define x_mu and x_sd values.----
core.preds <- dat$val$x_mu.val
core.sd    <- dat$val$x_sd.val

#run forecast over all phylo/functional levels.----
all.output <- list()
for(i in 1:length(all.mod)){
  mod <- all.mod[[i]]
  core.fit <- dmulti_ddirch_forecast(mod, cov_mu = core.preds, cov_sd = core.sd, names = core.preds$sampleID, make_it_work = T)
  #store output as a list and save.----
  output <- list(core.fit,core.preds,core.sd)
  names(output) <- c('core.fit','core.preds','core.sd')
  all.output[[i]] <- output
  cat(names(all.mod)[i],'forecast complete.',i,'of',length(all.mod),'forecasts completed.\n')
}
names(all.output) <- names(all.mod)

#Save output.----
saveRDS(all.output, output.path)
