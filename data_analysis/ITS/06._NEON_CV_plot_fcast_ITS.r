#Forecast using multinomial-dirichlet fit - NEON plot-level 50/50 cross validation.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/ddirch_forecast.r')

#set output path.----
output.path <- plot.CV_NEON_fcast.path

#load model and NEON site predictors..----
all.mod <- readRDS(plot.CV_NEON_dmulti.ddirch_all.path)
dat <- readRDS(plot.CV_NEON_cal.val_data.path) #calibration/validation dat core-level NEON.

#Define x_mu and x_sd values.----
plot.preds <- dat$val$x_mu.val
plot.sd    <- dat$val$x_sd.val

#run forecast over all phylo/functional levels.----
all.output <- list()
for(i in 1:length(all.mod)){
  mod <- all.mod[[i]]
  #drop covariate uncertainty.
 #plot.fit <- ddirch_forecast(mod, cov_mu = plot.preds, cov_sd = plot.sd, names = plot.preds$plotID)
  plot.fit <- ddirch_forecast(mod, cov_mu = plot.preds,                   names = plot.preds$plotID)
  #store output as a list and save.----
  output <- list(plot.fit,plot.preds,plot.sd)
  names(output) <- c('plot.fit','plot.preds','plot.sd')
  all.output[[i]] <- output
  cat(names(all.mod)[i],'forecast complete.',i,'of',length(all.mod),'forecasts completed.\n')
}
names(all.output) <- names(all.mod)

#Save output.----
saveRDS(all.output, output.path)
