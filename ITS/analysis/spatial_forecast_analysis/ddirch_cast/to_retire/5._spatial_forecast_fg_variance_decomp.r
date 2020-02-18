#spatial forecast variance decomposition: fungal functional groups.
#clear envrionment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/ddirch_forecast.r')


#set output path.----
output.path <- NEON_ddirch_var.decomp_fg.path

#load model results.----
#mod 1 is data from maps.
#mod 2 is site-specific data, no maps.
#mod 3 is all covariates.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)
mod <- mod[[3]] #just the all predictor case.

#get core-level covariate means and sd.----
dat <- readRDS(hierarch_filled.path)
core_mu <- dat$core.core.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu
#merge together.
plot_mu$siteID <- NULL
core.preds <- merge(core_mu   , plot_mu)
core.preds <- merge(core.preds, site_mu)
core.preds$relEM <- NULL
names(core.preds)[names(core.preds)=="b.relEM"] <- "relEM"

#get core-level SD.
core_sd <- dat$core.core.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together.
plot_sd$siteID <- NULL
core.sd <- merge(core_sd   , plot_sd)
core.sd <- merge(core.sd, site_sd)
core.sd$relEM <- NULL
names(core.sd)[names(core.sd)=="b.relEM"] <- "relEM"

#take median value core means and sd.
core_median_mu <- list()
core_median_sd <- list()
for(i in 1:ncol(core.preds)){
  if(is.numeric(core.preds[,i])){
    core_median_mu[[i]] <- median(core.preds[,i])
    core_median_sd[[i]] <- median(core.sd   [,i])
    names(core_median_mu)[[i]] <- colnames(core.preds)[i]
    names(core_median_sd)[[i]] <- colnames(core.sd   )[i]
  }
}
core_median_mu <- data.frame(t(unlist(core_median_mu)))
core_median_sd <- data.frame(t(unlist(core_median_sd)))

#get plot-level covariate means and sd.----
core_mu <- dat$core.plot.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu
#merge together, .
plot_mu$siteID <- NULL
plot.preds <- merge(core_mu,plot_mu)
plot.preds <- merge(plot.preds,site_mu)
plot.preds$relEM <- NULL
names(plot.preds)[names(plot.preds)=="b.relEM"] <- "relEM"

#get plot-level SD.
core_sd <- dat$core.plot.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together
plot_sd$siteID <- NULL
plot.sd <- merge(core_sd,plot_sd)
plot.sd <- merge(plot.sd,site_sd)
plot.sd$relEM <- NULL
names(plot.sd)[names(plot.sd)=='b.relEM'] <- "relEM"

#take median value core means and sd.
plot_median_mu <- list()
plot_median_sd <- list()
for(i in 1:ncol(plot.preds)){
  if(is.numeric(plot.preds[,i])){
    plot_median_mu[[i]] <- median(plot.preds[,i])
    plot_median_sd[[i]] <- median(plot.sd   [,i])
    names(plot_median_mu)[[i]] <- colnames(plot.preds)[i]
    names(plot_median_sd)[[i]] <- colnames(plot.sd   )[i]
  }
}
plot_median_mu <- data.frame(t(unlist(plot_median_mu)))
plot_median_sd <- data.frame(t(unlist(plot_median_sd)))


#get site-level covariate means and sd.----
core_mu <- dat$core.site.mu
plot_mu <- dat$plot.site.mu
site_mu <- dat$site.site.mu
#merge together
site.preds <- merge(core_mu, plot_mu)
site.preds <- merge(site.preds,site_mu)
names(site.preds)[names(site.preds)=='b.relEM'] <- "relEM"

#get site-level SD.
core_sd <- dat$core.site.sd
plot_sd <- dat$plot.site.sd
site_sd <- dat$site.site.sd
#merge together.
site.sd <- merge(core_sd,plot_sd)
site.sd <- merge(site.sd,site_sd)
names(site.sd)[names(site.sd)=='b.relEM'] <- "relEM"

#take median value core means and sd.
site_median_mu <- list()
site_median_sd <- list()
for(i in 1:ncol(site.preds)){
  if(is.numeric(site.preds[,i])){
    site_median_mu[[i]] <- median(site.preds[,i])
    site_median_sd[[i]] <- median(site.sd   [,i])
    names(site_median_mu)[[i]] <- colnames(site.preds)[i]
    names(site_median_sd)[[i]] <- colnames(site.sd   )[i]
  }
}
site_median_mu <- data.frame(t(unlist(site_median_mu)))
site_median_sd <- data.frame(t(unlist(site_median_sd)))


#core-scale variance decomposition.----
    core_ref <- ddirch_forecast(mod = mod, cov_mu = core_median_mu, cov_sd = core_median_sd, names = 'core_scale')
    core_cov <- ddirch_forecast(mod = mod, cov_mu = core_median_mu, cov_sd = core_median_sd, names = 'core_scale', zero_process_uncertainty = T, zero_parameter_uncertainty = T)
    core_par <- ddirch_forecast(mod = mod, cov_mu = core_median_mu, cov_sd = core_median_sd, names = 'core_scale', zero_process_uncertainty = T, zero_covariate_uncertainty = T)
    core_pro <- ddirch_forecast(mod = mod, cov_mu = core_median_mu, cov_sd = core_median_sd, names = 'core_scale', zero_parameter_uncertainty = T, zero_covariate_uncertainty = T)
core_cov.par <- ddirch_forecast(mod = mod, cov_mu = core_median_mu, cov_sd = core_median_sd, names = 'core_scale', zero_process_uncertainty = T)
core_cov.pro <- ddirch_forecast(mod = mod, cov_mu = core_median_mu, cov_sd = core_median_sd, names = 'core_scale', zero_parameter_uncertainty = T)
core_par.pro <- ddirch_forecast(mod = mod, cov_mu = core_median_mu, cov_sd = core_median_sd, names = 'core_scale', zero_covariate_uncertainty = T)
core_list <- list(core_ref, core_cov, core_par, core_pro, core_cov.par, core_cov.pro, core_par.pro)
decomp_names <- c('All','covariate','parameter','process','covariate + parameter','covariate + process','parameter + process')
core.decomp.list <- list()
for(i in 1:length(core_list)){
  core.decomp.list[[i]] <- core_list[[i]]$pi_0.975 - core_list[[i]]$pi_0.025
  core.decomp.list[[i]] <- core_list[[i]]$variance
}
core.decomp.list <- do.call(rbind, core.decomp.list)
rownames(core.decomp.list) <- decomp_names
core_variance_decomp <- sweep(core.decomp.list,2, core.decomp.list[1,],'/')

#plot-scale variance decomposition.----
    plot_ref <- ddirch_forecast(mod = mod, cov_mu = plot_median_mu, cov_sd = plot_median_sd, names = 'plot_scale')
    plot_cov <- ddirch_forecast(mod = mod, cov_mu = plot_median_mu, cov_sd = plot_median_sd, names = 'plot_scale', zero_process_uncertainty = T, zero_parameter_uncertainty = T)
    plot_par <- ddirch_forecast(mod = mod, cov_mu = plot_median_mu, cov_sd = plot_median_sd, names = 'plot_scale', zero_process_uncertainty = T, zero_covariate_uncertainty = T)
    plot_pro <- ddirch_forecast(mod = mod, cov_mu = plot_median_mu, cov_sd = plot_median_sd, names = 'plot_scale', zero_parameter_uncertainty = T, zero_covariate_uncertainty = T)
plot_cov.par <- ddirch_forecast(mod = mod, cov_mu = plot_median_mu, cov_sd = plot_median_sd, names = 'plot_scale', zero_process_uncertainty = T)
plot_cov.pro <- ddirch_forecast(mod = mod, cov_mu = plot_median_mu, cov_sd = plot_median_sd, names = 'plot_scale', zero_parameter_uncertainty = T)
plot_par.pro <- ddirch_forecast(mod = mod, cov_mu = plot_median_mu, cov_sd = plot_median_sd, names = 'plot_scale', zero_covariate_uncertainty = T)
plot_list <- list(plot_ref, plot_cov, plot_par, plot_pro, plot_cov.par, plot_cov.pro, plot_par.pro)
decomp_names <- c('All','covariate','parameter','process','covariate + parameter','covariate + process','parameter + process')
plot.decomp.list <- list()
for(i in 1:length(plot_list)){
  plot.decomp.list[[i]] <- plot_list[[i]]$pi_0.975 - plot_list[[i]]$pi_0.025
  plot.decomp.list[[i]] <- plot_list[[i]]$variance
}
plot.decomp.list <- do.call(rbind, plot.decomp.list)
rownames(plot.decomp.list) <- decomp_names
plot_variance_decomp <- sweep(plot.decomp.list,2, plot.decomp.list[1,],'/')

#site-scale variance decomposition.----
    site_ref <- ddirch_forecast(mod = mod, cov_mu = site_median_mu, cov_sd = site_median_sd, names = 'site_scale')
    site_cov <- ddirch_forecast(mod = mod, cov_mu = site_median_mu, cov_sd = site_median_sd, names = 'site_scale', zero_process_uncertainty = T, zero_parameter_uncertainty = T)
    site_par <- ddirch_forecast(mod = mod, cov_mu = site_median_mu, cov_sd = site_median_sd, names = 'site_scale', zero_process_uncertainty = T, zero_covariate_uncertainty = T)
    site_pro <- ddirch_forecast(mod = mod, cov_mu = site_median_mu, cov_sd = site_median_sd, names = 'site_scale', zero_parameter_uncertainty = T, zero_covariate_uncertainty = T)
site_cov.par <- ddirch_forecast(mod = mod, cov_mu = site_median_mu, cov_sd = site_median_sd, names = 'site_scale', zero_process_uncertainty = T)
site_cov.pro <- ddirch_forecast(mod = mod, cov_mu = site_median_mu, cov_sd = site_median_sd, names = 'site_scale', zero_parameter_uncertainty = T)
site_par.pro <- ddirch_forecast(mod = mod, cov_mu = site_median_mu, cov_sd = site_median_sd, names = 'site_scale', zero_covariate_uncertainty = T)
site_list <- list(site_ref, site_cov, site_par, site_pro, site_cov.par, site_cov.pro, site_par.pro)
decomp_names <- c('All','covariate','parameter','process','covariate + parameter','covariate + process','parameter + process')
site.decomp.list <- list()
for(i in 1:length(site_list)){
  site.decomp.list[[i]] <- site_list[[i]]$pi_0.975 - site_list[[i]]$pi_0.025
  site.decomp.list[[i]] <- site_list[[i]]$variance
}
site.decomp.list <- do.call(rbind, site.decomp.list)
rownames(site.decomp.list) <- decomp_names
site_variance_decomp <- sweep(site.decomp.list,2, site.decomp.list[1,],'/')

#wrap output and save.----
output <- list(core_variance_decomp, plot_variance_decomp, site_variance_decomp)
names(output) <- c('core_decomp','plot_decomp','site_decomp')
saveRDS(output,output.path)

#end script.
