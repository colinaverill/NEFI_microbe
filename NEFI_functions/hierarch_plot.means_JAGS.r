#' hierarch_plot.means_JAGS
#' Get site level means and sd from NEON data observed at the plot scale.
#' Also return global mean and sd for a given predictor.
#' depends on the runjags package.
#'
#' @param x_mu         #vector of observations to modeled.
#' @param x_sd         #vector of observation uncertainty. Default NA.
#' @param plot_site    #NEON plot IDs associated with x_mu observations.
#' @param n.adapt      #number of adaptive iterations. default 200.
#' @param n.burnin     #number of burnin iterations. default 1000.
#' @param n.sample     #number of sample iterations. default 1000.
#' @param n.chains     #number of MCMC chains. default 3.
#' @param parallel     #run chains in parallel or not? default F, notsetup to work yet.
#'
#' @return
#' @export
#'
#' @examples
#' rm(list=ls())
#' source('paths.r'); source('NEFI_functions/crib_fun.r')
#' dp.10098.plot <- readRDS(dp1.10098.00_plot.level.path)
#' to.ag <- dp.10098.plot[,c('plotID','siteID','relEM','basal_live')]
#' x_mu = to.ag$basal_live; plot_site=to.ag$siteID;n.adapt=200;n.burnin=1000;n.sample=1000;n.chains=3; parallel = F
hierarch_plot.means_JAGS <- function(x_mu, x_sd = NA, plot_site,
                                     n.adapt = 200, n.burnin = 1000, n.sample = 1000, n.chains = 3, 
                                     parallel = F, gamma_variance = F){
  #Other functions.
  source('NEFI_functions/precision_matrix_match.r')
  source('NEFI_functions/sd_to_precision.r')
  
  #Generate precision matrix.----
  #make sd object if not supplied.
  suppressWarnings(
    if(is.na(x_sd)){x_sd = rep(0.0001,length(x_mu))}
  )
  
  #covert sd to precision. output is matrix.
  x_precision <- sd_to_precision(as.matrix(x_sd))

  #### specify JAGS model ####
  jags.model = "
  model{
    #get site level means.
    for(i in 1:N.plot) {
      plot_la[i] ~ dnorm(site_mu[plot_site[i]], tau.site[plot_site[i]])
      plot_mu[i] ~ dnorm(plot_la[i], plot_precision[i]) #plot precision is observation uncertainty.
    }
    
    #get global mean and uncertainty.
    for(i in 1:N.site){
    site_mu[i] ~ dnorm(global_mu, tau.glob)
    }
    
    #Site level priors
    for(i in 1:N.site) {
      sigma.site[i] ~ dunif(0, 100)
        tau.site[i] <- pow(sigma.site[i], -2)
    }
  
    #Global level priors
    global_mu ~ dnorm(0,1E-3)  #I(0, ) #interval censored to be greater than zero.
    sigma.glob ~ dunif(0,100)
    tau.glob <- pow(sigma.glob, -2)
    
  }"
  
  #### JAGS model with gamma on variance prior. ####
  jags.gamma_variance.model = "
  model {
    #get site level means.
    for(i in 1:N.plot) {
      plot_mu[i] ~ dnorm(site_mu[plot_site[i]], tau.site[plot_site[i]])
    }
    
    #get global mean and uncertainty.
    for(i in 1:N.site){
      site_mu[i] ~ dnorm(global_mu, tau.glob)
    }
    
    #Site level priors
    for(i in 1:N.site) {
       tau.site[i]  ~ dgamma(1,1)
    }
  
    #Global level priors
    global_mu ~ dnorm(0,1E-3)  #I(0, ) #interval censored to be greater than zero.
    sigma.glob ~ dunif(0,100)
      tau.glob <- pow(sigma.glob, -2)
    
  }"
  
  #### setup JAGS data object. ####
  dat <- data.frame(x_mu, x_precision,plot_site)
  dat <- dat[complete.cases(dat),]
  site.names <- unique((plot_site))
  #reorder site.names so it matches what they are going to be used as.
  site.order <- order(site.names)
  site.names <- as.character(site.names[site.order])
  
  jags.data <- list(N.plot = length(dat$plot_site), N.site = length(unique(dat$plot_site)),
                    plot_mu = dat$x_mu,
                    plot_precision = dat$x_precision,
                    plot_site = droplevels(as.factor(dat$plot_site))
                    )
  
  #### fit JAGS model ####
  #runmode <- ifelse(parallel == T,T,F)
  j.model <- jags.model
  if(gamma_variance == T){j.model = jags.gamma_variance.model}
  mod  <- runjags::run.jags(model = j.model,
                            data = jags.data,
                            monitor = c('site_mu','global_mu'),
                            adapt = n.adapt,
                            burnin = n.burnin,
                            sample = n.sample,
                            n.chains = n.chains)
  out <- summary(mod)
  site.table <- data.frame(out[grep('site_mu', rownames(out)),])
  site.table$siteID <- site.names
  glob.table <- data.frame(t(out[grep('global_mu', rownames(out)),]))
  
  #### return output ####
  output.list <- list(out,site.table,glob.table)
  names(output.list) <- c('jags.summary','site.table','glob.table')
  return(output.list)
  
}
