#' ddirch_obs_uncertainty
#' This function takes dirichlet variances at the site, plot and core level.
#' It also requires a matrix of species abundances to represent "true" values at the site level.
#' It starts at the site level, and draws plots within site from the site mean, and within site-scale variance.
#' It then draws cores within plots, given the "true" plot means, and within plot-scale variance.
#' It then draws estimates of observed core values, given the "true" core mean and within core variance.
#' Once all this data is drawn, the function hierarchically estimates core, plot and site means (y_obs)
#' It then regresses y_obs at each scale against the "true" values that were drawn, and estimates an R2 value.
#' This is repeated some number of times to generate a distribution of potential R2 values.
#' This reflects the maximum predictive accuracy once could expect, given observation uncertainty, in principle a perfect model would predict the "true" values generated here.
#' Predictaibility likely varies as a function of abundance, as well as variance, and their interaction.
#' Predictiability will also vary as a function scale.
#' NOTE: larger variance values lead to less variance in the dirichlet.
#'
#' @param site_mu     #matrix of site-scale relative abundances. Number of rows = number of sites.
#' @param site.var    #within site-scale variance.
#' @param plot.var    #within plot-scale variance.
#' @param core.var    #within core-scale variance.
#' @param n.sim       #numbre of simulations to run.
#' @param n.plot      #Number of plots within site. Default 10 for NEON.
#' @param n.core      #Number of cores within plot. Default  3 for NEON.
#' @param n.threads   #Number of cores to run in parallel. Default will auto-detect.
#'
#' @return
#' @export
#'
#' @examples
ddirch_obs_uncertainty <- function(site_mu,
                                   site.var, plot.var, core.var,
                                   n.sim = 100, n.plot = 10, n.core = 3, n.threads = NA){
  #function tests to throw errors.----
  check <- "DirichletReg" %in% rownames(installed.packages())
  if(check == F){stop('The package DirichletReg is not installed. We need it to be.')}
  check <- "doParallel" %in% rownames(installed.packages())
  if(check == F){stop('The package doParallel is not installed. We need it to be.')}
  library(doParallel)
  
  #setup parallel.----
  if( is.na(n.threads)){
    n.cores <- detectCores()
    registerDoParallel(n.cores)
  }
  if(!is.na(n.threads)){
    registerDoParallel(n.threads)
  }
  
  
  #execute loop in parallel.----
  all.output <- 
    foreach(j = 1:n.sim) %dopar% {
      #begin try-catch loop (sometimes we fail some type of hessian matrix calc...)----
      attempt = 0
      mod.plot <- NULL
      mod.site <- NULL
      while(is.null(mod.plot) && is.null(mod.site) && attempt <= 10){
        attempt = attempt + 1
        #1. Generate plot and core values.----
        #assign site labels.
        n.site <- nrow(site_mu)
        sites <- letters[1:n.site]
        
        #generate plot values.
        y.plot <- list()
        for(i in 1:nrow(site_mu)){
          plots <- DirichletReg::rdirichlet(n.plot, site_mu[i,]*site.var)
          #kill zeros.
          plots[plots == 0] <- min(plots[plots > 0])
          plot.lab <- paste0(sites[i],c(1:n.plot))
          plots <- data.frame(plot.lab, plots)
          colnames(plots) <- c('plotID',colnames(site_mu))
          y.plot[[i]] <- plots
        }
        y.plot <- do.call(rbind, y.plot)
        y.plot$siteID <- substr(y.plot$plotID, 1,1)
        
        #Generate core values.
        y.core <- list()
        for(i in 1:nrow(y.plot)){
          cores <- DirichletReg::rdirichlet(n.core, as.numeric(y.plot[i,grep('y',colnames(y.plot))]*plot.var))
          #kill zeros
          cores[cores == 0] <- min(cores[cores > 0])
          core.lab <- paste0(y.plot[i,'plotID'], '.', c(1:n.core))
          cores <- data.frame(core.lab, cores)
          colnames(cores) <- c('coreID',colnames(site_mu))
          y.core[[i]] <- cores
        }
        y.core <- do.call(rbind, y.core)
        y.core$coreID <- as.character(y.core$coreID)
        y.core$plotID <- substr(y.core$coreID, 1, nchar(y.core$coreID) - 2)
        y.core$siteID <- substr(y.core$coreID, 1,1)
        
        #Generate observed data by drawing from intra-core variance.
        y.obs <- DirichletReg::rdirichlet(nrow(y.core), as.matrix(y.core[,grep('y',colnames(y.core))]) * core.var)
        y.obs <- y.obs + min(y.obs[y.obs > 0])
        y.obs <- data.frame(y.obs)
        colnames(y.obs) <- colnames(y.core)[grep('y',colnames(y.core))]
        y.obs$coreID <- y.core$coreID
        y.obs$plotID <- substr(y.obs$coreID, 1, nchar(y.obs$coreID) - 2)
        
        #2. Aggregate to plot and site scale.----
        #I should maybe use JAGS here instead...
        #Plot scale aggregation.
        suppressWarnings(
          y.obs$Y <- DR_data(y.obs[,grep('y',colnames(y.obs))])
        )
        try(mod.plot <- DirichletReg::DirichReg(Y ~ plotID - 1, data = y.obs))
        plot.obs <- matrix(mod.plot$coefficients, ncol = ncol(site_mu))
        plot.obs <- exp(plot.obs)
        plot.obs <- data.frame(plot.obs / rowSums(plot.obs))
        colnames(plot.obs) <- colnames(site_mu)
        plot.obs$plotID <- unique(y.obs$plotID)
        plot.obs$siteID <- substr(plot.obs$plotID,1,1)
        
        #Site scale aggregation.
        suppressWarnings(
          plot.obs$Y <- DR_data(plot.obs[,grep('y',colnames(plot.obs))])
        )
        try(mod.site <- DirichletReg::DirichReg(Y ~ siteID - 1, data = plot.obs))
        site.obs <- matrix(mod.site$coefficients, ncol = ncol(site_mu))
        site.obs <- exp(site.obs)
        site.obs <- data.frame(site.obs / rowSums(site.obs))
        colnames(site.obs) <- colnames(site_mu)
        site.obs$siteID <- unique(plot.obs$siteID)
        #end try-catch loop.----
      }
      #3. Calculate r-squared values of 'true' vs. observed at core-plot-site scales.----
      #core-scale rsq values.
      core.rsq <- list()
      for(i in 1:ncol(site_mu)){
        pred <-  y.obs[,grep('y',colnames(y.obs))]
        obs <- y.core[,grep('y',colnames(y.core))]
        mod <- lm(obs[,i] ~ pred[,i])
        core.rsq[[i]] <- summary(mod)$r.squared
      }
      core.rsq <- unlist(core.rsq)
      names(core.rsq) <- colnames(site_mu)
      
      #plot-scale rsq values.
      plot.rsq <- list()
      for(i in 1:ncol(site_mu)){
        pred <-  plot.obs[,grep('y',colnames(plot.obs))]
        obs <- y.plot[,grep('y',colnames(y.plot))]
        mod <- lm(obs[,i] ~ pred[,i])
        plot.rsq[[i]] <- summary(mod)$r.squared
      }
      plot.rsq <- unlist(plot.rsq)
      names(plot.rsq) <- colnames(site_mu)
      
      #site-scale rsq values.
      site.rsq <- list()
      for(i in 1:ncol(site_mu)){
        pred <- site.obs[,grep('y',colnames(site.obs))]
        mod <- lm(site_mu[,i] ~ pred[,i])
        site.rsq[[i]] <- summary(mod)$r.squared
      }
      site.rsq <- unlist(site.rsq)
      names(site.rsq) <- colnames(site_mu)
      
      #wrap rsq output into list and return.
      rsq.out <- list(core.rsq, plot.rsq, site.rsq)
      names(rsq.out) <- c('core','plot','site')
      return(rsq.out)
    }
  
  #clean up all rsq output.----
  core.output <- list()
  plot.output <- list()
  site.output <- list()
  for(i in 1:length(all.output)){
    core.output[[i]] <- all.output[[i]]$core
    plot.output[[i]] <- all.output[[i]]$plot
    site.output[[i]] <- all.output[[i]]$site
  }
  core.output <- do.call(rbind, core.output)
  plot.output <- do.call(rbind, plot.output)
  site.output <- do.call(rbind, site.output)
  output <- list(core.output,plot.output,site.output)
  names(output) <- c('core','plot','site')
  
  #return output.----
  return(output)
} #end function.
