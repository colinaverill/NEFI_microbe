#' ddirch_forecast_noLogMap.r
#'
#' @param mod                        #model output from ddirch fitting function.
#' @param cov_mu                     #dataframe of site-level covariates. Must include all covariates used in models.
#' @param cov_sd                     #associated sd's if present.
#' @param glob.covs                  #global level covariates and standard deviations to fill in missing data.
#' @param n.samp                     #number of times to sample covariates/parameters for forecast. Default 1000.
#' @param zero_parameter_uncertainty #turns off drawing from parameter distributions, keeps parameters fixed at means.
#' @param zero_predictor_uncertainty #turns off drawing from covariate distributions, keeps covaraites fixed at means.
#' @param zero_process_uncertainty   #turns off process draw from rdirichlet. Basically just makes predictive interval = credible interval.
#' @param zero_study_variance        #turns off adding variance associated with the study (drawn from normal dist with mean = 0, variance estimated from calibration)
#'
#' @return                           #returns a list of forecasts (mean, 95% CI, 95% PI) same length as model list.
#' @export
#'
#' @examples
source('NEFI_functions/precision_matrix_match.r')
ddirch_forecast_noLogMap <- function(mod, cov_mu, names, cov_sd = NA, n.samp = 1000,
                                     zero_parameter_uncertainty = F,
                                     zero_covariate_uncertainty = F,
                                     zero_process_uncertainty   = F,
                                     zero_study_variance = T){
  #run some tests.----
  if(is.list(mod) == F){
    stop("Your model object isn't a list. It really needs to be.")
  }
  
  #grab the model out of list
  j.mod <- mod$jags_model #separate out the jags model.
  
  #### organize your covariates and covariate standard deviations.----
  preds <- as.character(mod$species_parameter_output$other$predictor)
  #grab covariates that were predictors.
  covs <- cov_mu[,colnames(cov_mu) %in% preds]
  #add intercept.
  covs <- cbind(rep(1,nrow(covs)), covs)
  colnames(covs)[1] <- 'intercept'
  #grab uncertainties in sd, if present.
  if(is.data.frame(cov_sd)){
    cov.sd <- precision_matrix_match(covs, cov_sd)
  }
  
  #report iuf you have a mismatch in covariates vs predictors.----
  check <- preds[!(preds %in% colnames(covs))]
  if(length(check > 0)){
    cat(paste0('Warning: the predictor(s): ',paste(check, sep = ','), ' are absent from the supplied covariate matrix.\n'))
    preds <- preds[preds %in% colnames(covs)]
  }
  check <- colnames(covs)[!(colnames(covs) %in% preds)]
  if(length(check > 0)){
    cat(paste0('Warning: the covariate(s):',paste(check, sep = ','), ' are absent from the model predictor list.\n'))
  }
  
  #re-order to match predictor order. Conditional otherwise it breaks the intercept only case.
  if(ncol(covs) > 1){
    covs   <-   data.frame(covs[,preds])
    if(is.data.frame(cov_sd)){
      cov.sd <- data.frame(cov.sd[,preds]) 
    }
  }
  
  #### Sample parameter and covariate space, make predictions.----
  pred.out <- list()
  cred.out <- list()
  for(j in 1:n.samp){
    
    #Sample parameters from mcmc output.----
    mcmc <- do.call(rbind,j.mod$mcmc)
    mcmc.sample <- mcmc[sample(nrow(mcmc),1),]
    #grab x.m values, convert to matrix.
    x.m <- mcmc.sample[grep("^x\\.m\\[", names(mcmc.sample))]
    x.m <- matrix(x.m, nrow = ncol(covs), ncol = length(x.m)/ncol(covs))
    
    #grab precision of study random effects
    s.tau <-  mcmc.sample[grep("study_tau", names(mcmc.sample))]
    # convert study_tau to SD
    s.eff.sd <- sqrt(1/s.tau)
    
    #if we are fixing parameter uncertainty to zero, do something different.----
    if(zero_parameter_uncertainty == T){
      mcmc.sample <- colMeans(mcmc)
      x.m <- mcmc.sample[grep("^x\\.m\\[", names(mcmc.sample))]
      x.m <- matrix(x.m, nrow = ncol(covs), ncol = length(x.m)/ncol(covs))
    }
    
    #Sample from covariate distributions.----
    #fix covariate uncertainty? Do this by making sd zero.
    if(zero_covariate_uncertainty==T){
      cov_sd <- NA
    }
    
    #only sample if you supplied uncertainties.
    if(is.data.frame(cov_sd)){
      now.cov <- matrix(NA, ncol = ncol(covs), nrow = nrow(covs))
      for(k in 1:ncol(covs)){now.cov[,k] <- rnorm(nrow(covs),covs[,k], cov.sd[,k])}
      colnames(now.cov) <- preds
      now.cov <- data.frame(now.cov)
      #anti-logit relEM if this is one of your covariates.
      if('relEM' %in% colnames(now.cov)){now.cov$relEM <- boot::inv.logit(now.cov$relEM)}
      now.cov <- as.matrix(now.cov)
    }
    
    #If you did not supply covariate uncertainties then now.cov is just covs.
    if(!is.data.frame(cov_sd)){
      now.cov <- as.matrix(covs)
      colnames(now.cov) <- preds
      now.cov <- data.frame(now.cov)
      #anti-logit relEM if this is one of your covariates.
      if('relEM' %in% colnames(now.cov)){now.cov$relEM <- boot::inv.logit(now.cov$relEM)}
      now.cov <- as.matrix(now.cov)
    }
    
    #Combine covariates and parameters to make a prediction.----
    pred.x.m <- matrix(NA, ncol=ncol(x.m), nrow = nrow(covs))
    for(k in 1:ncol(x.m)){
      if (zero_study_variance==FALSE){
        # sample study effects from normal distribution with mean 0 and SD from calibration
        s.eff <- rnorm(1, 0, s.eff.sd)
        pred.x.m[,k] <- exp(now.cov %*% x.m[,k] + s.eff)
      } else {
        pred.x.m[,k] <- exp(now.cov %*% x.m[,k])
      }
    }
    #get mean prediction and then draw from dirichlet distribution.
    cred.out[[j]] <- pred.x.m / rowSums(pred.x.m)
    pred.out[[j]] <- DirichletReg::rdirichlet(nrow(pred.x.m),pred.x.m)
    #Turn off process uncertainty? If so pred.out is just cred.out.
    if(zero_process_uncertainty == T){
      pred.out[[j]] <- cred.out[[j]]
    }
  }
  
  #Summarize prediction mean and confidence intervals.----
  pred.mean       <- apply(simplify2array(cred.out), 1:2, mean)
  pred.ci.0.025   <- apply(simplify2array(cred.out), 1:2, quantile, probs = c(0.025))
  pred.ci.0.975   <- apply(simplify2array(cred.out), 1:2, quantile, probs = c(0.975))
  pred.pi.0.025   <- apply(simplify2array(pred.out), 1:2, quantile, probs = c(0.025))
  pred.pi.0.975   <- apply(simplify2array(pred.out), 1:2, quantile, probs = c(0.975))
  pred.var        <- apply(simplify2array(pred.out), 1:2, var)
  
  #batch into list, column names are group names, row names are siteIDs.
  output <- list(pred.mean, pred.ci.0.025, pred.ci.0.975, pred.pi.0.025, pred.pi.0.975,pred.var)
  for(k in 1:length(output)){
    colnames(output[[k]]) <- names(mod$species_parameter_output)
    rownames(output[[k]]) <- names
  }
  names(output) <- c('mean','ci_0.025','ci_0.975','pi_0.025','pi_0.975','variance')
  
  #return output.----
  return(output)
} #end function.----
