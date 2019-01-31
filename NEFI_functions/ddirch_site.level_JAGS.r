#' site.level_dirlichet_jags()
#' Fits a site.level dirlichet model with x-values ONLY OBSERVED AT SITE LEVEL.
#' Accounts for missing data.
#'
#' @param y species matrix as a dataframe. Column names are species names.
#' @param x_mu site level observation means as a dataframe. column names are predictor names.
#' @param x_sd site level observation sd as a dataframe. column names should match site_mu. optional.
#' @param adapt   number of adaptive iterations for JAGS simulation.
#' @param burnin  number of burnin iterations for JAGS simulation.
#' @param sample  number of sample iterations for JAGS simulation.
#' @param n.chains number of chains for JAGS simulation.
#' @param parallel whether or not to run JAGS chains in parallel.
#'
#' @return returns a list with the fitted model, a list of species by parameter tables, and matrices of predicted, observed and residual values on the observation scale (0,1)
#' @export
#'
#' @examples
site.level_dirlichet_jags     <- function(y,
                                          x_mu, 
                                          x_sd = NA,
                                          adapt = 500, burnin = 1000, sample = 2000, n.chains = 3, 
                                          parallel = F, silent.jags = F, parallel_method = 'rjparallel'){
  #Load some important dependencies.
  source('NEFI_functions/crib_fun.r')
  source('NEFI_functions/sd_to_precision.r')
  source('NEFI_functions/precision_matrix_match.r')
  
  #Some checks before we get started.
  y    <- as.data.frame(y)
  x_mu <- as.data.frame(x_mu)
  #core level. First column needs to be intercept, a vector of 1s.
  if(mean(x_mu[,1]) != 1){stop('First column in x_mu is not a vector of 1s for the intercept. This needs to be or this function doesnt work.')}
  if(ncol(x_mu) < 2){stop('We need an intercept and at least one predictor for this to work. x_mu has less than 2 columns. Try again buddy.')}
  
  #grab names
  y.names <- colnames(y)
  x.names <- colnames(x_mu)

  ###massage your data together.
  #deal with zero relative abundances.
  y <- data.frame(y)
  #y <- data.frame(lapply(y, crib_fun))
  
  #make sd objects if they were not supplied.
  if(is.na(x_sd)){x_sd = data.frame(rep(1,nrow(x_mu)))}

  #Match up predictors and their SD. if no SD supplied we assign ~perfect precision.
  x_sd <- precision_matrix_match(x_mu,x_sd)

  #covert sd to precision. output is matrix.
  x_precision <- sd_to_precision(x_sd)

  #make sure every else is a matrix.
  y <- as.matrix(y)
  x_mu <- as.matrix(x_mu)
  
  ###setup jags data object.----
  jags.data <- list(N = nrow(y), N.spp = ncol(y), #number of observations and number of species
                    N.preds = ncol(x_mu),    #number of x predictors
                    x_mu = x_mu,                  #x-value mean matrix
                    x_precision = x_precision,    #x-value precision matrix
                    y = y)                        #species matrix, y
  
  ###specify JAGS model.----
  jags.model = "
  model {
  #parameter priors for each species.
  alpha ~ dnorm(0, 1.0E-3) 
  for(i in 1:N.preds){
    x.mm[i,1] <- 0
    for (j in 2:N.spp) {
      x.mm[i,j] ~ dnorm(0, 1.0E-3)
    }
  }
  
  ### Begin missing data model ###
  #missing X data priors, site-level.
  for(m in 1:N.preds){
  x.global[m] ~ dnorm(0,1.0E-4) #global level parameter prior.
     x.tau[m] ~ dgamma(0.1,0.1)
  }
  
  #fill in any missing X values at site level.
  for(m in 1:N.preds){
  for(i in 1:N){   x_mu[i,m] ~ dnorm(x.global[m], x.tau[m])}   ## estimate global mean
  }
  
  #predictor (x) values drawn from distributions.
  for(j in 1:N.preds){for(i in 1:N){x[i,j] ~ dnorm(x_mu[i,j], x_precision[i,j])}} #x values

  #mean center all predictors (except intercept).
  for(i in 1:N){
    x.center[i,1] <- 1
    for(j in 2:N.preds){
      x.center[i,j] <- x[i,j] - mean(x[,j])
    }
  }

  #save mean values for back transforming intercept values.
  for(j in 1:N.preds){x.center.save[j] <- mean(x[,j])}
  
  #fit species abundances as a linear combination of predictors and parameters.
  for(i in 1:N){
    for(j in 1:N.spp){
      log(a0[i,j]) <- alpha + inprod(x.mm[,j], x.center[i,])
    }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
  }
  
  #map to original parameterization, assuming first column of predictors is intercept.
  for (j in 1:N.spp) {
    x.m[1,j] <- alpha + x.mm[1,j] - inprod(x.mm[2:N.preds,j], x.center.save[2:N.preds])
    for (i in 2:N.preds){
      x.m[i,j] <- x.mm[i,j]
    }
  }
  
  } #close model loop.
  "
  
  ###Fit JAGS model.----
  #parallel or not parallel.
  run.method <- 'rjags'
  if(parallel == T){
    run.method = parallel_method
  }
  #run jags model.
  jags.out <- runjags::run.jags(   model = jags.model,
                                   data = jags.data,
                                   adapt = adapt,
                                   burnin = burnin,
                                   sample = sample,
                                   n.chains = n.chains,
                                   method = run.method,
                                   silent.jags = silent.jags,
                                   monitor = c('x.m','x.mm','alpha','deviance'))
  #summarize output
  out <- summary(jags.out)
  
  #grab parmeters by species, make a list of species-parameter dataframes----
  output.list <- list()
  for(i in 1:ncol(y)){
    z <- out[grep("^x\\.m\\[",rownames(out), value = T),]
    z <-   z[grep(paste0(',',i,']'),rownames(z)),]
    names <- c(x.names)
    z <- cbind(names,z)
    colnames(z)[1] <- 'predictor'
    output.list[[i]] <- data.frame(z)
  }
  names(output.list) <- y.names
  
  #get x.mm and alpha means for calculating pD later.
  x.mm <- out[grep('x.mm', rownames(out)),4]
  dim(x.mm) <- c(length(x.names), ncol(y))
  alpha <- out[grep('alpha', rownames(out)), 4]
  
  #get the matrix of predicted y values.
  super.x <- x_mu
  
  pred.list <- list()
  for(i in 1:ncol(y)){
    pred <- exp(as.matrix(super.x) %*% as.numeric(as.character(output.list[[i]][,5])))
    pred.list[[i]] <- pred
  }
  predicted <- do.call('cbind',pred.list)
  predicted <- predicted / rowSums(predicted)
  colnames(predicted) <- colnames(y)
  #get matrix of residuals
  resid <- y - predicted
  
  #get deviance score.
  deviance <- out[grep('deviance',rownames(out)),]
  
  #make a super output that also returns model.-----
  super.list <- list(jags.out, output.list,predicted,y,resid,deviance,x.mm,alpha)
  names(super.list) <- c('jags_model','species_parameter_output','predicted','observed','residual','deviance','x.mm','alpha')
  
  #return model and output
  return(super.list)
  
}