#' convert_K.r
#' For potassium (K), but general two step linear converter.
#' Lets say you want to convert a value, A, to a different, value, B.
#' You have data relating A to B.
#' This will take those datasets, and a vector of A values you want to convert.
#' It will return a list, where "mean" is the mean predicted value of B, given A, and "sd" is the predicted standard deviation.
#'
#' @param to_convert vector of values to convert from a->b
#' @param d1         dataframe that relates a to b, a in the first column, b in the second.
#' @param n.sim      number of simulations to generate mean and sd from. 1000 is default and fine.
#'
#' @return           returns a list where mean is the mean predicted value of b, given a, and sd is the predicted standard deviation.
#' @export
#'
#' @examples
#' #phosphorus data converter.
#' a <- runif(100, 0, 1)
#' b <- a*1.2 + rnorm(length(a), sd = 0.1)
#' d1 <- data.frame(a,b)
#' to_convert <- runif(25, 0, 1)
#' convert_K(to_convert,d1)
convert_K <- function(to_convert,d1,n.sim=1000, log10_flip = F){
  #setup predictors for JAGS.----
  y1 <- d1[,2]
  x1 <- d1[,1]
  x1 <- cbind(rep(1,length(x1)),x1)

  #JAGS data object.
  jd <- list(N = length(y1), N.pred = ncol(x1),
             y1=y1, x1=x1)
  
  #Define JAGS model.----
  jags.model = "
  model{
  #covariate parameter priors.
  for(k in 1:N.pred){
  m1[k] ~ dnorm(0,0.001)
  }
  #data model precision priors.
  tau1 <- pow(sigma1, -2)
  sigma1 ~ dunif(0, 100)

  #regression models.
  for(i in 1:N){
  mu1[i] <- inprod(x1[i,], m1)
  y1[i] ~ dnorm(mu1[i], tau1)
  }
  } #end JAGS model loop.
  "
  
  #Fit JAGS model.----
  cat('Fitting JAGS model...\n')
  suppressWarnings(
    jags.out <- run.jags(jags.model,
                         data=jd,
                         n.chains=3,
                         adapt = 500,
                         burnin = 1000,
                         sample = 1000,
                         monitor=c('m1'),
                         method = 'rjags')
  )
  jout <- summary(jags.out)
  
  #grab mcmc of model parameters.
  mcmc <- do.call(rbind, jags.out$mcmc)
  
  #convert to_convert data from a to b, propagating uncertainty!----
  convert.out <- list()
  for(i in 1:n.sim){
    par <- mcmc[sample(nrow(mcmc),1),]
    b <- par[2]*to_convert + par[1]
    convert.out[[i]] <- b
  }
  
  #wrap up output to return.----
  convert.out <- do.call(cbind, convert.out)
  #log10 variable?
  suppressWarnings(
    if(log10_flip == T){
      convert.out <- log10(convert.out)
    }
  )
  mu <- rowMeans(convert.out, na.rm = T)
  mu.sd <- apply(convert.out, 1, sd, na.rm = T)
  output <- list(mu, mu.sd)
  names(output) <- c('mean','sd')
  
  #return output, end function.----
  return(output)
} #end function.
