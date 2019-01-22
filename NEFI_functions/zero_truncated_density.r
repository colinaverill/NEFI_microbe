#' zero_truncated_density.r
#' For when you want to make a density plot, but your values are bounded to [0,]
#'
#' @param x vector of values to make a density object out ouf
#'
#' @return  returns a zero truncated density object.
#' @export
#'
#' @examples
zero_truncated_density <- function(x){
  h <- density(x)$bw                                                #get dansity bandwith
  w <- 1 / pnorm(0, mean=x, sd=h, lower.tail=FALSE)                 #Compute edge weights.
  suppressWarnings(
    h <- density(x, bw=h, kernel="gaussian", weights=w / length(x)) #Generate truncated distribution.
  )
  h$y[h$x < 0] <- 0                                                 #set negative values to zero.
  #Check: the integral ought to be close to 1:
  if(sum(h$y * diff(h$x)[1]) > 1.1 | sum(h$y * diff(h$x)[1]) < 0.9){
    warning('Zero-truncated density sum is not close to 1. Consider doing something else...')
  }
  return(h)
}