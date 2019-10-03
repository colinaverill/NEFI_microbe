#' Title
#'
#' @param observed   observed values of dependent variable. 
#' @param predicted  predicted values of dependent variable. 
#' @param norm       Normalize RMSE score to the mean abundance of the observed values? Default is FALSE.
#'
#' @return           returns RMSE score.
#' @export
#'
#' @examples
rmse <- function(observed, predicted, norm = F){
  #compatability tests.
  if(sum(is.na(observed )) > 0){stop('Error: NAs in observed vector.' )}
  if(sum(is.na(predicted)) > 0){stop('Error: NAs in predicted vector.')}
  if(length(observed) != length(predicted)){stop('Error: predicted and observed vectors must be the same length.')}
  
  #get sample size.
  N <- length(observed)
  
  #calculate residual sums of squares relative to actual prediction (i.e. "1:1 line")
  SS.res <- sum((observed - predicted)^2)
  
  #divide by sample size, take square root.
  obs.rmse <- sqrt(SS.res / N)
  
  #Normalize RMSE score, if desired.
  if(norm == T){
    obs.rmse <- obs.rmse / mean(observed)
  }
  
  #return RMSE.
  return(obs.rmse)
}