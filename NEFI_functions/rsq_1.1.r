#' rsq_1.1 - generate r-squared for observations relative to 1:1 line
#'
#' @param observed   observed values of dependent variable.
#' @param predicted  predicted values of dependent variable.
#'
#' @return returns r-squared value for observed vs. predicted. Can be negative.
#' @export
#'
#' @examples
rsq_1.1 <- function(observed, predicted){
  if(sum(is.na(observed )) > 0){stop('Error: NAs in observed vector.' )}
  if(sum(is.na(predicted)) > 0){stop('Error: NAs in predicted vector.')}
  #calculate residual sums of squares relative to actual prediction (i.e. "1:1 line")
  SS.res <- sum((observed - predicted)^2)
  #calculate total sum of squares.
  SS.tot <- sum((observed - mean(observed))^2)
  #calculate rsq and return.
  rsq.1 <- 1 - SS.res/SS.tot
  return(rsq.1)
}
