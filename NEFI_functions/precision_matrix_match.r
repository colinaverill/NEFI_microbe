#' precision_matrix_match.r
#' takes a a dataframe of x values, and a dataframe of precision
#' Renames columns of precision matrix to match x_mu matrix, so long as some of x label is in precision label (mat to mat_SD)
#' Add columns to precision matrix for x_mu columns not present, sets precision to 10000 (v high).
#' makes sure columns of x_mu and x_precision are in the same order.
#'
#' @param a x_mu matrix as a dataframe with column names.
#' @param b x_precision matrix as a dataframe with column names.
#'
#' @return returns precision matrix as dataframe with columns added for x_mu columns not already in it.
#' @export
#'
#' @examples
#' a <- data.frame(rnorm(100),rnorm(100),rnorm(100))
#' b <- data.frame(rnorm(100),rnorm(100))
#' colnames(a) <- c('intercept','mat','map')
#' colnames(b) <- c('mat_sd','map_sd')
#' precision_matrix_match(a,b)
#' 
precision_matrix_match <- function(a,b){
  #gotta drop any data.table formatting if present.
  #if supplied as matrices need to be dataframes.
  a <- as.data.frame(a)
  b <- as.data.frame(b)
  
  
  #rename columns to match x_mu
  for(i in 1:ncol(a)){
    position <- grep(colnames(a)[i], colnames(b))
    if(length(position) > 0){
      position <- position[1]
      colnames(b)[position] <- colnames(a)[i]
    }
  }
  
  #add columns for x values not present in precision matrix.
  #give these columns a sd value of 0.01, which corresponds to a precision of 10000
  to_add <- colnames(a)[!(colnames(a) %in% colnames(b))]
  new.precision <- matrix(ncol = length(to_add), nrow = nrow(b))
  new.precision[,] <- 0.0001
  colnames(new.precision) <- to_add
  b <- cbind(b,new.precision)
  #order the precision matrix columns to match x_mu matrix.
  b <- b[,colnames(a)]
  b <- as.data.frame(b)
  colnames(b) <- colnames(a)
  
  #return new precision matrix.
  return(b)
  
} ##end function.
