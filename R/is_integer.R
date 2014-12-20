#' Determine if number is integer
#' 
#' This function determines if x is an integer.
#' 
#' @param x vector of numbers
#' @return describe output of function
#' @author Ross Whippo
#' @examples 
#' #description of examples here
#' #actual R code to execute examples
#' 
#built-in function to determine if number is integer (needed for zipf)
is_integer <- function(x, val = .Machine$double.eps^0.5){
  abs(x - round(x)) < val
}