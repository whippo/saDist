#' Quantile of the Zipf distribution
#' 
#' This function generates quantiles from the Zipf distribution via the VGAM pzipf() function.
#' 
#' @param p vector of probabilities.
#' @param N the number of elements.
#' @param s the exponent characterizing the distribution
#' @param logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @return describe output of function
#' @author piklprado
#' @examples 
#' #description of examples here
#' #actual R code to execute examples
#' @importFrom VGAM pzipf
#' 
qzipf <- function(p, N, s, lower.tail = TRUE, log.p = FALSE){
  if (s <= 0)
    stop("s must be greater than zero")
  if (N < 1||!any(is_integer(N)))
    stop("N must be positive integer")
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  ## Ugly: just to make qgs(1, ...) = N
  p[p==1] <- 1+1e-12
  Y <- 1:N
  X <- pzipf(Y, N=N, s=s)
  f1 <- approxfun(x=X, y=Y, method="constant", rule=2)
  f1(p)
}  
#built-in function from 'piklprado' (github.com) to create zipf quantile