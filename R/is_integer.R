#built-in function to determine if number is integer (needed for zipf)
is_integer <- function(x, val = .Machine$double.eps^0.5){
  abs(x - round(x)) < val
}