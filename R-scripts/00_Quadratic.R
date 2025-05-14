## Rate Summation Project
## Functions to load for use in other scripts

## Write our own quadratic TPC function (since the one in rTPC is not what we want)


quadratic <- function(temp, tmin, tmax, a){
  est <- -a * (temp - tmin) * (temp - tmax)
  return(est)
}

quadratic.starting_vals <- function(x, y){
  tmin = min(x, na.rm = TRUE)
  tmax = max(x, na.rm = TRUE)
  a = 2 * 10^-4
  return(c(tmin = tmin, tmax = tmax, a = a))
}

quadratic.lower_lims <- function(x, y){
  tmin = -50
  tmax = min(x, na.rm = TRUE)
  a = 0
  return(c(tmin = tmin, tmax = tmax, a = a))
}

quadratic.upper_lims <- function(x, y){
  tmin = max(x, na.rm = TRUE)
  tmax = max(x, na.rm = TRUE) * 10
  a = Inf
  
  return(c(tmin = tmin, tmax = tmax, a = a))
}
