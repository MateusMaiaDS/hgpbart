# Getting simulation function
#' @export
sim_het <- function(n = 100){

  # Generating random samples
  x <- as.matrix(sort(stats::runif(n = n)))
  fx = 4*(x[,1]^2) #quadratric function f
  sx = .2*exp(2*x[,1]) # exponential function s
  y = as.matrix(fx + sx*stats::rnorm(n)) # y = f(x) + s(x) Z
  colnames(x) <- "x"
  colnames(y) <- "y"

  # Returning the list
  return(list(x = x,
              y = y))
}

# Goldberg et. al 1998
#' @export
goldberg_sim <- function(n, seed){

  # Setting a seed
  set.seed(seed)

  # Importing data
  x <- sort(stats::runif(n = n,min = 0,max = 1))
  y <- 2*sin(2*pi*x) + stats::rnorm(n = n,mean = 0,sd = 1.5*x)
  y_true <- 2*sin(2*pi*x)

  return(list( noise_data = data.frame(x = x, y  = y, sd = 1.5*x),
               true_data = data.frame(x = x, y_true = y_true) ) )

}


# Yuhan and and Waba
#' @export
yuhan_waba_sim <- function(n, seed){


  # Setting a seed
  set.seed(seed)

  # Importing data
  x <- sort(stats::runif(n = n,min = 0,max = 1))

  # Sampling the gaussina mean
  gaussian_mean <- 2*exp((-30*(x-0.25)^2+sin(pi*x^2)))-2

  # Gaussian sd
  sd <- exp(sin(2*pi*x))

  y <- stats::rnorm(n = n,
             mean = gaussian_mean,
             sd = sd)

  y_true <- gaussian_mean

  return(list( noise_data = data.frame(x = x, y  = y, sd = sd),
               true_data = data.frame(x = x, y_true = y_true) ) )

}

# Willians 2006
#' @export
willian_sim <- function(n, seed){

  # Setting a seed
  set.seed(seed)

  # Importing data
  x <- sort(stats::runif(n = n,min = 0,max = pi))

  # Retrieving the gaussian mean
  gaussian_mean <- sin(2.5*x)*sin(1.5*x)

  # Retrieving the gaussian sd
  gaussian_sd <- 0.01 + 0.25*(1-sin(2.5*x))^2

  y <- stats::rnorm(n = n,
             mean = gaussian_mean,
             sd = gaussian_sd)

  y_true <- gaussian_mean


  return(list( noise_data = data.frame(x = x, y  = y, sd = gaussian_sd),
               true_data = data.frame(x = x, y_true = y_true) ) )

}

# LIDAR dataset
# lidar_sim <- function(){
#
#   library(SemiPar)
#
#   # Importing the data.
#   lidar_data <- data("lidar")
#
#   x <- lidar$range
#   y <- lidar$logratio
#
#   return(list( noise_hdata = data.frame(x = x, y  = y)) )
#
# }


