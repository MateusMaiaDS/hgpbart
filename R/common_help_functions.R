### All helper functions which are common to bart.R & gpbart.R & some other helper functions
# other common functions related to tree structures are found in

tree_prior <- function(tree, alpha, beta) {

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting internal nodes names
  names_internal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 0))

  # Selecting the depth of the terminal nodes
  depth_terminal <- vapply(tree[names_terminal_nodes], "[[", numeric(1), "depth_node")

  # Selecting the depth of the internal nodes
  depth_internal <- vapply(tree[names_internal_nodes], "[[", numeric(1), "depth_node")

  # Case for stump (No internal node)
  if(length(depth_internal) == 0) {
    log_p <- log(1 - alpha)
  } else {
    # Calculating the log-likelihood
    log_p <- sum(log(1 - alpha * (1 + depth_terminal)^(-beta))) + sum(log(alpha) - beta * log1p(depth_internal))
  }
    return(log_p)
}

# Update tau_j values
update_tau <- function(x,
                       y,
                       a_tau,
                       d_tau,
                       predictions) {

  # Calculating the values of a and d
  n <- nrow(x)

  # Getting the shape parameter from the posterior
  shape_tau_post <- 0.5 * n + a_tau

  # Getting the ratio parameter
  rate_tau_post <- 0.5 * crossprod(y - predictions) + d_tau

  # Updating the \tau
  tau_sample <- stats::rgamma(n = 1, shape = shape_tau_post, rate = rate_tau_post)
    return(tau_sample)
}

# Functions to find the zero for tau
zero_tau_prob <- function(x, naive_tau_value, prob, shape) {

  # Find the zero to the function P(tau < tau_ols) = 0.1, for a defined
  return(stats::pgamma(naive_tau_value,
                shape = shape,
                rate = x) - (1 - prob))
}

zero_tau_prob_squared <- function(rate, naive_tau_value, prob, shape) {

  # Find the zero to the function P(tau < tau_ols) = 0.1, for a defined
  return((stats::pgamma(naive_tau_value,
                 shape = shape,
                 rate = rate) - (1 - prob))^2)
}

# Naive tau_estimation
naive_tau <- function(x, y) {

  # Getting the valus from n and p
  n <- length(y)

  df <- data.frame(x,y)
  colnames(df) <- c(colnames(x),"y")

  # Naive lm_mod
  lm_mod <- stats::lm(formula = y ~ ., data =  df)

  sigma <- stats::sigma(lm_mod)

  tau <- sigma^(-2)

  return(tau)
}

# Naive sigma_estimation
naive_sigma <- function(x,y){

  # Getting the valus from n and p
  n <- length(y)

  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))

  # Naive lm_mod
  lm_mod <- stats::lm(formula = y ~ ., data =  data.frame(y,x))

  sigma <- sqrt(sum((lm_mod$residuals)^2)/(n - p))
    return(sigma)
}

# Return rate parameter from the tau prior
rate_tau <- function(x, # X value
                     y, # Y value
                     prob = 0.9,
                     shape) {
  # Find the tau_ols
  tau_ols <- naive_tau(x = x,
                       y = y)

  # Getting the root
  min_root <-  try(stats::uniroot(f = zero_tau_prob, interval = c(1e-2, 100),
                           naive_tau_value = tau_ols,
                           prob = prob, shape = shape)$root, silent = TRUE)

  if(inherits(min_root, "try-error")) {
    # Verifying the squared version
    min_root <- stats::optim(par = stats::runif(1), fn = zero_tau_prob_squared,
                      method = "L-BFGS-B", lower = 0,
                      naive_tau_value = tau_ols,
                      prob = prob, shape = shape)$par
  }
    return(min_root)
}

# Normalize BART function (Same way as theOdds code)
normalize_bart <- function(y) {

  # Defining the a and b
  a <- min(y)
  b <- max(y)

  # This will normalize y between -0.5 and 0.5
  y  <- (y - a)/(b - a) - 0.5
    return(y)
}

# Now a function to return everything back to the normal scale

unnormalize_bart <- function(z, a, b) {
  # Just getting back to the regular BART
  y <- (b - a) * (z + 0.5) + a
    return(y)
}

# Exporting the RMSE function
#' @export
rmse <- function(obs, pred) {
  return(sqrt(mean((obs - pred)^2)))
}

rMVN_var <- function(mean, Sigma) {
  if(is.matrix(Sigma)) {
    drop(mean + crossprod(PD_chol(Sigma), stats::rnorm(length(mean))))
  } else {
    mean + sqrt(Sigma) * stats::rnorm(length(mean))
  }
}

is_diag_matrix <- function(m) all(m[!diag(nrow(m))] == 0)

PD_chol  <- function(x, ...) tryCatch(chol(x, ...), error=function(e) {
    d    <- nrow(x)
    eigs <- eigen(x, symmetric = TRUE)
    eval <- eigs$values
    evec <- eigs$vectors
      return(chol(x + evec %*% tcrossprod(diag(pmax.int(1e-8, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec), ...))
  }
)

# Calculating CRPS from (https://arxiv.org/pdf/1709.04743.pdf)
#' @export
crps <- function(y,means,sds){

  # scaling the observed y
  z <- (y-means)/sds

  crps_vector <- sds*(z*(2*stats::pnorm(q = z,mean = 0,sd = 1)-1) + 2*stats::dnorm(x = z,mean = 0,sd = 1) - 1/(sqrt(pi)) )

  return(list(CRPS = mean(crps_vector), crps = crps_vector))
}

# Getting e-statistic (Based on Pratola et. al 2018)
#' @export
e_statistic <- function(mean_one,sd_one,
                        mean_two,sd_two,
                        seed = NULL){

  # Setting a seed
  set.seed(seed)

  # Getting the samples and values
  n1 <- length(mean_one)
  n2 <- length(mean_two)
  u_sample <- stats::rnorm(n = n1,mean = mean_one,sd = sd_one)
  v_sample <- stats::rnorm(n = n2, mean = mean_two, sd = sd_two)



  # Calculating the norm of a vector
  norm_vec <- function(x) sqrt(sum(x^2))

  # Creating a distance matrix
  distance_matrix <- as.matrix(stats::dist(rbind(matrix(u_sample),matrix(v_sample))))

  u_dist <- sum(distance_matrix[1:n1,1:n1])
  v_dist <- sum(distance_matrix[(n1+1):sum(n1+n2),(n1+1):sum(n1+n2)])
  u_v_dist <- sum(distance_matrix[1:n1,(n1+1):sum(n1+n2)])

  e_stat <- ((n1*n2)/(n1+n2))*( 2*(1/(n1*n2))*u_v_dist-
                                  ((1/(n1^2))*u_dist)-
                                  (1/(n2^2))*v_dist)


  # Returning the e-statistic
  return(e_stat)

}

# Calculating the
#' @export
#'
hgpbart_qqplot <- function(y, hgp_bart_pred,n_unif=1000,...){

  # Calculating the empirical CDF
  cdf_emp <- function(y_value, y_sample){
    return( mean(y_sample<=y_value))
  }

  # Getting the vector of CDF's for each column for a posterior sample object
  cdf_sample <- function(y,sample_matrix){
    n_obs <- ncol(sample_matrix) # Number of observation
    cdf_vec <- numeric(n_obs)
    # Iterating over all columns
    for(i in 1:n_obs){
      cdf_vec[i] <- cdf_emp(y_value = y[i],
                            y_sample = sample_matrix[,i])
    }
    return(cdf_vec)
  }

  n_post <- nrow(hgp_bart_pred$out$pred)

  # Generating each of the the samples
  posterior_draw <- hgp_bart_pred$out$pred + hgp_bart_pred$out$sd*matrix(stats::rnorm(n = length(hgp_bart_pred$out$pred)),
                                                                        nrow=n_post)

  cdf_y <- cdf_sample(y = y,sample_matrix = posterior_draw)

  # Sampling the unif
  uniform_samples <- stats::runif(n_unif)

  # Plotting the qqplot
  stats::qqplot(cdf_y,uniform_samples,...)
  graphics::abline(0,1,col = "red", lwd = 2)
  return(cdf_y)

}
