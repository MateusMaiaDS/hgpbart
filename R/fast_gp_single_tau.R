# GP-function main
gp_main <- function(x_train, y_train, x_star, precision_vector, phi, nu, distance_matrix_train, get_cov_star = FALSE) {

  # Getting the distance matrix from x_train and x_star
  distance_matrix_K_star <- distance_matrix(m1 = x_train, m2 = x_star)
  distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star)


  # Calculating the K elements from the covariance structure
  n_train <- nrow(x_train)
  K_y <- kernel_function(squared_distance_matrix = distance_matrix_train,
                         nu = nu,
                         phi = phi) + diag(x = 1/precision_vector, nrow = n_train)
  K_diag <- is_diag_matrix(K_y)
  K_star <- kernel_function(squared_distance_matrix = distance_matrix_K_star,
                            nu = nu, phi = phi)

  # Calculating \alpha
  if(K_diag) {
    L <- diag(K_y)
    alpha <- y_train/L
  } else {
    L <- chol(K_y)
    alpha <- backsolve(L, backsolve(L, y_train, transpose = TRUE, k = n_train), k = n_train)
  }
  mu_star <- crossprod(K_star, alpha)


  # HERE I WILL NOT WORK WITH THE VARIANCE OF EACH TREE
  results <- list(mu_pred = mu_star)


  # ===============#
    return(results)
}

# Function to create the the function K that will be used
# in a Gaussian process (Andrew's Version)
kernel_function <- function(squared_distance_matrix, nu, phi) {

  # Calculating the square matrix
  kernel_matrix <- exp(-squared_distance_matrix / (2 * phi^2)) / nu

  # Case nu = 0
  if(nu == 0 || nu > 1e13){
    kernel_matrix <- matrix(0, nrow = dim(squared_distance_matrix)[1],
                               ncol = dim(squared_distance_matrix)[2])
  }
  # Getting the kernel matrix
    return(kernel_matrix)
}
