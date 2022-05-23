## GP-Bart
#' @useDynLib hgpbart
#' @importFrom Rcpp sourceCpp
# ==================================#
# Objects to test the tree_complete_conditional function
# ==================================#

# Function to calculate the tree complete conditional using het-GP-BART model (FOR \tau)
tree_complete_conditional_het_bart_tau <- function(x,
                                                   precision_sq_residuals_values,
                                                   tree,
                                                   a_tau,
                                                   d_tau) {

  # Getting the number of observations of the data
  n <- nrow(x)

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))

  # Retrieving Residuals terminal nodes
  sum_precision_sq_residuals_terminal_nodes <- unlist(lapply(terminal_nodes, function(x) {
    sum(precision_sq_residuals_values[x$observations_index])
  }))


  # Retrieve all nodes values and calculate all of them
  log_posterior <- lgamma(0.5*nodes_size+a_tau) - (0.5*nodes_size+a_tau)*log(0.5*sum_precision_sq_residuals_terminal_nodes+d_tau)

  return(log_posterior)
}

# Function to calculate the tree complete conditional using het-BART model (FOR \tau)
update_tau_het_bart <- function(x,
                                precision_sq_residuals_values,
                                tree,
                                a_tau,
                                d_tau) {


  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Getting the number of observations of the data
  n <- nrow(x)

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))

  # Retrieving Residuals terminal nodes
  sum_precision_sq_residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    sum(precision_sq_residuals_values[x$observations_index])
  })

  # Sampling tau values
  tau_sample <- mapply(nodes_size, sum_precision_sq_residuals_terminal_nodes,
                       FUN = function(nodes_size, sum_res_sq){
                         stats::rgamma(n = 1,shape = nodes_size*0.5+a_tau,
                                rate = 0.5*sum_res_sq+d_tau)
                       })

  # Adding the mu values calculated
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$tau <- tau_sample[[names_terminal_nodes[i]]]
  }

  return(tree)
}

tree_complete_conditional_het_gpbart_mu <- function(tree, x,
                                                    residuals,
                                                    nu = 1,
                                                    phi = 1,
                                                    tau_mu,
                                                    precision_vector,
                                                    number_trees_mu = number_trees_mu) {

  # Getting the number of observations of data
  n <- nrow(x)

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Retrieving Precision Residuals terminal nodes
  precision_terminal_nodes <- lapply(terminal_nodes, function(x) {
    precision_vector[x$observations_index]
  })


  # Calculating value S
  S <- unlist(mapply(terminal_nodes,
                     FUN=function(z, x=z$Omega_plus_P_inv) {
    sum(x)
  }, SIMPLIFY = TRUE)) + tau_mu

  # Log Omega
  log_det_Omega_plus_P_inv <- vapply(terminal_nodes, function(z, x=z$Omega_plus_P_inv) {
    if(z$is_Omega_diag) sum(log(diag(x))) else determinant(x, logarithm = TRUE)$modulus
  }, numeric(1))

  # Defining RT_Omega_I_R
  RTR <- unlist(mapply(terminal_nodes, residuals_terminal_nodes,
                       FUN = function(nodes, resid, x=nodes$Omega_plus_P_inv) {
    if(nodes$is_Omega_diag) sum(resid^2 * diag(x)) else crossprod(resid, crossprod(x, resid))
  }, SIMPLIFY = FALSE))

  # The term R^{T} solve(Omega + I ) 1
  R_Omega_I_one <- unlist(mapply(terminal_nodes, residuals_terminal_nodes,
                                 FUN = function(nodes, residuals, x=nodes$Omega_plus_P_inv) {
    if(nodes$is_Omega_diag) sum(residuals * diag(x)) else rowSums(crossprod(residuals, x))
  }, SIMPLIFY = FALSE))

  log_posterior <- 0.5 * sum(log_det_Omega_plus_P_inv) -
    0.5 * sum(log(S) - log(tau_mu)) - 0.5 * sum(RTR) + 0.5 * sum((R_Omega_I_one^2) / S)

    return(list(log_posterior = log_posterior,
                S = S,
                RTR = RTR,
                R_Omega_I_one = R_Omega_I_one))
}

# Generate mu_j values
update_mu <- function(tree,
                      x,
                      residuals,
                      likelihood_object,
                      seed = NULL) {

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Mu mean value
  mu_mean <- likelihood_object$R_Omega_I_one / likelihood_object$S

  # Calculating mu values
  mu <- mapply(mu_mean, likelihood_object$S,
               FUN = function(x, y) {
                 stats::rnorm(
                   n = 1,
                   mean = x,
                   sd = sqrt(1/y)
                 )
               }, SIMPLIFY = FALSE)

  # Adding the mu values calculated
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$mu <- mu[[names_terminal_nodes[i]]]
  }
    return(tree)
}

update_residuals <- function(tree, x, nu, phi, residuals, seed = NULL) {
  # set.seed(seed)

  # New g (new vector prediction for g)
  residuals_new <- rep(NA, length(residuals))

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Getting the \mu_{j} vector
  mu_values <- vapply(terminal_nodes, "[[", numeric(1), "mu")

  # Calculating Omega matrix
  Omega_matrix <- lapply(terminal_nodes, "[[", "Omega_matrix")

  # Checking if diagonal
  is_Omega_diag <- lapply(terminal_nodes, "[[", "is_Omega_diag")

  # Calculating Omega matrix plus I INVERSE
  Omega_matrix_plus_P_inv  <- lapply(terminal_nodes, "[[", "Omega_plus_P_inv")

  # Calculating g_mean posterior
  residuals_mean <- mapply(Omega_matrix,
                           Omega_matrix_plus_P_inv,
                           residuals_terminal_nodes,
                           mu_values,
                           is_Omega_diag,
                           FUN = function(omega, omega_p_inv, residuals, mu, omega_is_diag) {
                             # Getting the values
                             if(omega_is_diag) {
                               mu + diag(omega) * diag(omega_p_inv) * (residuals - mu)
                             } else {
                               mu + crossprod(omega, crossprod(omega_p_inv, residuals - mu))
                             }
                           }, SIMPLIFY = FALSE)

  # Calculating g_mean posterior
  residuals_variance <- mapply(Omega_matrix,
                               Omega_matrix_plus_P_inv,
                               residuals_terminal_nodes,
                               mu_values,
                               is_Omega_diag,
                               FUN = function(omega, omega_p_inv, residuals, mu, omega_is_diag) {
                                 # Getting the Omega value
                                 if(omega_is_diag) {
                                   diag(omega) - diag(omega)^2 * diag(omega_p_inv)
                                 } else {
                                   omega - crossprod(omega, crossprod(omega_p_inv, omega))
                                 }
                               }, SIMPLIFY = FALSE)

  # Calculating g_mean posterior
  residuals_sample <- mapply(FUN=rMVN_var, residuals_mean, residuals_variance, SIMPLIFY = FALSE)

  # Adding the mu values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving g
    residuals_new[terminal_nodes[[i]]$observations_index] <- residuals_sample[[i]]
  }
    return(residuals_new)
}


# ==============#
# Bart-hGP FUNCTION
#' @export
hgp_bart <- function(x, y,
                     number_trees_mu = 2, # Setting the number of trees mu
                     number_trees_tau = 2, # Setting the number of trees for tau
                     prior = list(mu = 0,
                        alpha_mu = 0.5, # Alpha from prior
                        beta_mu = 10, # Beta from prior
                        alpha_tau = 0.95, # Alpha from prior
                        beta_tau = 2, # Beta from prior

                        tau = 1, # Tau from prior,

                        a_tau = 10, # Prior from a_v_ratio gamma
                        d_tau = 3, # Prior from d_v_ratio gamma,
                        K_bart = 2,
                        prob_tau = 0.9,
                        kappa = 0.5),
                     mcmc = list(n_iter = 2000, # Number of iterations
                        burn = 500, # Number of burn
                        bart_boolean = TRUE,
                        bart_number_iter = 250), #
                     control = list( # Controlling and defining other parameters
                        node_min_size_mu = 15, # Min node size,
                        node_min_size_tau = 5,
                        scale_boolean = TRUE,
                        x_scale =  TRUE,
                        seed = NULL, # Setting a fixed seed
                        rotation = TRUE, # If rotated lon. and lat. will be used in tree building
                        theta = NULL, # If theta is NULL, then the rotation angle will be randomly selected
                        discrete_phi_boolean = FALSE,
                        gp_variables = colnames(x))   # Selecting the GP-Variables

                     ) {



  # Setting the parameters from the lists to actual objects

  # Prior
  mu <- prior$mu
  alpha_mu <- prior$alpha_mu
  beta_mu <- prior$beta_mu
  alpha_tau <- prior$alpha_tau
  beta_tau <- prior$beta_tau
  tau <- prior$tau
  a_tau <- prior$a_tau
  d_tau <- prior$d_tau
  K_bart <- prior$K_bart
  prob_tau <- prior$prob_tau
  kappa <- prior$kappa

  # MCMC
  n_iter <- mcmc$n_iter
  burn <- mcmc$burn
  bart_boolean <- mcmc$bart_boolean
  bart_number_iter <- mcmc$bart_number_iter

  # Control -- other parameters
  node_min_size_mu <- control$node_min_size_mu
  node_min_size_tau <- control$node_min_size_tau

  scale_boolean <- control$scale_boolean
  x_scale <- control$x_scale
  seed <- control$seed
  rotation <- control$rotation
  theta <- control$theta
  discrete_phi_boolean <- control$discrete_phi_boolean
  gp_variables <- control$gp_variables


  # Changing the node_min_size
  if(node_min_size_mu>=nrow(x)){
    stop("Node min size is greater than the number of observation")
  }

  # Error of the matrix
  if(is.null(colnames(x))) {
    stop("Insert a valid NAMED matrix")
  }

  if(node_min_size_mu == 0) {
    stop("Node Minimum Size need to be greater than 0")
  }



  # Recommendation about the min_node_size
  if(node_min_size_mu < 15) {
    warning("\n It is recommended that the node_min_size should be of at least 15 observations.", immediate.=TRUE)
  }




  # This will be defining the nu the default value
  nu_vector = NULL

  # This parameter is a "scale paramter" to the GP
  phi_vector = rep(0.1 / (sqrt(number_trees_mu)), number_trees_mu)


  # Creating the prediction elements to be stored
  predictions = matrix(0, nrow = number_trees_mu, ncol = nrow(x))
  predictions_list =  NULL

  # If there's only one covariate
  rotation <- ncol(x) != 1

  # Scaling the x
  if(x_scale){
    x_original <- x

    # Scaled version
    xscale <- scale(x)
    mean_x <- attr(xscale,"scaled:center")
    sd_x <- attr(xscale,"scaled:scale")
    x <- as.matrix(xscale)
  }

  # Adjusting the kappa (avoiding the Infinity error)
  if(kappa == 1 ){
    kappa <- kappa - 2*.Machine$double.eps
  }

  if(kappa == 0 ){
    kappa <- kappa + 2*.Machine$double.eps
  }

  # Getting the maximum and minimum values from a distance matrix
  distance_matrix_x <- symm_distance_matrix(m1 = x[,gp_variables, drop = FALSE])
  distance_range <- range(distance_matrix_x[upper.tri(distance_matrix_x)])
  distance_min <- sqrt(distance_range[1])
  distance_max <- sqrt(distance_range[2])

  # Setting seed
  set.seed(seed)
  acc_ratio <- 0
  acc_ratio_phi <- 0
  acc_ratio_nu <- 0

  # Saving a_min and b_max
  a_min <- NULL
  b_max <- NULL

  # Scale values
  if(scale_boolean) {
    # Normalizing y
    a_min <- min(y)
    b_max <- max(y)

    y_scale <- normalize_bart(y = y)

    # Defing the nu vector if not in default
    if(is.null(nu_vector)) {
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4 * number_trees_mu * K_bart^2)/(1 - kappa), number_trees_mu)
    } else if(length(nu_vector) == 1) {
      nu_vector <- rep(nu_vector, number_trees_mu)
    }

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu_bart <- (4 * number_trees_mu * K_bart^2)
    tau_mu_gpbart <- tau_mu_bart/kappa

    # Saving a_tau original
    a_tau_original <- a_tau
    a_tau <- a_tau_original^(1/number_trees_tau)

    # Getting the optimal tau values
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau_original)^(1/number_trees_tau)

  } else {

    # Not scaling the y
    y_scale <- y

    if(is.null(nu_vector)) {
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep(4 * number_trees_mu * K_bart^2/((1 - kappa) * (max(y_scale) - min(y_scale))^2), number_trees_mu)
    } else if(length(nu_vector) == 1) {
      nu_vector <- rep(nu_vector, number_trees_mu)
    }

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu_bart <- (4 * number_trees_mu * K_bart^2)/((max(y_scale) - min(y_scale))^2)
    tau_mu_gpbart <- tau_mu_bart/kappa

    # Saving a_tau original
    a_tau_original <- a_tau
    a_tau <- a_tau_original^(1/number_trees_tau)

    # Getting the optimal tau values
    d_tau <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau_original)^(1/number_trees_tau)
  }

  # Getting the number of observations
  n <- length(y)

  # Creating the stor of accepted tree verbs and which split it was
  verb_store_list_mu <- list()
  verb_store_list_tau <- list()


  # Creating the mu and tau list of store
  tree_mu_store <- list()
  tree_tau_store <- list()

  # Getting the current trees
  current_trees_mu <- list()
  current_trees_tau <- list()

  current_trees_proposal_mu <- list()
  current_trees_proposal_tau <- list()


  # Creating the current_partial_residuals
  current_partial_residuals_matrix <-
  current_predictions_matrix <- matrix(NA, nrow = number_trees_mu, ncol = nrow(x))

  # Creating the precision matrix
  precisions <- matrix(1,
                       ncol = nrow(x),
                       nrow = number_trees_tau)

  # Creating the predictions saving
  current_partial_residuals_list <- list()
  current_predictions_list <- list()
  current_precision_list <- list()


  # Storage containers
  store_size <- (n_iter - burn)

  y_hat_store <-
  y_precision_hat_store <-
  y_hat_store_proposal <- matrix(NA, ncol = length(y), nrow = store_size)

  # Saving full conditional stores
  full_cond_store <-
  phi_store <-
  phi_proposal_store <- matrix(NA, ncol = number_trees_mu, nrow = store_size)
  phi_vector_proposal <- rep(0.1, number_trees_mu)

  # Creating the list of trees stumps
  for(i in seq_len(number_trees_mu)) {
    # Creating the fixed two split trees
    current_trees_mu[[i]] <- stump_mu(x = x, mu = mu)
    current_trees_proposal_mu[[i]] <- stump_mu(x = x, mu =  mu)
  }

  # Creating the list of trees stumps
  for(i in seq_len(number_trees_tau)) {
    # Creating the fixed two split trees
    current_trees_tau[[i]] <- stump_tau(x = x, tau = tau)
    current_trees_proposal_tau[[i]] <- stump_tau(x = x, tau = tau)
  }

  # Giving the names of trees mu
  names(current_trees_mu) <-
  names(current_trees_proposal_mu) <- vapply(seq_len(number_trees_mu), function(x) paste0("tree_", x), character(1)) # Naming each tree

  # Giving the names of trees tau
  names(current_trees_tau) <-
    names(current_trees_proposal_tau) <- vapply(seq_len(number_trees_tau), function(x) paste0("tree_", x), character(1)) # Naming each tree


  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running GP-Sum-Sampler..."
  )

  # Getting the parameters to unormalize the data
  a <- min(y)
  b <- max(y)

  # Having the partial residuals matrix and precisions
  current_partial_residuals_matrix <- predictions
  current_partial_precision_matrix <- precisions

  # Setting initial values for phi vector
  for(i in seq_len(n_iter)) {

    utils::setTxtProgressBar(progress_bar, i)

    # (FIX THIS LATTER) Changing the bart boolean, when reach the maximum
    if(i >= bart_number_iter){
      bart_boolean <- FALSE
      tau_mu <- tau_mu_gpbart
    } else tau_mu <- tau_mu_bart

    if((i > burn)) {

      # Saving the store of the other ones
      curr <- (i - burn)
      tree_mu_store[[curr]] <- current_trees_mu
      tree_tau_store[[curr]] <- current_trees_tau

      y_hat_store[curr, ] <- colSums(predictions)
      y_precision_hat_store[curr, ] <- apply(precisions,2,function(x){exp(sum(log(x)))})


      # Saving the current partial
      current_partial_residuals_list[[curr]] <- current_partial_residuals_matrix

      # Saving the predictions
      current_predictions_list[[curr]] <- predictions
      current_precision_list[[curr]] <- precisions

      phi_store[curr, ] <- phi_vector
      verb_store_list_mu[[curr]] <- verb_store_mu
      verb_store_list_tau[[curr]] <- verb_store_tau
    }

    # Creating a boolean to create the first trees only using BART model
    if(bart_boolean) {

      # Verb Store
      verb_store_mu <- data.frame(verb = rep(NA, number_trees_mu),
                               accepted = rep(NA, number_trees_mu),
                               identical = rep(NA, number_trees_mu))

      for(j in seq_len(number_trees_mu)) {

        # Getting the verb list

        # Calculating the residuals for each tree
        if(number_trees_mu > 1) {

          # Calculating what Chipman called as R(j) = y - g_others_trees
          if(number_trees_mu > 2) {
            # current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }

        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow", "grow_projection", "prune", "change", "change_projection", "swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change", "swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }

        # Case of rotation
        if(rotation){
          if(i < max(floor(0.1 * burn), 10) || length(current_trees_mu[[j]]) == 1) verb <- sample(c("grow", "grow_projection"),
                                                                                               size = 1) # Grow the tree for the first few iterations
        } else {
          if(i < max(floor(0.1 * burn), 10) || length(current_trees_mu[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }

        new_trees_mu <- current_trees_mu # Creating new trees to updated as candidate

        new_trees_mu[[j]] <- update_tree_verb(
          tree = current_trees_mu[[j]],
          x = x,
          node_min_size = node_min_size_mu,
          verb = verb
        )

        # Calculating the likelihood of the new tree
        likelihood_new <- tree_complete_conditional_het_bart_mu(
          tree = new_trees_mu[[j]], # Calculate the full conditional
          residuals_values = current_partial_residuals,
          x = x,tau_mu = tau_mu,
          precision_vector = apply(precisions,
                                   2,
                                   function(x){exp(sum(log(x)))})
        )

        # Calculating the likelihood of the old tree
        likelihood_old <- tree_complete_conditional_het_bart_mu(
          tree = current_trees_mu[[j]], # Calculate the full conditional
          residuals_values = current_partial_residuals,
          x = x,tau_mu = tau_mu,
          precision_vector = apply(precisions,
                                   2,
                                   function(x){exp(sum(log(x)))})
        )

        # Extracting only the likelihood
        l_new <- sum(likelihood_new) +
          tree_prior(
            tree = new_trees_mu[[j]], # Calculate the tree prior
            alpha = alpha_mu,
            beta = beta_mu
          )

        # Extracting only the likelihood
        l_old <- sum(likelihood_old) +
          tree_prior(
            tree = current_trees_mu[[j]], # Calculate the tree prior
            alpha = alpha_mu,
            beta = beta_mu
          )

        # Getting the log of transitin prob
        log_transition <- log_transition_prob(current_tree = current_trees_mu[[j]],
                                              new_tree = new_trees_mu[[j]],verb = verb)

        # (log) Probability of accept the new proposed tree
        acceptance <- (l_new - l_old + log_transition)

        # In case of acceptance
        if(acceptance > 0 || acceptance > -stats::rexp(1)) {

          # Counting acc ratio
          acc_ratio <- acc_ratio + 1

          # Make changes if accept
          current_trees_mu <- new_trees_mu

          # Create a data.frame with the verb and if it was accepted or not
          verb_store_mu[j,"verb"] <- verb
          verb_store_mu[j,"accepted"] <- TRUE

        } else {

          # Create a data.frame with the verb and if it was accepted or not
          verb_store_mu[j,"verb"] <- verb
          verb_store_mu[j,"accepted"] <- FALSE

        } # End of accept for MH sample

        # # # To update the mu values
        current_trees_mu[[j]] <- update_mu_het_bart(
          tree = current_trees_mu[[j]],
          x = x,
          residuals = current_partial_residuals,
          tau_mu = tau_mu,
          precision_vector = apply(precisions,
                                   2,
                                   function(x){exp(sum(log(x)))})
        )

        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar
        predictions[j,] <- get_mu_bart(
          tree = current_trees_mu[[j]], x = x
        )
        current_partial_residuals_matrix[j,] <- current_partial_residuals

      } # End of iterations over mu trees

      # ========
      # ITERATING NOW OVER THE \TAU TREES
      # ========
      # Iterating over the trees
      for(j in seq_len(number_trees_tau)) {

        # Verb Store
        verb_store_tau <- data.frame(verb = rep(NA, number_trees_tau),
                                    accepted = rep(NA, number_trees_tau),
                                    identical = rep(NA, number_trees_tau))

        # Making a exception for the case of only one tree
        if(number_trees_tau == 1) {

          # Getting the current partial values
          current_partial_residuals_precision_squared <- matrix((y_scale^2)/c(precisions), ncol = length(y_scale))
          # MAYBE NEED TO REVIEW THIS LINE

        } else {

          # Getting the current partial values
          current_partial_residual_precision_squared <- ((y_scale - colSums(predictions[,,drop = FALSE]))^2)*apply(precisions[-j,,drop = FALSE],2, function(x){exp(sum(log(x)))})
        }

        # Propose a new tree based on the verbs: grow/prune/change/swap
        verb <- sample(c("grow", "prune", "change", "swap"),
                       prob = c(0.25,0.25,0.4,0.1), size = 1)


        # Case of rotation
        if (i < max(floor(0.1 * burn), 10) || length(current_trees_tau[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations


        # GETTING A NEW TREE
        new_trees_tau <- current_trees_tau # Creating new trees to updated as candidate

        new_trees_tau[[j]] <- update_tree_verb(
          tree = current_trees_tau[[j]],
          x = x,
          node_min_size = node_min_size_tau,
          verb = verb
        )

        # Calculating the likelihood of the new tree
        likelihood_new <- tree_complete_conditional_het_bart_tau(
          tree = new_trees_tau[[j]], # Calculate the full conditional
          precision_sq_residuals_values = current_partial_residual_precision_squared,
          x = x,a_tau = a_tau, d_tau = d_tau
        )

        # Calculating the likelihood of the old tree
        likelihood_old <- tree_complete_conditional_het_bart_tau(
          tree = current_trees_tau[[j]], # Calculate the full conditional
          precision_sq_residuals_values = current_partial_residual_precision_squared,
          x = x,a_tau = a_tau, d_tau = d_tau
        )

        # Extracting only the likelihood
        l_new <- sum(likelihood_new) +
          tree_prior(
            tree = new_trees_tau[[j]], # Calculate the tree prior
            alpha = alpha_tau,
            beta = beta_tau
          )

        # Extracting only the likelihood
        l_old <- sum(likelihood_old) +
          tree_prior(
            tree = current_trees_tau[[j]], # Calculate the tree prior
            alpha = alpha_tau,
            beta = beta_tau
          )

        # Getting the log of transitin prob
        log_transition <- log_transition_prob(current_tree = current_trees_tau[[j]],
                                              new_tree = new_trees_tau[[j]],verb = verb)

        # (log) Probability of accept the new proposed tree
        acceptance <- (l_new - l_old + log_transition)

        if(acceptance > 0 || acceptance > -stats::rexp(1)) {

          # Counting acc ratio
          acc_ratio <- acc_ratio + 1

          # Make changes if accept
          current_trees_tau <- new_trees_tau

          # Create a data.frame with the verb and if it was accepted or not
          verb_store_tau[j,"verb"] <- verb
          verb_store_tau[j,"accepted"] <- TRUE

        } else {

          # Create a data.frame with the verb and if it was accepted or not
          verb_store_tau[j,"verb"] <- verb
          verb_store_tau[j,"accepted"] <- FALSE

        } # End of accept for MH sample

        # # # To update the mu values
        current_trees_tau[[j]] <- update_tau_het_bart(
          tree = current_trees_tau[[j]],
          x = x,
          precision_sq_residuals_values = current_partial_residual_precision_squared,
          a_tau = a_tau,d_tau = d_tau)

        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar
        precisions[j,] <- get_tau_bart(
          tree = current_trees_tau[[j]], x = x
        )
        current_partial_precision_matrix[j,] <- current_partial_residual_precision_squared

      }

      # =================
      # ATTENTION HERE!!!
      # =================
    } else { # Going over the case where the BART-boolean is no more valid

      # Verb Store
      verb_store_mu <- data.frame(verb = rep(NA,number_trees_mu),
                               accepted = rep(NA,number_trees_mu),
                               identical = rep(NA,number_trees_mu))

      for(j in seq_len(number_trees_mu)) {

        # Calculating the residuals for each tree
        if(number_trees_mu > 1) {

          # Calculating what Chipman called as R(j) = y - g_others_trees
          if(number_trees_mu > 2) {
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }

        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow", "grow_projection", "prune", "change", "change_projection", "swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change","swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }

        # Case of rotation
        if(rotation){
          if(i < max(floor(0.1 * burn), 10) | length(current_trees_mu[[j]]) == 1) verb <- sample(c("grow","grow_projection"),
                                                                                              size = 1) # Grow the tree for the first few iterations
        } else {
          if(i < max(floor(0.1 * burn), 10) || length(current_trees_mu[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }

        # GETTING A NEW TREE
        new_trees_mu <- current_trees_mu # Creating new trees to updated as candidate

        new_trees_mu[[j]] <- update_tree_verb(
          tree = current_trees_mu[[j]],
          x = x,
          node_min_size = node_min_size_mu,
          verb = verb, rotation = rotation, theta = theta
        )

        # ==================== #
        # Getting the Omega Inverse the current and the future tree
        # ==================== #

        # Getting the inverse for the current terminal nodes
        current_trees_mu[[j]] <- inverse_omega_plus_P(tree = current_trees_mu[[j]],
                                                   x = x,
                                                   precision_vector = apply(precisions,
                                                                            2,
                                                                            function(x){exp(sum(log(x)))}),
                                                   nu = nu_vector[j],
                                                   phi = phi_vector[j])

        # Getting the inverse for the new tree terminal nodes
        new_trees_mu[[j]] <- inverse_omega_plus_P(tree = new_trees_mu[[j]],
                                               x = x,
                                               precision_vector = apply(precisions,
                                                                        2,
                                                                        function(x){exp(sum(log(x)))}),
                                               nu = nu_vector[j],
                                               phi = phi_vector[j])

        # Calculating the likelihood of the new tree
        likelihood_new <- tree_complete_conditional_het_gpbart_mu(
          tree = new_trees_mu[[j]], # Calculate the full conditional
          residuals = current_partial_residuals,
          x = x, tau_mu = tau_mu,
          precision_vector  = apply(precisions,
                                                      2,
                                                      function(x){exp(sum(log(x)))}),
          nu = nu_vector[j], phi = phi_vector[j],
          number_trees_mu = number_trees_mu
        )

        # Calculating the likelihood of the old tree
        likelihood_old <- tree_complete_conditional_het_gpbart_mu(
          tree = current_trees_mu[[j]], # Calculate the full conditional
          residuals = current_partial_residuals,
          x = x, tau_mu = tau_mu,
          precision_vector = apply(precisions,
                                   2,
                                   function(x){exp(sum(log(x)))}),
          nu = nu_vector[j], phi = phi_vector[j],
          number_trees_mu = number_trees_mu
        )

        # Extracting only the likelihood
        l_new <- likelihood_new$log_posterior +
          tree_prior(
            tree = new_trees_mu[[j]], # Calculate the tree prior
            alpha = alpha_mu,
            beta = beta_mu
          )

        # Extracting only the likelihood
        l_old <- likelihood_old$log_posterior +
          tree_prior(
            tree = current_trees_mu[[j]], # Calculate the tree prior
            alpha = alpha_mu,
            beta = beta_mu
          )

        # Getting the log of transitin prob
        log_transition <- log_transition_prob(current_tree = current_trees_mu[[j]],
                                              new_tree = new_trees_mu[[j]],verb = verb)

        # (log) Probability of accept the new proposed tree
        acceptance <- (l_new - l_old + log_transition)

        # If Storage or not based on thin and burn parameters
        if((i > burn) ) {
          full_cond_store[curr, j] <- l_old
        }

        if(acceptance > 0 || acceptance > -stats::rexp(1)) { #
          acc_ratio <- acc_ratio + 1

          # Checking whether the trees are identical
          if(identical(current_trees_mu[[j]], new_trees_mu[[j]])){
            verb_store_mu[j,"identical"] <- TRUE
          } else {
            verb_store_mu[j,"identical"] <- FALSE
          }

          # Make changes if accept
          current_trees_mu <- new_trees_mu

          # Create a data.frame with the verb and if it was accepted or not
          verb_store_mu[j,"verb"] <- verb
          verb_store_mu[j,"accepted"] <- TRUE

          # Storing likelihood matrix objects
          likelihood_object <- likelihood_new

        } else {
          # Storing likelihood matrix objects
          likelihood_object <- likelihood_old

          # Create a data.frame with the verb and if it was accepted or not
          verb_store_mu[j,"verb"] <- verb
          verb_store_mu[j,"accepted"] <- FALSE
          verb_store_mu[j,"identical"] <- FALSE

        } # End of accept if statement

        # # # To update the mu values
        current_trees_mu[[j]] <- update_mu(
          tree = current_trees_mu[[j]],
          x = x,
          residuals = current_partial_residuals,
          likelihood_object = likelihood_object)

        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar
        predictions[j, ] <- update_residuals(
          tree = current_trees_mu[[j]], x = x,
          residuals = current_partial_residuals,
          phi = phi_vector[j], nu = nu_vector[j],
        )

        # To update phi
        mh_update_phi <- update_phi_marginal(current_tree_iter = current_trees_mu[[j]],
                                               residuals = current_partial_residuals,
                                               x = x,nu = nu_vector[j],phi = phi_vector[j],
                                               gp_variables = gp_variables,
                                               likelihood_object = likelihood_object,
                                               number_trees_mu = number_trees_mu,
                                               discrete_phi = discrete_phi_boolean,
                                               tau_mu = tau_mu,
                                               precisions = precisions,
                                               distance_min = distance_min,
                                               distance_max = distance_max)

          # In case of accept the update over \phi update everything
          if(mh_update_phi$phi_boolean) {

            # Updating the tree and the \phi object from the tree
            current_trees_mu[[j]] <- mh_update_phi$tree

            # Updating the likelihood objects
            likelihood_object <- mh_update_phi$likelihood_object

            # Updating the phi value
            phi_vector[j] <- mh_update_phi$phi_proposal

          } # If doesn't accept, nothing changes.



        # current_partial_residuals_matrix<-
        current_partial_residuals_matrix[j, ] <- current_partial_residuals
        current_predictions_matrix[j, ] <- predictions[j, ]

        current_trees_mu[[j]] <- remove_omega_plus_P_inv(current_tree_iter = current_trees_mu[[j]])
      } # End of Loop through the trees

    }# End of the IF statement of BART boolean

  # =======
  # HERE I GONNA START THE LOOP OVER THE TREES TAU
  # ======

    # Verb Store
    verb_store_tau <- data.frame(verb = rep(NA, number_trees_tau),
                                accepted = rep(NA, number_trees_tau),
                                identical = rep(NA, number_trees_tau))

    for(j in 1:number_trees_tau){

      # Making a exception for the case of only one tree
      if(number_trees_tau == 1) {

        # Getting the current partial values
        current_partial_residuals_precision_squared <- matrix((y_scale^2)/c(precisions), ncol = length(y_scale))
        # MAYBE NEED TO REVIEW THIS LINE

      } else {

        # Getting the current partial values

        current_partial_residual_precision_squared <- ((y_scale - colSums(predictions[,,drop = FALSE]))^2)*apply(precisions[-j,,drop = FALSE],2, function(x){exp(sum(log(x)))})

      }

      # Propose a new tree based on the verbs: grow/prune/change/swap
      verb <- sample(c("grow", "prune", "change", "swap"),
                     prob = c(0.25,0.25,0.4,0.1), size = 1)


      # Case of rotation
      if (i < max(floor(0.1 * burn), 10) || length(current_trees_tau[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations


      # GETTING A NEW TREE
      new_trees_tau <- current_trees_tau # Creating new trees to updated as candidate

      new_trees_tau[[j]] <- update_tree_verb(
        tree = current_trees_tau[[j]],
        x = x,
        node_min_size = node_min_size_tau,
        verb = verb
      )

      # Calculating the likelihood of the new tree
      likelihood_new <- tree_complete_conditional_het_bart_tau(
        tree = new_trees_tau[[j]], # Calculate the full conditional
        precision_sq_residuals_values = current_partial_residual_precision_squared,
        x = x,a_tau = a_tau, d_tau = d_tau
      )

      # Calculating the likelihood of the old tree
      likelihood_old <- tree_complete_conditional_het_bart_tau(
        tree = current_trees_tau[[j]], # Calculate the full conditional
        precision_sq_residuals_values = current_partial_residual_precision_squared,
        x = x,a_tau = a_tau, d_tau = d_tau
      )

      # Extracting only the likelihood
      l_new <- sum(likelihood_new) +
        tree_prior(
          tree = new_trees_tau[[j]], # Calculate the tree prior
          alpha = alpha_tau,
          beta = beta_tau
        )

      # Extracting only the likelihood
      l_old <- sum(likelihood_old) +
        tree_prior(
          tree = current_trees_tau[[j]], # Calculate the tree prior
          alpha = alpha_tau,
          beta = beta_tau
        )

      # Getting the log of transitin prob
      log_transition <- log_transition_prob(current_tree = current_trees_tau[[j]],
                                            new_tree = new_trees_tau[[j]],verb = verb)

      # (log) Probability of accept the new proposed tree
      acceptance <- (l_new - l_old + log_transition)


      # In case of acceptance
      if(acceptance > 0 || acceptance > -stats::rexp(1)) {

        # Counting acc ratio
        acc_ratio <- acc_ratio + 1

        # Make changes if accept
        current_trees_tau <- new_trees_tau

        # Create a data.frame with the verb and if it was accepted or not
        verb_store_tau[j,"verb"] <- verb
        verb_store_tau[j,"accepted"] <- TRUE

      } else {

        # Create a data.frame with the verb and if it was accepted or not
        verb_store_tau[j,"verb"] <- verb
        verb_store_tau[j,"accepted"] <- FALSE

      } # End of accept for MH sample

      # # # To update the mu values
      current_trees_tau[[j]] <- update_tau_het_bart(
        tree = current_trees_tau[[j]],
        x = x,
        precision_sq_residuals_values = current_partial_residual_precision_squared,
        a_tau = a_tau,d_tau = d_tau)

      # EQUATION FROM SECTION 4
      # ==== Using the prediction from R_star_bar
      precisions[j,] <- get_tau_bart(
        tree = current_trees_tau[[j]], x = x
      )

      current_partial_precision_matrix[j,] <- current_partial_residual_precision_squared

    }

  } # End of Loop through the n_inter
  cat("\n")

  # Returning X to its original scale
  if(x_scale) {
    x <- x_original
  }

  results <- list( posterior = list(trees_mu = tree_mu_store,
                                    trees_tau  = tree_tau_store,
                                    phi_store = phi_store,
                                    verb_store_list_mu = verb_store_list_mu,
                                    verb_store_list_tau = verb_store_list_tau,
                                    current_predictions_list = current_predictions_list,
                                    current_precision_list = current_precision_list,
                                    current_partial_residuals_list = current_partial_residuals_list,
                                    y_hat_store = y_hat_store,
                                    y_precision_hat_store = y_precision_hat_store),
                   prior = list(nu_vector = nu_vector,
                                a_tau = a_tau,
                                d_tau = d_tau,
                                alpha_mu = alpha_mu,
                                beta_mu = beta_mu,
                                alpha_tau = alpha_tau,
                                beta_tau = beta_tau,
                                tau_mu = tau_mu,
                                kappa = kappa),
                   mcmc = list(iter = n_iter,
                               burn = burn,
                               bart_boolean = bart_boolean,
                               bart_number_iter = bart_number_iter,
                               acc_ratio = acc_ratio,
                               acc_ratio_phi = acc_ratio_phi,
                               store_size = store_size),
                   control = list(y = y_scale,
                                X = x,
                                x_scale = x_scale,
                                mean_x = mean_x,
                                sd_x = sd_x,
                                scale_boolean = scale_boolean,
                                number_trees_mu = number_trees_mu,
                                number_trees_tau = number_trees_tau,
                                node_min_size_mu = node_min_size_mu,
                                node_min_size_tau = node_min_size_tau,
                                a_min = a_min,
                                b_max = b_max))
      class(results) <- "hgpbart_hGPBART"
    return(results)
}



# #Do a MH for PHI
update_phi_marginal <- function(x, current_tree_iter,
                                residuals,
                                seed = NULL,
                                tau_mu,
                                phi, nu,number_trees_mu,
                                precisions,
                                likelihood_object, p, gp_variables,
                                discrete_phi = TRUE,
                                distance_min,
                                distance_max) {

  # Increased the range of tree proposal
  if(discrete_phi){
    phi_proposal <- sample(c(0.1,0.5,1,5,10), size = 1)
  } else {
    phi_proposal <- stats::runif(1, min = distance_min, max = distance_max)
  }

  # Calculating the likelihood from the new step
  tree_from_phi_proposal <- inverse_omega_plus_P(tree = current_tree_iter,
                                                 x = x,nu = nu,
                                                 precision_vector = apply(precisions,
                                                                          2,
                                                                          function(x){exp(sum(log(x)))}),
                                                 phi = phi_proposal, gp_variables = gp_variables)

  likelihood_phi_proposal <- tree_complete_conditional_het_gpbart_mu(tree = tree_from_phi_proposal,
                                                              x = x,
                                                              residuals = residuals,
                                                              precision_vector = apply(precisions,
                                                                                       2,
                                                                                       function(x){exp(sum(log(x)))}),
                                                              nu = nu, tau_mu = tau_mu,
                                                              phi = phi_proposal,
                                                              number_trees_mu = number_trees_mu)

  # Old phi likelhood
  l_old_phi <- likelihood_object$log_posterior

  # Proposal likelihood
  l_proposal_phi <- likelihood_phi_proposal$log_posterior

  # (log) Probability of accept the new proposed tree
  acceptance_phi <- l_proposal_phi - l_old_phi

  # If storage for phi

  if(acceptance_phi > 0 || acceptance_phi > -stats::rexp(1)) { #

    # Nu boolean to see if was accepted or not
    phi_boolean <- TRUE
    return(list(phi_boolean = phi_boolean,
                likelihood_object = likelihood_phi_proposal,
                tree = tree_from_phi_proposal,
                phi_proposal = phi_proposal)) # Returning the proposal value for phi
  } else {
    # Case of not accepting
    phi_boolean <- FALSE
    return(list(phi_boolean = phi_boolean)) # Returning the old value for phi
  } #
}


# Function to return the depth trees

tree_depth_hist <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)

  for(k in seq_along(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for (i in seq_len(gpbart_model$number_trees)) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- max(vapply(tree, "[[", numeric(1), "depth_node"))
    }
  }
    return(tree_depth)
}



# ORIGINAL PREDICT GAUSSIAN FROM MULTIPLE TREES
predict_gaussian_from_multiple_trees <- function(multiple_trees, # A list of trees
                                                 phi_vector, # A vector of phi values
                                                 nu_vector, # The nu value to be used
                                                 x_train, # The x of the training model
                                                 x_new, # The x that will be predicted
                                                 partial_residuals, # The partial_residual values
                                                 precision_vector_trees,
                                                 pred_bart_only # Boolean argument to predict a BART object
                                                 ) {
  # Defining objects
  y_pred_final <-

  # Calculating the sd from the prediction interval
  y_pred_sd_final <- matrix(0, nrow = length(multiple_trees), ncol = nrow(x_new))

  # print(precision_vector_trees)
  # print(seq)
  # Iterating over all trees
  for(m in seq_along(multiple_trees)) {

    # Creating the list to be predicted (The if is just in case of of just one tree)
    new_tree <- multiple_trees[[m]]
    phi <- phi_vector[m]
    nu <- nu_vector[m]

    # Creating the pred  vector
    y_pred <- numeric(nrow(x_new))
    y_sd_pred <- numeric(nrow(x_new))


    # Setting the root node with the new observations
    new_tree[["node_0"]]$test_index <- seq_len(nrow(x_new))

    # Creating the list of nodes
    list_nodes <- names(new_tree)[-1]

    # IN CASE OF JUST ROOT NODE
    if(length(new_tree) == 1) {
      list_nodes <- "node_0"
    }

    # Updating all nodes
    for(i in seq_along(list_nodes)) {

      current_node_aux <- new_tree[[list_nodes[i]]]

      # In case of more than one node
      if(length(list_nodes) > 1) {
        # Veryfing the type of the current node
        if (is.list(current_node_aux$node_var)) {

          # Rotated Lon and Lat

          rotated_x <- tcrossprod((A(current_node_aux$theta)), x_new[,current_node_aux$node_var$node_var_pair])
          rownames(rotated_x) <- current_node_aux$node_var$node_var_pair

          # Updating observations from the left node
          if(current_node_aux$left == 1) {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index]  < current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index] >= current_node_aux$node_var_split)]
          }
        } else { # Evaluating the case where is not a rotated lat/lon

          # To continous covariates
          if(is.numeric(x_new[, current_node_aux$node_var])) {

            # Updating observations from the left node
            if(current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var]  < current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] >= current_node_aux$node_var_split)]
            }

            # To categorical covariates
          } else {
            # Updating observations from the left node
            if(current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
            }
          }
        }
      }
      # Here I will calculate the predicted values based on the terminal nodes AND only on terminal nodes which have test

      if(new_tree[[list_nodes[[i]]]]$terminal == 1 && length(new_tree[[list_nodes[[i]]]]$test_index) > 0) {

        # Selecting the observations from the current node
        # CHANGE HERE, TO SELECT WHICH ONE WILL BE USED
        x_current_node <- matrix(x_train[new_tree[[list_nodes[[i]]]]$observations_index, ],
                                 nrow = length(new_tree[[list_nodes[[i]]]]$observations_index))

        # Selecting the observations from the test node
        # CHANGE HERE, TO SELECT WHICH ONE WILL BE USED

        x_star <- matrix(x_new[new_tree[[list_nodes[[i]]]]$test_index, ], nrow = length(new_tree[[list_nodes[[i]]]]$test_index))

        # Seleting the precision from this current node from this current tree
        precision_vector <- apply(precision_vector_trees,2,function(x){exp(sum(log(x)))})[new_tree[[list_nodes[[i]]]]$observations_index]

        # Calcualting the distance matrix
        distance_matrix_current_node <- symm_distance_matrix(m1 = x_current_node)

        # Generating the scenario where there are only BART predictions
        if(isFALSE(pred_bart_only)) {

          # Getting the GP from a terminal node
          gp_process <- gp_main(
            x_train = x_current_node, distance_matrix_train = distance_matrix_current_node,
            y_train = matrix((partial_residuals[m,new_tree[[list_nodes[[i]]]]$observations_index]) - new_tree[[list_nodes[[i]]]]$mu,
                              nrow = nrow(x_current_node)),
            x_star = x_star, precision_vector = precision_vector,
            nu = nu, phi = phi
          )

          # Creating the mu vector
          y_pred[new_tree[[list_nodes[i]]]$test_index] <- gp_process$mu_pred + new_tree[[list_nodes[[i]]]]$mu

        } else {

          # Attributing the predicted value as the sampled \mu from the terminal node
          y_pred[new_tree[[list_nodes[i]]]$test_index] <- new_tree[[list_nodes[[i]]]]$mu
        }

      }
    }

    y_pred_final[m, ] <- y_pred
  }
  # Return the new tree
  return(list(y_pred = colSums(y_pred_final),
              all_tree_pred = y_pred_final

              ))
}

# ORIGINAL PREDICT GAUSSIAN FROM MULTIPLE TREES
predict_tau_test <- function(multiple_trees, # A list of trees
                             precision_vector_trees,
                             x_train, # The x of the training model
                             x_new # The x that will be predicted
) {
  # Defining objects
  y_pred_precision_final <- matrix(0, nrow = length(multiple_trees), ncol = nrow(x_new))

  # print(seq_along(multiple_trees))
  # Iterating over all trees
  for(m in seq_along(multiple_trees)) {

    # Creating the list to be predicted (The if is just in case of of just one tree)
    new_tree <- multiple_trees[[m]]

    # Creating the pred  vector
    y_pred_precision <- numeric(nrow(x_new))


    # Setting the root node with the new observations
    new_tree[["node_0"]]$test_index <- seq_len(nrow(x_new))

    # Creating the list of nodes
    list_nodes <- names(new_tree)[-1]

    # IN CASE OF JUST ROOT NODE
    if(length(new_tree) == 1) {
      list_nodes <- "node_0"
    }

    # Updating all nodes
    for(i in seq_along(list_nodes)) {

      current_node_aux <- new_tree[[list_nodes[i]]]

      # In case of more than one node
      if(length(list_nodes) > 1) {
        # Veryfing the type of the current node
        if (is.list(current_node_aux$node_var)) {

          # Rotated Lon and Lat

          rotated_x <- tcrossprod((A(current_node_aux$theta)), x_new[,current_node_aux$node_var$node_var_pair])
          rownames(rotated_x) <- current_node_aux$node_var$node_var_pair

          # Updating observations from the left node
          if(current_node_aux$left == 1) {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index]  < current_node_aux$node_var_split)] # Updating the left node
          } else {
            new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(rotated_x[current_node_aux$node_var$node_var, new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index] >= current_node_aux$node_var_split)]
          }
        } else { # Evaluating the case where is not a rotated lat/lon

          # To continous covariates
          if(is.numeric(x_new[, current_node_aux$node_var])) {

            # Updating observations from the left node
            if(current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var]  < current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] >= current_node_aux$node_var_split)]
            }

            # To categorical covariates
          } else {
            # Updating observations from the left node
            if(current_node_aux$left == 1) {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
            } else {
              new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x_new[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
            }
          }
        }
      }
      # Here I will calculate the predicted values based on the terminal nodes AND only on terminal nodes which have test

      if(new_tree[[list_nodes[[i]]]]$terminal == 1 && length(new_tree[[list_nodes[[i]]]]$test_index) > 0) {

        # Selecting the observations from the current node
        # CHANGE HERE, TO SELECT WHICH ONE WILL BE USED
        x_current_node <- matrix(x_train[new_tree[[list_nodes[[i]]]]$observations_index, ],
                                 nrow = length(new_tree[[list_nodes[[i]]]]$observations_index))

        # Selecting the observations from the test node
        x_star <- matrix(x_new[new_tree[[list_nodes[[i]]]]$test_index, ], nrow = length(new_tree[[list_nodes[[i]]]]$test_index))

        # Seleting the precision from this current node from this current tree
        precision_vector <- precision_vector_trees[m,new_tree[[list_nodes[[i]]]]$observations_index]

        # Attributing the predicted value as the sampled \tau from the terminal node
        y_pred_precision[new_tree[[list_nodes[i]]]$test_index] <- new_tree[[list_nodes[[i]]]]$tau

      }
    }

    y_pred_precision_final[m, ] <- y_pred_precision
  }
  # Return the new tree
  return(list(precision_new = apply(y_pred_precision_final,2,function(x){exp(sum(log(x)))}),
              all_tree_pred_precision = y_pred_precision_final

  ))
}

# Function to count the number of terminal nodes
count_terminal_nodes <- function(tree) {
  return(sum(vapply(tree, "[[", numeric(1), "terminal") == 1))
}

# Predicting a gpbart
#' @method predict "hgpbart_hGPBART"
#' @rdname hgpbart_hGPBART
#' @param x_test the test set
#' @param pred_bart_only boolean if there are only bart predictions
#' @param type select the prediction outputs among 'c("all", "mean","median"))'
#' @param ... other parameters
#' @usage
#' \method{predict}{hgpbart_hGPBART}(object,
#'         x_test,
#'         pred_bart_only = FALSE,
#'         type = c('all', 'median', 'mean'),
#'         ...)
#' @export
predict.hgpbart_hGPBART <- function(object, x_test,
                                  pred_bart_only = FALSE,
                                  type = c("all", "mean", "median"),...) {
  # Matching the type
  type <- match.arg(type)
  # Loading x_train
  if(object$control$x_scale){
    x_train <- as.matrix(scale(object$control$X,center = object$control$mean_x,scale = object$control$sd_x))
    x_test  <- as.matrix(scale(x_test,center = object$control$mean_x,scale = object$control$sd_x))

    # Retrieving the col.names
    colnames(x_train) <- colnames(x_test)  <- colnames(object$control$X)
  } else {
    x_train <- object$X
  }

  # Number of iters of bayesian simulation
  n_iter <- object$mcmc$store_size

  # The number of columns is the number of test observations and the rows are the iterations
  y_hat_matrix <-
  y_hat_precision_matrix <- matrix(0, nrow = n_iter, ncol = nrow(x_test))

  # Getting the training objects
  y_hat_matrix_train <-
  y_hat_precision_matrix_train <- matrix(0, nrow = n_iter, ncol = nrow(x_train))

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running rBART..."
  )
  y_list_matrix <- list()

  # Creating the final vector
  y_pred_final <- matrix(0, nrow = object$control$number_trees_mu, ncol = nrow(x_test))
  cov_pred_final <- list()
  variance <- matrix(0, nrow = nrow(x_test), ncol = nrow(x_test))


  # Having the partial residuals matrix and precisions
  # current_partial_residuals_matrix <- predictions
  # current_partial_precision_matrix <- precisions

  # Looping around the trees
  for(i in seq_len(n_iter)) {

    utils::setTxtProgressBar(progress_bar, i)

    # Selecting one tree from BART model
    current_tree_mu <- object$posterior$trees_mu[[i]]
    current_tree_tau <- object$posterior$trees_tau[[i]]

    # Getting the predictions from the test observations
    y_pred_aux <- predict_gaussian_from_multiple_trees(
      multiple_trees = current_tree_mu,
      x_train = x_train,
      x_new = x_test,
      partial_residuals = object$posterior$current_partial_residuals_list[[i]],
      phi_vector = object$posterior$phi_store[i, ],
      nu_vector = object$prior$nu_vector,
      precision_vector_trees = object$posterior$current_precision_list[[i]],
      pred_bart_only = pred_bart_only
    )

    y_pred_precision_aux <- predict_tau_test(multiple_trees = current_tree_tau,
                                              x_train = x_train,
                                              precision_vector_trees = object$posterior$current_precision_list[[i]],
                                              x_new = x_test)

    # Iterating over all trees (test)
    y_pred_final <- y_pred_aux$all_tree_pred
    y_pred_precision_final <- y_pred_precision_aux$all_tree_pred_precision

    if(object$control$scale_boolean) {
      # Recovering the prediction interval from test
      y_hat_matrix[i, ] <- unnormalize_bart(colSums(y_pred_final), a = object$control$a_min, b = object$control$b_max)
      y_hat_precision_matrix[i,] <-  apply(y_pred_precision_final,2, function(x) { exp(sum(log(x)))})/( (object$control$b_max- object$control$a_min)^2)
    } else {
      y_hat_matrix[i, ] <- colSums(y_pred_aux$all_tree_pred)
      y_hat_precision_matrix[i,] <- apply(y_pred_precision_final,2, function(x) { exp(sum(log(x)))})
    }

    # Storing all
    # y_list_matrix[[i]] <- y_pred_final

  }

  out <- list(
    pred = switch(type,
                  all = y_hat_matrix,
                  mean = colMeans(y_hat_matrix),
                  median = apply(y_hat_matrix, 2, "median")
    ),
    sd = switch(type,
                all = sqrt(1/y_hat_precision_matrix),
                mean = sqrt(1/(colMeans(y_hat_precision_matrix))),
                median = stats::median(1/colMeans(y_hat_precision_matrix), sqrt, numeric(1)))
    )


  return(list(out = out))
}

# A function to count the mean values of observations in terminal nodes
gpbart_count_terminal_nodes <- function(mod_gpbart){

  # Auxiliar matrix
  all_tree_terminal_nodes <- array(NA,dim = c(mod_gpbart$number_trees, 50, mod_gpbart$store_size),
                                   dimnames = list(paste0("tree_", seq_len(mod_gpbart$number_trees)),
                                                   paste0("node_", seq_len(50)),
                                                   paste0("iter_", seq_len(mod_gpbart$store_size))))

  # Iterating over MH
  for(k in seq_len(mod_gpbart$store_size)) {
    tree_iter <- mod_gpbart$trees[[k]]

    # Iterating with for over the terminal nodes
    for(i in seq_len(mod_gpbart$number_trees)) {

      # Gathering the node number
      node_number <- unlist(lapply(tree_iter[[i]], function(x) { x[x$terminal == 1]$node_number}))
      all_tree_terminal_nodes[i,node_number,k] <- unlist(lapply(tree_iter[[i]][names(node_number)],function(x) {
        length(x[x$terminal == 1]$observations_index)
      }))

    }
  }
  # Three splits
  split_nodes <- unique(colMeans(all_tree_terminal_nodes, na.rm = TRUE))
    return(split_nodes[!is.na(split_nodes)])
}



# Retrieve new points on the cubic scale
unit_cube_scale_new_points_uni <- function(x, x_new) {

  # Scaling to -1 to 1 function
  scaling <-
    (2 * x_new - (max(x) + min(x))) / (max(x) - min(x))

  # Applying on all covariates
    return(scaling)
}

# Getting the likelihood from the Gaussian Processes
neg_loglike <- function(prediction_object,
                        y_test) {

  # Getting the prediction
  y_pred <- colMeans(prediction_object$mcmc_pi_mean)

  # Getting the variance
  sd_pred <- colMeans(prediction_object$mcmc_pi_sd)

  return(-stats::dnorm(x = y_test,
                mean = y_pred,
                sd = sqrt(sd_pred), log = TRUE))
}

# sum(neg_loglike(prediction_object = pred_gpbart,y_test = y_test))

# Get tau from terminal nodes
get_tau_values_from_single_tree <- function(gpbart_mod, tree_number, mh_iter = 100){
  # Getting the vector of tau values from the terminal nodes
  return(unlist(lapply(gpbart_mod$trees[[mh_iter]][[tree_number]],function(x) x[x$terminal == 1]$tau)))
}

# Getting Omega Inverse argument list
# get_omega_inverse_terminal_nodes <- function(tree,
#                                              x = x,
#                                              nu = nu_vector[j], phi = phi_vector[j]){
#
#   # Selecting terminal nodes names
#   names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))
#
#   # Selecting the terminal nodes
#   terminal_nodes <- tree[names_terminal_nodes]
#
#   # Number of nodes
#   n_node <- length(terminal_nodes)
#
#   # Picking each node size
#   nodes_size <- sapply(terminal_nodes, function(x) {
#     length(x$observations_index)
#   })
#
#   # Calculating Omega matrix INVERSE
#   Omega_matrix_INV <- mapply(terminal_nodes, FUN = function(y) {
#     chol2inv(  chol(kernel_function(
#       x = matrix(x[y$observations_index, ], nrow = length(y$observations_index)),
#       nu = nu, phi = phi
#     )) )
#   }, SIMPLIFY = FALSE)
#
#   # Adding the mu values calculated
#   for(i in seq_along(names_terminal_nodes)) {
#     tree[[names_terminal_nodes[i]]]$Omega_inv <- Omega_matrix_INV[[names_terminal_nodes[i]]]
#   }
#     return(tree)
# }

# Getting Omega Inverse + Diag Inverse
inverse_omega_plus_P <- function(tree,
                                 x = x,
                                 nu, phi,
                                 precision_vector,
                                 gp_variables = colnames(x)  # Selecting which gp-variables to use

) {
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- sapply(terminal_nodes, function(x) {
    length(x$observations_index)
  })

  # Calculating Omega matrix INVERSE
  distance_matrices <- mapply(terminal_nodes, FUN = function(y) {
    symm_distance_matrix(matrix(x[y$observations_index, gp_variables], nrow = length(y$observations_index)))
  }, SIMPLIFY = FALSE)

  # Calculating Omega
  Omega_matrix <- mapply(distance_matrices, FUN = function(dist_m) {
    kernel_function(
    squared_distance_matrix = dist_m,
    nu = nu, phi = phi)
  }, SIMPLIFY = FALSE)

  # Retrieving Precision from terminal nodes
  precision_terminal_nodes <- lapply(terminal_nodes, function(x) {
    precision_vector[x$observations_index]
  })

  # Checking if diagonal
  is_Omega_diag <- lapply(Omega_matrix, is_diag_matrix)

  # Calculating Omega_plus_I*tau^(-1)
  Omega_matrix_plus_P <- mapply(Omega_matrix,precision_terminal_nodes, FUN = function(omega,prec) {
    omega + diag(1/prec, nrow = nrow(omega)) # Remember that P is p^(-1)
  }, SIMPLIFY = FALSE)

  # Calculating Omega matrix plus I INVERSE
  Omega_matrix_plus_P_INV <- lapply(Omega_matrix_plus_P, FUN = function(omega_plus_p) { # p is the shrinkage factor
    chol2inv(PD_chol(omega_plus_p))
  })

  # Adding the Omega_matrix_plus_I_Inv
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$Omega_plus_P_inv <- Omega_matrix_plus_P_INV[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$distance_matrix <- distance_matrices[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$Omega_matrix <- Omega_matrix[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$is_Omega_diag <- is_Omega_diag[[names_terminal_nodes[i]]]
  }
    return(tree)
}

# # Removing the Omega_plus_I_inv object
remove_omega_plus_P_inv <- function(current_tree_iter) {

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(current_tree_iter, "[[", numeric(1), "terminal") == 1))

  for(i in names_terminal_nodes) {
    current_tree_iter[[i]]$Omega_plus_P_inv <-
    current_tree_iter[[i]]$distance_matrix <-
    current_tree_iter[[i]]$Omega_matrix <-
    current_tree_iter[[i]]$is_Omega_diag <- NULL
  }
    return(current_tree_iter)
}


# Getting the average values from the training predictions
gpbart_train_mean <- function(gpbart_mod) {
  colSums(Reduce("+", gpbart_mod$current_predictions_list)/length(gpbart_mod$current_predictions_list))
}

# Calculating the variance from the training set.
gpbart_training_var <-  function(gpbart_mod) {

  # Number of MCMC samples
  n_mcmc <- length(gpbart_mod$trees)
  # Creating the var vector
  var_train <- matrix(0, nrow = n_mcmc, ncol = length(gpbart_mod$y))

  # Iterating over all trees and getting the terminal nodes
  for(i in seq_len(n_mcmc)) {

    for(m in seq_len(gpbart_mod$number_trees)) {
      # Selecting the current tree
      tree <- gpbart_mod$trees[[i]][[m]]

      # Selecting the terminal nodes
      terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

      # Iterating over the terminal nodes
      for(node in terminal_nodes){

        # In this line I adding up the quantity of 1/tau for each tree
        var_train[i,node$observations_index] <- var_train[i,node$observations_index] + 1/node$tau
      }
    }
  }

  # Returning the var_train vector
    return(var_train)
}


# Calculating the residuals log-likelihood
loglike_residuals <- function(tree,
                              x,
                              current_partial_residuals,
                              phi,
                              nu,p){

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    current_partial_residuals[x$observations_index]
  })

  # Getting the \tau values
  tau_terminal <- lapply(terminal_nodes, "[[", "tau")

  # Getting the \mu values
  mu_terminal <- lapply(terminal_nodes, "[[", "mu")

  Omega_matrix_plus_I <- lapply(terminal_nodes, function(nodes){
    kernel_function(
      squared_distance_matrix = nodes$distance_matrix,
      nu = nu, phi = phi) + diag(p, nrow = nrow(nodes$distance_matrix))
  })

  # Creating the loglikeresiduals_vec
  loglike_residuals_vec <- numeric()
  # Doing a quick for
  for(i in seq_along(residuals_terminal_nodes)){
    loglike_residuals_vec[i] <-  mvtnorm::dmvnorm(x = residuals_terminal_nodes[[i]],
                                                  mean = rep(mu_terminal[[i]], length(residuals_terminal_nodes[[i]])),
                                                  sigma = (1/tau_terminal[[i]]) * Omega_matrix_plus_I[[i]], log = TRUE)
  }
  # Returning the sum of the log_like_residual_vec
    return(sum(loglike_residuals_vec))
}


# NEW UPDATE G
# Update tau_j values
update_g <- function(tree, x, nu, phi, residuals, seed = NULL, p) {
  # set.seed(seed)

  # New g (new vector prediction for g)
  g_new <- rep(NA, length(residuals))

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal")))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Getting the \mu_{j} vector
  mu_values <- vapply(terminal_nodes, "[[", numeric(1), "mu")

  # Selecting the tau values
  tau_j <- vapply(terminal_nodes, "[[", numeric(1), "tau")

  # Calculating Omega matrix
  Omega_matrix_inverse <- mapply(terminal_nodes, FUN = function(nodes) {
    chol2inv(PD_chol(
      kernel_function(squared_distance_matrix = nodes$distance_matrix,
                      nu = nu,
                      phi = phi)
    ))
  }, SIMPLIFY = FALSE)

  # Getting the A matrix inverse
  A_matrix_inv <- mapply(Omega_matrix_inverse, FUN = function(x) {
    chol2inv(PD_chol(diag(p, nrow = nrow(x)) + x))
  }, SIMPLIFY = FALSE)

  # Calculating g_mean posterior
  g_mean <- mapply(A_matrix_inv,
                   residuals_terminal_nodes,
                   mu_values,Omega_matrix_inverse, FUN = function(A_inv, res, mu, omg) {
                     crossprod(A_inv, p * res + mu * rowSums(omg))
                   }, SIMPLIFY = FALSE)

  # Putting in the Keefe's speed order
  g_sd <- mapply(tau_j, Omega_matrix_inverse, FUN = function(tau, omg) {
    tau * (omg + diag(p, nrow=nrow(omg)))
  }, SIMPLIFY = FALSE)

  g_sample <- mapply(FUN=rMVN_var, g_mean, g_sd, SIMPLIFY = FALSE)

  # Adding the mu values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving g
    g_new[terminal_nodes[[i]]$observations_index] <- g_sample[[i]]
  }
    return(g_new)
}

# Some tests over nu parameter
calculate_nu <- function(nu) {
  (nu + 1)/nu
}

calculate_p_nu <- function(nu, p) {
  (p * nu + 1)/nu
}

# Get train predictions
get_train_predictions <- function(gpbart_mod) {

  # Getting the quantile
  gpbart_sum_pred <- do.call(rbind, lapply(gpbart_mod$current_predictions_list, colSums))

  # Returning the matrix of final predictions
    return(gpbart_sum_pred)
}

# Calculating a PI coverage
pi_coverage <- function(y, y_pred, sd_pred, prob = 0.5){

  # CI boundaries
  low_ci <- y_pred + stats::qnorm(p = prob/2)*sd_pred
  up_ci <- y_pred + stats::qnorm(p = 1-prob/2)*sd_pred

  pi_cov <- sum(y<=up_ci & y>=low_ci)/length(y)

  return(pi_cov)
}
