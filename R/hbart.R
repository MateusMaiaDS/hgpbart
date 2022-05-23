# Function to calculate the tree complete conditional using BART model
tree_complete_conditional_het_bart_mu <- function(x,
                                                  residuals_values,
                                                  tree,
                                                  tau_mu,
                                                  precision_vector) {

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
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals_values[x$observations_index]
  })

  # Retrieving Residuals terminal nodes
  precision_terminal_nodes <- lapply(terminal_nodes, function(x) {
    precision_vector[x$observations_index]
  })

  # Retrieving Residuals terminal nodes
  precision_sum_terminal_nodes <- unlist(lapply(terminal_nodes, function(x) {
    sum(precision_vector[x$observations_index])
  }))


  # Calculating the  log(p(x_{i}))
  sum_log_precision_terminal_nodes <- unlist(lapply(terminal_nodes, function(x){
    sum(log(precision_vector[x$observations_index]))
  }))

  # \sum_p(x_i)*r^2
  sum_precision_residuals_squared <- unlist(mapply(precision_terminal_nodes,
                                                   residuals_terminal_nodes,
                                                   FUN = function(precision,residuals) {
                                                     sum(precision*residuals^2)
                                                   },SIMPLIFY = FALSE))

  # Calculating the "last term" \frac{(\sum_{p*ri})^2}{\sum_p + \tau_\mu}
  squared_sum_precision_residuals <- unlist( mapply(precision_sum_terminal_nodes,
                                                    precision_terminal_nodes,
                                                    residuals_terminal_nodes, FUN = function(precision_sum,
                                                                                             precision,
                                                                                             residuals){
                                                      (sum(precision*residuals)^2)/(precision_sum+tau_mu)
                                                    }, SIMPLIFY = FALSE))

  # Retrieve all nodes values and calculate all of them
  log_posterior <- 0.5*sum_log_precision_terminal_nodes+0.5*log(precision_sum_terminal_nodes+tau_mu)-
    0.5*sum_precision_residuals_squared-0.5*squared_sum_precision_residuals



  return(log_posterior)
}

# Function to calculate the tree complete conditional using het-BART model (FOR \tau)
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


# Update \mu using the BART simple version
update_mu_het_bart <- function(tree,
                               x,
                               tau,
                               tau_mu,
                               residuals,
                               precision_vector,
                               seed = NULL) {

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
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })

  # Retrieving Residuals terminal nodes
  precision_terminal_nodes <- lapply(terminal_nodes, function(x) {
    precision_vector[x$observations_index]
  })

  # Retrieving Residuals terminal nodes
  precision_sum_terminal_nodes <- unlist(lapply(terminal_nodes, function(x) {
    sum(precision_vector[x$observations_index])
  }))


  # Calculating the mean of \mu_post
  mu_post_mean <- mapply(residuals_terminal_nodes,
                         precision_terminal_nodes,
                         precision_sum_terminal_nodes, FUN = function(resid,prec,prec_sum){
                           sum(resid*prec)/(prec_sum+tau_mu)
                         },SIMPLIFY = FALSE)

  mu_post_var <- (precision_sum_terminal_nodes+tau_mu)^(-1)

  # Calculating mu values
  mu <- mapply(mu_post_mean, mu_post_var,
               FUN = function(x, y) {
                 stats::rnorm(
                   n = 1,
                   mean = x,
                   sd = sqrt(mu_post_var)
                 )
               }, SIMPLIFY = FALSE)

  # Adding the mu values calculated
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$mu <- mu[[names_terminal_nodes[i]]]
  }
  return(tree)
}

# Getting the \mu vector from terminal nodes
get_mu_bart <- function(tree, x) {

  # New g (new vector prediction for g)
  predictions_new <- rep(NA, nrow(x))

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Getting the \mu_{j} vector
  mu_values <- vapply(terminal_nodes, "[[", numeric(1), "mu")

  # Adding the mu values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving g
    predictions_new[terminal_nodes[[i]]$observations_index] <- mu_values[[i]]
  }

  return(predictions_new)
}

# Getting the \tau vector from terminal nodes
get_tau_bart <- function(tree, x) {

  # New g (new vector prediction for g)
  predictions_new <- rep(NA, nrow(x))

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Getting the \tau_{j} vector
  tau_values <- vapply(terminal_nodes, "[[", numeric(1), "tau")

  # Adding the tau values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving tau
    predictions_new[terminal_nodes[[i]]$observations_index] <- tau_values[[i]]
  }

  return(predictions_new)
}
