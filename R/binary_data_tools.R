# R/binary_data_tools.R

#' Functions for following the BEstat framework
#'
#' Power operator function
#' @param X a matrix or vector of binary data, in our case nodes in a graph
#' @return a vector of p probabilities for P(A_i = 1) for i in 1,...,p
#' @export
power_expand <- function(X) { 
    n <- nrow(X)
    p <- ncol(X)

    X_pow <- matrix(1, nrow = n, ncol = 1)

  for(i in 1:p) {
    X_pow <- cbind(X_pow, X_pow * X[, i])
  }
  return(X_pow)
}

#' Select the indices for a RV A_i of the vector A to capture the indices of all interactions that include A_i
#' We only use the dimension of data 
#' @param p the number of nodes in the graph or RVs in a BELIEF regression
#' @param i the index of the node or RV for which we want to capture interactions
#' @return a vector of indices for the interactions that include A_i
#' @export
J_i <- function(j, p){ 
    return(which(bitwAnd(0:(2^p - 1), 2^(j-1)) != 0))
} 




#' Compute Least Squares Estimates
#'
#' Computes least squares estimates for regression coefficients using power expansion of features.
#' @param X an n x (p-1) matrix of predictors
#' @param Y an n-length vector of the response variable (left out node)
#' @return a vector of coefficient estimates obtained from the least squares solution
#' @export 
compute_lse <- function(X, Y) {
  # X is n x p-1, Y is left out node
  cal_X <- power_expand(X)
  cal_Y <- Y
  beta_hat <- ginv(cal_X) %*% cal_Y
  return(beta_hat)
}

# internal function for computing log-likelihood given the MLE 
#' @param beta_hat a vector of coefficient estimates obtained from the least squares solution
#' @param X an n x (p-1) matrix of predictors
#' @param Y an n-length vector of the response variable (left out node)
#' @keywords internal

compute_log_likelihood <- function(beta_hat, X, Y) {
  cal_X <- power_expand(X)
  n <- nrow(cal_X)
  eta <- as.vector(cal_X %*% as.numeric(beta_hat))
  terms <- 1 + as.numeric(Y) * eta

  if (any(terms <= 0)) {
    warning("Some terms in the log-likelihood are non-positive, returning -Inf")
    return(-Inf)
  }

  -n * log(2) + sum(log(terms))
}

#' Run the Wilks Likelihood Ratio Test for Binary Data
#' @param A an n x p binary data matrix
#' @param node_j the index of node we are regressing on
#' @return a matrix of p-values for the likelihood ratio test comparing the full model (with all predictors) to the null model (without each predictor) for node_j. The diagonal entry corresponding to node_j is NA since we do not test the node against itself.
#' The matrix should be symmetric and have NA on the diagonal, with p-values for the tests in the off-diagonal entries. This hopefully should prevent multiple testing control on repeated tests 
#' 

wilks_LRT_test <- function(A) {
  p_values <- matrix(NA, nrow = ncol(A), ncol = ncol(A))
  p_values[node_j, node_j] <- NA  # No test for the node against itself
  n <- nrow(A)
  p <- ncol(A)
  for (node_j in 1:p) {
    # Define response and predictors
    A_j <- A[, node_j]
    A_minus_j <- A[, -node_j, drop = FALSE]
    
    # Fit full model
    full_beta <- compute_lse(A_minus_j, A_j)
    full_log_likelihood <- compute_log_likelihood(full_beta, A_minus_j, A_j)
    
    # Test over all nodes neq j
    # To get a upper triangular matrix of p-values, we only test for k > j
    for (k in (1:p)[-node_j]) {
      if (k <= node_j) {
        next  # Skip tests for k <= j to maintain upper triangular structure
      }

      # Fit null model without node k
      A_minus_jk <- A_minus_j[, -which((1:(p-1)) == k), drop = FALSE]
      null_beta <- compute_lse(A_minus_jk, A_j)
      null_log_likelihood <- compute_log_likelihood(null_beta, A_minus_jk, A_j)
    
      # Compute test statistic and p-value
      test_statistic <- 2 * (full_log_likelihood - null_log_likelihood)
      p_value <- pchisq(test_statistic, df = 2^p - 2^(p-1), lower.tail = FALSE)
      p_values[node_j, k] <- p_value 
    }
  }
  return(p_values)
}