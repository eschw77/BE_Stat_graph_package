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
#' @return a matrix of p-values for the likelihood ratio test comparing the full model (with all predictors) to the null model (without each predictor) for each node. The diagonal entries are NA since we do not test nodes against themselves.
#' The matrix is symmetric with NA on the diagonal, with p-values for the tests in the off-diagonal entries to allow for multiple testing control on repeated tests 
#' @export

wilks_LRT_test <- function(A) {
  p_values <- matrix(NA_real_, nrow = ncol(A), ncol = ncol(A))
  p <- ncol(A)

  for (node_j in 1:p) {
    A_j <- A[, node_j]
    A_minus_j <- A[, -node_j, drop = FALSE]

    full_beta <- compute_lse(A_minus_j, A_j)
    full_log_likelihood <- compute_log_likelihood(full_beta, A_minus_j, A_j)

    for (k in (1:p)[-node_j]) {
      if (k <= node_j) {
        next
      }

      # Convert global column index k to local index in A_minus_j, accounting for the removed node_j
      local_k <- k - as.integer(k > node_j)
      A_minus_jk <- A_minus_j[, -local_k, drop = FALSE]

      null_beta <- compute_lse(A_minus_jk, A_j)
      null_log_likelihood <- compute_log_likelihood(null_beta, A_minus_jk, A_j)

      test_statistic <- 2 * (full_log_likelihood - null_log_likelihood)
      p_value <- pchisq(test_statistic, df = (2^p - 2^(p - 1)), lower.tail = FALSE)
      p_values[node_j, k] <- p_value
    }
  }

  return(p_values)
}

#' Construct split likelihood-ratio e-values for binary data
#' The full model is fit on one split and evaluated on the other split, while
#' the null model is re-fit on the evaluation split so the denominator is the
#' null likelihood maximized on the held-out data.
#' @param A an n x p binary data matrix
#' @return an upper-triangular matrix of e-values comparing the full model
#'   (with all predictors) to the null model (without each predictor) for each
#'   node. The diagonal and lower-triangular entries are `NA`.
#' @export

wilks_LRT_test_e_val <- function(A) {
  e_values <- matrix(NA_real_, nrow = ncol(A), ncol = ncol(A))
  n <- nrow(A)
  p <- ncol(A)

  # Split the data in half into two matrices, randomly
  n_half <- floor(n / 2)
  idx <- sample.int(n)
  idx1 <- idx[1:n_half]
  idx2 <- idx[(n_half + 1):n]
  A0 <- A[idx1, , drop = FALSE]
  A1 <- A[idx2, , drop = FALSE]

  split_e_value <- function(train_A, eval_A, node_j, k) {
    train_j <- train_A[, node_j]
    train_minus_j <- train_A[, -node_j, drop = FALSE]
    eval_j <- eval_A[, node_j]
    eval_minus_j <- eval_A[, -node_j, drop = FALSE]

    local_k <- k - as.integer(k > node_j)
    eval_minus_jk <- eval_minus_j[, -local_k, drop = FALSE]
    
    # estimators fitted on opposite halves of the data
    full_beta <- compute_lse(eval_minus_j, eval_j)           # D1, alternative
    null_beta <- compute_lse(train_minus_jk, train_j)        # D0, null

    #log-likelihoods evaluated on A0 for both models
    full_log_likelihood <- compute_log_likelihood(full_beta, train_minus_j, train_j)
    null_log_likelihood <- compute_log_likelihood(null_beta, train_minus_jk, train_j)

    exp(full_log_likelihood - null_log_likelihood)
  }

  for (node_j in 1:p) {
    for (k in (1:p)[-node_j]) {
      if (k <= node_j) {
        next
      }
      
      e_01 <- split_e_value(A0, A1, node_j, k)
      e_10 <- split_e_value(A1, A0, node_j, k)

      e_values[node_j, k] <- mean(c(e_01, e_10))
    }
  }
  return(e_values)
}