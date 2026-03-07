# R/binary_data_tools.R

#' Functions for following the BEstat framework
#'
#' Power operator function
#' @param X a matrix or vector of binary data, in our case nodes in a graph
#' @return a vector of p probabilities for P(A_i = 1) for i in 1,...,p
#' 
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
J_i <- function(j, p){ 
    return(which(bitwAnd(0:(2^p - 1), 2^(j-1)) != 0))
} 




#' Compute Least Squares Estimates
#'
#' Computes least squares estimates for regression coefficients using power expansion of features.
#' @param X an n x (p-1) matrix of predictors
#' @param Y an n-length vector of the response variable (left out node)
#' @return a vector of coefficient estimates obtained from the least squares solution
#'
compute_lse <- function(X, Y) {
  # X is n x p-1, Y is left out node
  cal_X <- power_expand(X)
  cal_Y <- Y
  beta_hat <- ginv(cal_X) %*% cal_Y
  return(beta_hat)
}



#' Estimate beta coefficients for a BELIEF regression
#' The BELIEF model is E[B|A] = beta^T A^\otimes, where A^\otimes is the power expansion of A. This function estimates the beta coefficients using least squares estimation.
#' @param A an n x p matrix of binary data for the nodes in the graph
#' return