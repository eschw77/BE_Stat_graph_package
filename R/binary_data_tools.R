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


# internal function for computing set of indice permutations and sampling one uniformly at random for the permutation test
#' @param n the number of samples in the data
#' @param seed an optional random seed for reproducibility of the permutation
#' @keywords internal


sample_perm <- function(n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  return(sample(1:n, n, replace = FALSE))
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
#' @param n_splits the number of random splits to perform (default is 2 for a single split, but can be increased for more stable estimates)
#' @param seed an optional random seed for reproducibility of splits
#' @return an upper-triangular matrix of e-values comparing the full model
#'   (with all predictors) to the null model (without each predictor) for each
#'   node. The diagonal and lower-triangular entries are `NA`.
#' @export

wilks_LRT_test_e_val <- function(A, n_splits = 2, seed = NULL) {
  n <- nrow(A)
  p <- ncol(A)
  
  split_e_value <- function(train_A, eval_A, node_j, k) {
    train_j        <- train_A[, node_j]
    train_minus_j  <- train_A[, -node_j, drop = FALSE]
    eval_j         <- eval_A[, node_j]
    eval_minus_j   <- eval_A[, -node_j, drop = FALSE]
    
    local_k        <- k - as.integer(k > node_j)
    train_minus_jk <- train_minus_j[, -local_k, drop = FALSE]
    
    # full model fitted on D1 (eval), null fitted on D0 (train)
    full_beta <- compute_lse(eval_minus_j, eval_j)
    null_beta <- compute_lse(train_minus_jk, train_j)
    
    # both likelihoods evaluated on D0 (train)
    full_log_likelihood <- compute_log_likelihood(full_beta, train_minus_j, train_j)
    null_log_likelihood <- compute_log_likelihood(null_beta, train_minus_jk, train_j)
    
    exp(full_log_likelihood - null_log_likelihood)
  }
  
  # compute cross-fit e-values for a single split
  compute_split_e_values <- function(A0, A1, p) {
    split_e_vals <- matrix(NA_real_, nrow = p, ncol = p)
    
    for (node_j in 1:p) {
      for (k in (1:p)[-node_j]) {
        if (k <= node_j) next
        
        # each direction is a valid e-value per equation (9)
        # average over both directions per equation (10)
        e_01 <- split_e_value(A0, A1, node_j, k)
        e_10 <- split_e_value(A1, A0, node_j, k)
        
        split_e_vals[node_j, k] <- mean(c(e_01, e_10))
      }
    }
    return(split_e_vals)
  }
  
  # accumulate e-value matrices across multiple random splits
  e_values_list <- list()
  n_half <- floor(n / 2)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  for (split in seq_len(n_splits)) {
    idx    <- sample.int(n)
    idx1   <- idx[1:n_half]
    idx2   <- idx[(n_half + 1):n]
    A0     <- A[idx1, , drop = FALSE]  # D0 = train
    A1     <- A[idx2, , drop = FALSE]  # D1 = eval
    
    # compute cross-fit e-values for this split
    e_values_list[[split]] <- compute_split_e_values(A0, A1, p)
  }
  
  # average cross-fit e-value matrices across splits — valid e-value by linearity of expectation
  e_values <- Reduce("+", e_values_list) / n_splits
  return(e_values)
}


#' Run an LRT permutation test for binary data
#' @param A an n x p binary data matrix
#' @param n_permutations the number of permutations to perform for estimating
#'   the null distribution of the test statistic (default is 1000)
#' @param seed an optional random seed for reproducibility of the permutations
#' @return a symmetric p x p matrix of p-values for the likelihood ratio test
#'   comparing the full model (with all predictors) to the null model (without
#'   predictor k) for each response node j. The diagonal and lower triangle are
#'   NA; upper-triangle entries [j, k] give the p-value for the edge (j, k).
#' @export
 
wilks_LRT_permutation_test <- function(A, n_permutations = 1000, seed = NULL) {
  if (!is.matrix(A)) {
    stop("A must be a matrix.")
  }
  if (!is.numeric(n_permutations) || length(n_permutations) != 1 ||
      n_permutations <= 0 || n_permutations %% 1 != 0) {
    stop("n_permutations must be a positive integer.")
  }
 
  if (!is.null(seed)) {
    set.seed(seed)
  }
 
  n <- nrow(A)
  p <- ncol(A)
 
  # ------------------------------------------------------------------
  # 1. Compute the observed test statistics (upper triangle only)
  # ------------------------------------------------------------------
  observed_stats <- matrix(NA_real_, nrow = p, ncol = p)
 
  for (node_j in 1:p) {
    A_j        <- A[, node_j]
    A_minus_j  <- A[, -node_j, drop = FALSE]
    full_beta  <- compute_lse(A_minus_j, A_j)
    full_ll    <- compute_log_likelihood(full_beta, A_minus_j, A_j)
 
    for (k in (1:p)[-node_j]) {
      if (k <= node_j) next                           # upper triangle only
 
      local_k    <- k - as.integer(k > node_j)        # index of col k in A_minus_j
      A_minus_jk <- A_minus_j[, -local_k, drop = FALSE]
      null_beta  <- compute_lse(A_minus_jk, A_j)
      null_ll    <- compute_log_likelihood(null_beta, A_minus_jk, A_j)
 
      observed_stats[node_j, k] <- 2 * (full_ll - null_ll)
    }
  }
 
  # ------------------------------------------------------------------
  # 2. Permutation loop
  #
  #  For each pair (i, k) we keep column k fixed and permute all other
  #  columns (via the rows) jointly.  Crucially, *both* the full and null likelihoods
  #  must be re-estimated from the same permuted data.
  # ------------------------------------------------------------------
  count_exceeds <- matrix(0L, nrow = p, ncol = p)
 
  for (i in 1:p) {
    for (k in (1:p)[-i]) {
      if (k <= i) next                                 # upper triangle only
 
      obs_stat <- observed_stats[i, k]
 
      for (perm in 1:n_permutations) {
        perm_idx <- sample_perm(n)
 
        # Keep column k fixed; permute all other columns jointly.
        A_perm       <- A
        A_perm[, -k] <- A[perm_idx, -k, drop = FALSE]
 
        # --- full model on permuted data ---
        A_j_perm       <- A_perm[, i]
        A_minus_j_perm <- A_perm[, -i, drop = FALSE]
        full_beta_perm <- compute_lse(A_minus_j_perm, A_j_perm)
        full_ll_perm   <- compute_log_likelihood(
                            full_beta_perm, A_minus_j_perm, A_j_perm)
 
        # --- null model on permuted data (drop column k from predictors) ---
        local_k_perm    <- k - as.integer(k > i)      # index of col k in A_minus_j_perm
        A_minus_jk_perm <- A_minus_j_perm[, -local_k_perm, drop = FALSE]
        null_beta_perm  <- compute_lse(A_minus_jk_perm, A_j_perm)
        null_ll_perm    <- compute_log_likelihood(
                            null_beta_perm, A_minus_jk_perm, A_j_perm)
 
        # Both likelihoods come from the same permuted data set.
        perm_stat <- 2 * (full_ll_perm - null_ll_perm)
 
        if (perm_stat >= obs_stat) {
          count_exceeds[i, k] <- count_exceeds[i, k] + 1L
        }
      }
    }
  }
 
  # ------------------------------------------------------------------
  # 3. Assemble p-value matrix
  # ------------------------------------------------------------------
  p_values <- matrix(NA_real_, nrow = p, ncol = p)
 
  upper_idx            <- upper.tri(p_values, diag = FALSE)
  p_values[upper_idx]  <- (count_exceeds[upper_idx] + 1) / (n_permutations + 1)
 
  # Symmetrise
  p_values[lower.tri(p_values, diag = FALSE)] <- t(p_values)[lower.tri(p_values, diag = FALSE)]
  
  return(p_values)
}
