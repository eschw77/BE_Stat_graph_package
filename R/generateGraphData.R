# R/generateGraphData.R    
#' Generate Binary Graph Data with CI Constraints (Ising Model, ±1 encoding)
#'
#' Enforces CI conditions X_i ⊥ X_j | X_rest iff (i,j) is not an edge.
#' Uses an Ising MRF: P(X) ∝ exp(X^T Θ X / 2 + b^T X), X ∈ {-1, +1}^n.
#' Bias terms b_i are calibrated so that P(X_i = +1) ≈ target_probs[i].
#' Similar to the Cai et al. 2019 paper on MRFs and the gut microbiome
#' @name isinggraph-package
#' @docType package
NULL

# ── Internal helpers ──────────────────────────────────────────────────────────
#' Build the coupling matrix from an edge list
#' @keywords internal
.build_theta <- function(n_nodes, edges, weak_edges = NULL, weights = NULL) {
  theta <- matrix(0, n_nodes, n_nodes)
  # If not in edges, encode conditional independence (zero coupling)
  if (is.null(edges) || nrow(edges) == 0) return(theta)
  # If weights not provided, assign random weights in [-0.8, -0.3] ∪ [0.3, 0.8]
  # stay away from poles of 0,1 to ideally have a model with less extreme conditional relations
  if (is.null(weights)) {
    weights <- runif(nrow(edges), 0.3, 0.8) * sample(c(-1, 1), nrow(edges), replace = TRUE)
  }
  # Fill in the symmetric coupling matrix Theta
  for (k in seq_len(nrow(edges))) {
    i <- edges[k, 1]; j <- edges[k, 2]
    theta[i, j] <- weights[k]
    theta[j, i] <- weights[k]
  }
  # Handle weak edges if provided
  if (!is.null(weak_edges)) {
    weak_weights <- runif(nrow(weak_edges), 0.05, 0.15 ) * sample(c(-1, 1), nrow(weak_edges), replace = TRUE)
    for (k in seq_len(nrow(weak_edges))) {
      i <- weak_edges[k, 1]; j <- weak_edges[k, 2]
      theta[i, j] <- weak_weights[k]
      theta[j, i] <- weak_weights[k]
    }
  }
  theta
}


#' Single Gibbs sweep over all nodes (±1 encoding)
#' @keywords internal
#' For each node i, sample X_i given the current state of all other nodes and 
#' the parameters theta and bias 
.gibbs_sweep <- function(x, theta, bias) {
  n <- length(x)
  for (i in seq_len(n)) {
    # Conditional: log-odds of x_i = +1 vs -1
    # Delta energy when flipping x_i from s to +1: 2*(sum_j theta[i,j]*x[j] + bias[i])
    neighbor_sum <- sum(theta[i, ] * x) - theta[i, i] * x[i]
    log_odds     <- 2.0 * (neighbor_sum + bias[i])
    p_plus       <- 1.0 / (1.0 + exp(-log_odds))
    x[i]         <- ifelse(runif(1) < p_plus, 1.0, -1.0)
  }
  x
}


#' Run Gibbs sampler
#' @keywords internal
#' Draw samples from the Ising model defined by theta and bias
.gibbs_sample <- function(theta, bias, n_samples, burn_in) {
  n  <- length(bias)
  x  <- sample(c(-1, 1), n, replace = TRUE)
  mat <- matrix(0L, nrow = n_samples, ncol = n)

  for (s in seq_len(burn_in)) x <- .gibbs_sweep(x, theta, bias)
  for (s in seq_len(n_samples)) {
    x       <- .gibbs_sweep(x, theta, bias)
    mat[s, ] <- as.integer(x)
  }
  mat
}


#' Estimate marginal P(X_i = +1) via short Gibbs run
#' @keywords internal
.estimate_marginals <- function(theta, bias, n_est = 7000, burn_in = 500) {
  samp <- .gibbs_sample(theta, bias, n_est, burn_in)
  colMeans(samp == 1L)
}


#' Calibrate bias vector so P(X_i=+1) ≈ target_probs
#'
#' Uses coordinate-wise Newton updates on the marginal mismatch with adaptive step size.
#' @keywords internal
.calibrate_bias <- function(theta, target_probs, n_est = 10000, burn_in = 1000,
                             tol = 0.01, max_iter = 30, verbose = FALSE) {
  n    <- length(target_probs)
  bias <- qlogis(target_probs)   # initialise with marginal log-odds (ignores coupling)

  for (iter in seq_len(max_iter)) {
    est <- .estimate_marginals(theta, bias, n_est, burn_in)
    err <- est - target_probs
    if (verbose) {
      cat(sprintf("Iter %2d | max|err| = %.4f\n", iter, max(abs(err))))
    }
    if (max(abs(err)) < tol) break
    # Adaptive Newton step with damping and gradient clipping
    # d(logit(p))/dp = 1/(p*(1-p))
    grad <- err / pmax(est * (1 - est), 1e-3)
    # Clip gradient to prevent extreme updates
    grad <- pmax(pmin(grad, 2), -2) 
    # Use adaptive step size: larger when error is large, smaller when close
    max_err <- max(abs(err))
    step_size <- ifelse(max_err > 0.1, 0.3, 0.7)
    bias <- bias - step_size * grad
  }
  list(bias = bias, final_est = est, converged = max(abs(err)) < tol)
}


# ── Main exported functions ───────────────────────────────────────────────────

#' Generate binary (±1) samples from an Ising MRF with specified graph structure
#' and target marginal probabilities.
#'
#' @param n_nodes   Integer. Number of nodes.
#' @param edges     Two-column integer matrix of edges (1-indexed). Each row is
#'                  one undirected edge (i, j) with i < j. Can be NULL or empty
#'                  if \code{weak_edges} is provided; otherwise required.
#' @param weak_edges Two-column integer matrix of edges (1-indexed). Each row is
#'                  one undirected edge (i, j) with i < j, expected to have weaker coupling than those in edges, 
#'                  default NULL. If provided, these edges will be assigned smaller random weights than those in edges.
#' @param target_probs Numeric vector of length \code{n_nodes} with target
#'                  P(X_i = +1) for each node. Values must be in (0, 1).
#'                  Default: all 0.5.
#' @param weights   Optional numeric vector of edge weights (same row order as
#'                  \code{edges}). If NULL, random weights in [-0.8, -0.3] ∪
#'                  [0.3, 0.8] are used.
#' @param n_samples Integer. Number of samples to return.
#' @param burn_in   Integer. Gibbs burn-in steps.
#' @param calibrate Logical. If TRUE (default), bias terms are tuned to match
#'                  \code{target_probs}.
#' @param tol       Convergence tolerance for marginal calibration.
#' @param max_iter  Maximum calibration iterations.
#' @param verbose   If TRUE, print calibration progress.
#' @param seed      Optional random seed for reproducibility.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{samples}}{Integer matrix (n_samples × n_nodes) of ±1 values.}
#'     \item{\code{theta}}{Coupling matrix used.}
#'     \item{\code{bias}}{Calibrated bias vector.}
#'     \item{\code{target_probs}}{Target marginal probabilities.}
#'     \item{\code{estimated_probs}}{Estimated marginals from calibration run.}
#'     \item{\code{edges}}{Edge matrix.}
#'     \item{\code{colnames}}{Column names for the sample matrix.}
#'   }
#'
#' @examples
#' edges <- matrix(c(1,2, 2,3, 3,4, 4,5, 1,5), ncol = 2, byrow = TRUE) for CI 
#' probs <- c(0.7, 0.4, 0.6, 0.3, 0.8)
#' result <- ising_generate(
#'   n_nodes      = 5,
#'   edges        = edges,
#'   target_probs = probs,
#'   n_samples    = 1000,
#'   seed         = 42
#' )
#' head(result$samples)
#' result$estimated_probs
#'
#' @export
ising_generate <- function(n_nodes,
                           edges,
                           weak_edges = NULL,
                           target_probs = rep(0.5, n_nodes),
                           weights      = NULL,
                           n_samples    = 1000,
                           burn_in      = 500,
                           calibrate    = TRUE,
                           tol          = 0.01,
                           max_iter     = 30,
                           verbose      = FALSE,
                           seed         = NULL) {

  # ── Input validation ──────────────────────────────────────────────────────
  stopifnot(
    is.numeric(n_nodes), n_nodes >= 2,
    length(target_probs) == n_nodes,
    all(target_probs > 0 & target_probs < 1),
    n_samples >= 1,
    burn_in   >= 0
  )

  # Allow edges to be NULL or empty if weak_edges is non-empty
  if (is.null(edges) || (is.matrix(edges) && nrow(edges) == 0)) {
    if (is.null(weak_edges) || (is.matrix(weak_edges) && nrow(weak_edges) == 0)) {
      stop("Either edges or weak_edges must be non-empty.")
    }
    edges <- matrix(integer(0), ncol = 2)
  } else {
    stopifnot(is.matrix(edges), ncol(edges) == 2)
    edges <- matrix(as.integer(edges), ncol = 2)
    if (any(edges < 1 | edges > n_nodes))
      stop("Edge indices must be between 1 and n_nodes.")
    if (any(edges[, 1] == edges[, 2]))
      stop("Self-loops are not allowed.")
  }

  if (!is.null(weak_edges)) {
    weak_edges <- matrix(as.integer(weak_edges), ncol = 2)
    if (any(weak_edges < 1 | weak_edges > n_nodes))
      stop("Weak edge indices must be between 1 and n_nodes.")
    if (any(weak_edges[, 1] == weak_edges[, 2]))
      stop("Self-loops are not allowed in weak edges.")
    # Ensure weak edges do not overlap with strong edges
    strong_set <- apply(edges, 1, function(r) paste(sort(r), collapse = "-"))
    weak_set   <- apply(weak_edges, 1, function(r) paste(sort(r), collapse = "-"))
    if (length(intersect(strong_set, weak_set)) > 0)
      stop("Weak edges cannot overlap with strong edges.")
  }



  # ── Set seed for reproducibility ──────────────────────────────────────────
  if (!is.null(seed)) set.seed(seed)

  # ── Build coupling matrix ─────────────────────────────────────────────────
  theta <- .build_theta(n_nodes, edges, weak_edges, weights)

  # ── Calibrate bias terms, i.e. adjust for target probabilities ────────────
  if (calibrate) {
    if (verbose) cat("Calibrating bias terms...\n")
    cal   <- .calibrate_bias(theta, target_probs,
                             tol = tol, max_iter = max_iter,
                             verbose = verbose)
    bias  <- cal$bias
    est_p <- cal$final_est
    if (!cal$converged && verbose)
      warning("Bias calibration did not fully converge. Try increasing max_iter or n_samples.")
  } else {
    bias  <- qlogis(target_probs)
    est_p <- .estimate_marginals(theta, bias)
  }

  # ── Draw samples ──────────────────────────────────────────────────────────
  samp <- .gibbs_sample(theta, bias, n_samples, burn_in)
  colnames(samp) <- paste0("A", seq_len(n_nodes))

  list(
    samples        = samp,
    theta          = theta,
    bias           = bias,
    target_probs   = target_probs,
    estimated_probs = est_p,
    edges          = edges,
    weak_edges     = weak_edges
  )
}


#' Empirically test all pairwise CI conditions in generated data
#'
#' For each pair (i, j), computes a conditional log-odds ratio (stratified
#' over all configurations of the remaining nodes). A value near zero
#' indicates conditional independence.
#'
#' @param samples Integer matrix of ±1 values (n_samples × n_nodes).
#' @param edges   Two-column integer matrix of graph edges (1-indexed).
#' @param min_stratum_size Minimum stratum count to include in averaging.
#'
#' @return A data frame with columns: node_i, node_j, is_edge,
#'         mean_log_OR, sd_log_OR, n_strata.
#'
#' @examples
#' edges  <- matrix(c(1,2, 2,3, 3,4), ncol = 2, byrow = TRUE)
#' result <- ising_generate(4, edges, seed = 1)
#' ising_ci_test(result$samples, edges)
#' 
#'              X_j = +1    X_j = -1
#' X_i = +1      n11         n1m
#' X_i = -1      nm1         nmm
#' OR = (n11 * nmm) / (n1m * nm1)  => log(OR) = log((n11 * nmm) / (n1m * nm1))
#' this smooths for cells that are zero by adding 0.5 to prevent log(0) issues
#' @export
ising_ci_test <- function(samples, edges, min_stratum_size = 10) {
  n      <- ncol(samples)
  pairs  <- combn(n, 2, simplify = FALSE)
  edge_set <- apply(edges, 1, function(r) paste(sort(r), collapse = "-"))

  results <- lapply(pairs, function(p) {
    i <- p[1]; j <- p[2]
    others <- setdiff(seq_len(n), c(i, j))
    is_edge <- paste(sort(c(i, j)), collapse = "-") %in% edge_set

    strata      <- unique(samples[, others, drop = FALSE])
    log_ors     <- numeric(0)

    for (r in seq_len(nrow(strata))) {
      cond_val <- strata[r, ]
      mask <- apply(samples[, others, drop = FALSE], 1,
                    function(row) all(row == cond_val))
      sub  <- samples[mask, , drop = FALSE]
      if (nrow(sub) < min_stratum_size) next

      n11 <- sum(sub[, i] ==  1 & sub[, j] ==  1)
      n1m <- sum(sub[, i] ==  1 & sub[, j] == -1)
      nm1 <- sum(sub[, i] == -1 & sub[, j] ==  1)
      nmm <- sum(sub[, i] == -1 & sub[, j] == -1)
      lor <- log((n11 + 0.5) * (nmm + 0.5) / ((n1m + 0.5) * (nm1 + 0.5)))
      log_ors <- c(log_ors, lor)
    }

    data.frame(
      node_i       = paste0("A", i),
      node_j       = paste0("A", j),
      is_edge      = is_edge,
      mean_log_OR  = if (length(log_ors)) mean(log_ors) else NA_real_,
      sd_log_OR    = if (length(log_ors)) sd(log_ors)   else NA_real_,
      n_strata     = length(log_ors)
    )
  })

  do.call(rbind, results)
}


#' Print a summary of an ising_generate result
#'
#' @param result List returned by \code{ising_generate}.
#' @export
ising_summary <- function(result) {
  n <- ncol(result$samples)
  cat("=== Ising Graph Data Summary ===\n")
  cat(sprintf("Nodes: %d   |   Edges: %d   |   Samples: %d\n",
              n, nrow(result$edges), nrow(result$samples)))
  cat("\nMarginal P(A_i = +1):\n")
  df <- data.frame(
    Node        = paste0("A", seq_len(n)),
    Target      = round(result$target_probs, 3),
    Estimated   = round(result$estimated_probs, 3),
    Empirical   = round(colMeans(result$samples == 1L), 3)
  )
  print(df, row.names = FALSE)
  cat("\nEdge list:\n")
  apply(result$edges, 1, function(r)
    cat(sprintf("  A%d -- A%d  (θ = %.3f)\n", r[1], r[2],
                result$theta[r[1], r[2]]))
  )
  invisible(result)
}

#' Generate Target Marginal Probabilities
#'
#' Generates a vector of target marginal probabilities P(X_i = +1) for use with
#' \code{\link{ising_generate}}. Supports several distribution types.
#'
#' @param p      Integer. Number of nodes (length of probability vector).
#' @param type   Character. Distribution type. Options:
#'   \describe{
#'     \item{\code{"unif"}}{Uniform random from [0, 1].}
#'     \item{\code{"equidistant"}}{Evenly spaced from 0.1 to 0.9.}
#'     \item{\code{"bimodal"}}{First half uniform in [0.1, 0.3], second half in [0.7, 0.9].}
#'     \item{other}{Default: all values set to 0.5.}
#'   }
#'
#' @return Numeric vector of length \code{p} with values in (0, 1).
#'
#' @examples
#' probs_unif <- generate_probs(5, "unif")
#' probs_bimodal <- generate_probs(10, "bimodal")
#' probs_default <- generate_probs(3, "other")
#'
#' @export
generate_probs <- function(p, type){
    if(type == "unif"){
        probs <- runif(p)
    } else if(type == "equidistant"){
        probs <- seq(0.1, 0.9, length.out = p)
    }
    else if (type == "bimodal"){
        probs <- c(runif(p/2, 0.1, 0.3), runif(p/2, 0.7, 0.9))
    } else {
       probs <- rep(0.5, p) # default to 0.5 if no type is specified
    }
    return(probs)
}