# tests/test-isinggraph.R
test_that("ising_generate returns correct dimensions", {
  edges  <- matrix(c(1,2, 2,3, 3,4), ncol = 2, byrow = TRUE)
  result <- ising_generate(4, edges, n_samples = 200, burn_in = 100, seed = 1)

  expect_equal(nrow(result$samples), 200)
  expect_equal(ncol(result$samples), 4)
  expect_true(all(result$samples %in% c(-1L, 1L)))
})

test_that("samples are exactly ±1", {
  edges  <- matrix(c(1,2, 2,3), ncol = 2, byrow = TRUE)
  result <- ising_generate(3, edges, n_samples = 100, seed = 42)
  expect_setequal(unique(as.vector(result$samples)), c(-1L, 1L))
})

test_that("marginal calibration reaches target within tolerance", {
  edges  <- matrix(c(1,2, 2,3, 1,3), ncol = 2, byrow = TRUE)
  probs  <- c(0.7, 0.3, 0.6)
  result <- ising_generate(3, edges,
                            target_probs = probs,
                            n_samples    = 3000,
                            burn_in      = 500,
                            tol          = 0.02,
                            seed         = 7)
  emp <- colMeans(result$samples == 1L)
  expect_true(all(abs(emp - probs) < 0.08),
              info = paste("Empirical:", round(emp, 3), "Target:", probs))
})

test_that("ising_ci_test flags edge vs non-edge pairs correctly", {
  edges  <- matrix(c(1,2, 2,3, 3,4, 4,5, 1,5), ncol = 2, byrow = TRUE)
  result <- ising_generate(5, edges, n_samples = 3000, burn_in = 500, seed = 99)
  ci     <- ising_ci_test(result$samples, edges)

  edge_lor    <- mean(abs(ci$mean_log_OR[ ci$is_edge]), na.rm = TRUE)
  nonedge_lor <- mean(abs(ci$mean_log_OR[!ci$is_edge]), na.rm = TRUE)
  expect_gt(edge_lor, nonedge_lor)
})

test_that("invalid inputs throw errors", {
  expect_error(ising_generate(1, matrix(c(1,2), ncol=2)))          # n_nodes < 2
  expect_error(ising_generate(3, matrix(c(1,4), ncol=2)))          # out-of-range edge
  expect_error(ising_generate(3, matrix(c(1,2), ncol=2),
                               target_probs = c(0, 0.5, 0.5)))     # prob = 0
})

test_that("wilks_LRT_test_e_val maximizes the null on held-out data", {
  A <- matrix(
    c(
      0, 0, 0, 0,
      0, 0, 1, 1,
      0, 1, 0, 1,
      0, 1, 1, 0,
      1, 0, 0, 1,
      1, 0, 1, 0,
      1, 1, 0, 0,
      1, 1, 1, 1,
      0, 0, 0, 1,
      0, 1, 0, 0,
      1, 0, 0, 0,
      1, 1, 1, 0
    ),
    ncol = 4,
    byrow = TRUE
  )

  compute_manual_log_likelihood <- function(beta_hat, X, Y) {
    cal_X <- power_expand(X)
    n <- nrow(cal_X)
    eta <- as.vector(cal_X %*% as.numeric(beta_hat))
    terms <- 1 + as.numeric(Y) * eta

    if (any(terms <= 0)) {
      return(-Inf)
    }

    -n * log(2) + sum(log(terms))
  }

  split_e_value_manual <- function(train_A, eval_A, node_j, k) {
    train_j <- train_A[, node_j]
    train_minus_j <- train_A[, -node_j, drop = FALSE]
    eval_j <- eval_A[, node_j]
    eval_minus_j <- eval_A[, -node_j, drop = FALSE]

    local_k <- k - as.integer(k > node_j)
    eval_minus_jk <- eval_minus_j[, -local_k, drop = FALSE]

    full_beta <- compute_lse(train_minus_j, train_j)
    null_beta <- compute_lse(eval_minus_jk, eval_j)

    full_log_likelihood <- compute_manual_log_likelihood(full_beta, eval_minus_j, eval_j)
    null_log_likelihood <- compute_manual_log_likelihood(null_beta, eval_minus_jk, eval_j)

    exp(full_log_likelihood - null_log_likelihood)
  }

  set.seed(2026)
  actual <- wilks_LRT_test_e_val(A)

  set.seed(2026)
  idx <- sample.int(nrow(A))
  n_half <- floor(nrow(A) / 2)
  A0 <- A[idx[1:n_half], , drop = FALSE]
  A1 <- A[idx[(n_half + 1):nrow(A)], , drop = FALSE]

  expected <- matrix(NA_real_, nrow = ncol(A), ncol = ncol(A))
  for (node_j in seq_len(ncol(A))) {
    for (k in seq_len(ncol(A))[-node_j]) {
      if (k <= node_j) {
        next
      }

      e_01 <- split_e_value_manual(A0, A1, node_j, k)
      e_10 <- split_e_value_manual(A1, A0, node_j, k)
      expected[node_j, k] <- mean(c(e_01, e_10))
    }
  }

  expect_equal(actual, expected)
})
