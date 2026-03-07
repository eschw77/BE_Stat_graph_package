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
