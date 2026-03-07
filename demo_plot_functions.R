# Demo: Plotting Ising Model Graphs
# This script demonstrates the two new plotting functions:
# 1. plot_ising_true() - plots the true graph structure
# 2. plot_ising_estimated() - plots estimated graph from CI tests

# Install igraph if needed
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

# Load the package (or source functions during development)
library(BEstatGraphTools)

# ── Example 1: Simple 5-node graph ────────────────────────────────────────────

# Define a pentagon graph structure
edges <- matrix(c(1,2, 2,3, 3,4, 4,5, 1,5), ncol = 2, byrow = TRUE)
probs <- c(0.7, 0.4, 0.6, 0.3, 0.8)

# Generate data from Ising model
set.seed(42)
result <- ising_generate(
  n_nodes      = 5,
  edges        = edges,
  target_probs = probs,
  n_samples    = 2000,
  seed         = 42,
  verbose      = TRUE
)

# Summary of the generated data
ising_summary(result)

# Plot 1: True graph structure
par(mfrow = c(1, 2))
plot_ising_true(result, 
                main = "True Graph Structure",
                show.weights = TRUE,
                layout = "circle")

# Plot 2: Estimated graph from CI tests
plot_ising_estimated(result$samples, 
                     edges = edges,  # compare with true edges
                     threshold = 0.15,
                     main = "Estimated from CI Tests")

# ── Example 2: More complex graph with different thresholds ──────────────────

# Define a more complex graph
edges2 <- matrix(c(
  1,2, 1,3, 1,4,
  2,3, 2,5,
  3,4, 3,5,
  4,5
), ncol = 2, byrow = TRUE)

# Generate data
set.seed(123)
result2 <- ising_generate(
  n_nodes   = 5,
  edges     = edges2,
  n_samples = 3000,
  seed      = 123
)

# Compare different thresholds
par(mfrow = c(2, 2))

plot_ising_true(result2, 
                main = "True Graph", 
                layout = "kamada.kawai")

plot_ising_estimated(result2$samples, 
                     edges = edges2,
                     threshold = 0.05,
                     main = "Estimated (threshold = 0.05)",
                     layout = "kamada.kawai")

plot_ising_estimated(result2$samples, 
                     edges = edges2,
                     threshold = 0.15,
                     main = "Estimated (threshold = 0.15)",
                     layout = "kamada.kawai")

plot_ising_estimated(result2$samples, 
                     edges = edges2,
                     threshold = 0.3,
                     main = "Estimated (threshold = 0.3)",
                     layout = "kamada.kawai")

# ── Example 3: Examine CI test results ──────────────────────────────────────

# Get detailed CI test results
ci_output <- plot_ising_estimated(result$samples, 
                                   edges = edges,
                                   threshold = 0.15)

# View the CI test results
print(ci_output$ci_results)

# Summary statistics
cat("\nSummary of |mean_log_OR| for edges vs non-edges:\n")
with(ci_output$ci_results, {
  cat("True edges:     mean =", mean(abs(mean_log_OR[is_edge]), na.rm = TRUE), "\n")
  cat("Non-edges (CI): mean =", mean(abs(mean_log_OR[!is_edge]), na.rm = TRUE), "\n")
})
