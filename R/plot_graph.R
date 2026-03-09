# R/plot_graph.R

#' Plot True Ising Model Graph Structure
#'
#' Visualizes the true graph structure used to generate Ising model data.
#' Edges represent conditional dependencies (non-independence) between nodes.
#'
#' @param result   List returned by \code{\link{ising_generate}}, containing the
#'                 true edge structure and theta coupling matrix.
#' @param layout   Character or function specifying the layout algorithm. Options
#'                 include "circle", "kamada.kawai", "fruchterman.reingold", etc.
#'                 Default: "circle".
#' @param vertex.size Numeric. Size of the vertices. Default: 25.
#' @param vertex.color Character. Color of the vertices. Default: "lightblue".
#' @param edge.color Character. Color of the edges. Default: "gray50".
#' @param edge.width Numeric. Width of the edges. Default: 2.
#' @param weak.edge.color Character. Color of the weak edges. Default: "gray80".
#' @param weak.edge.lty Numeric. Line type for weak edges (1=solid, 2=dashed, etc.). Default: 2.
#' @param main Character. Main title for the plot. Default: "True Ising Model Graph".
#' @param show.weights Logical. If TRUE, edge labels show the theta coupling values.
#'                     Default: FALSE.
#'
#' @return Invisibly returns the igraph object.
#'
#' @examples
#' edges <- matrix(c(1,2, 2,3, 3,4, 4,5, 1,5), ncol = 2, byrow = TRUE)
#' result <- ising_generate(n_nodes = 5, edges = edges, n_samples = 1000, seed = 42)
#' plot_ising_true(result)
#' plot_ising_true(result, layout = "kamada.kawai", show.weights = TRUE)
#'
#' @export
plot_ising_true <- function(result,
                            layout = "circle",
                            vertex.size = 25,
                            vertex.color = "lightblue",
                            edge.color = "gray50",
                            edge.width = 2,
                            weak.edge.color = "gray80",
                            weak.edge.lty = 2,
                            main = "True Ising Model Graph",
                            show.weights = FALSE) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for plotting. Install it with: install.packages('igraph')")
  }
  
  # Extract graph information
  n_nodes <- ncol(result$samples)
  edges <- result$edges
  weak_edges <- result$weak_edges
  theta <- result$theta
  
  # Create igraph object
  g <- igraph::graph.empty(n = n_nodes, directed = FALSE)
  igraph::V(g)$name <- paste0("A", seq_len(n_nodes))
  
  # Add strong edges
  n_strong <- 0
  if (!is.null(edges) && nrow(edges) > 0) {
    edge_list <- as.vector(t(edges))
    g <- igraph::add_edges(g, edge_list)
    n_strong <- nrow(edges)
  }
  
  # Add weak edges
  n_weak <- 0
  if (!is.null(weak_edges) && nrow(weak_edges) > 0) {
    weak_edge_list <- as.vector(t(weak_edges))
    g <- igraph::add_edges(g, weak_edge_list)
    n_weak <- nrow(weak_edges)
  }
  
  # Extract edge weights (theta values) for all edges
  total_edges <- n_strong + n_weak
  edge_weights <- numeric(total_edges)
  edge_colors <- character(total_edges)
  edge_ltys <- numeric(total_edges)
  edge_widths <- numeric(total_edges)
  
  # Strong edges
  if (n_strong > 0) {
    for (k in seq_len(n_strong)) {
      edge_weights[k] <- theta[edges[k, 1], edges[k, 2]]
      edge_colors[k] <- edge.color
      edge_ltys[k] <- 1
      edge_widths[k] <- edge.width
    }
  }
  
  # Weak edges
  if (n_weak > 0) {
    for (k in seq_len(n_weak)) {
      idx <- n_strong + k
      edge_weights[idx] <- theta[weak_edges[k, 1], weak_edges[k, 2]]
      edge_colors[idx] <- weak.edge.color
      edge_ltys[idx] <- weak.edge.lty
      edge_widths[idx] <- edge.width * 0.7
    }
  }
  
  # Set layout
  if (is.character(layout)) {
    layout_func <- switch(layout,
      "circle" = igraph::layout_in_circle,
      "kamada.kawai" = igraph::layout_with_kk,
      "fruchterman.reingold" = igraph::layout_with_fr,
      "grid" = igraph::layout_on_grid,
      igraph::layout_in_circle
    )
    coords <- layout_func(g)
  } else {
    coords <- layout(g)
  }
  
  # Prepare edge labels
  edge_labels <- if (show.weights && total_edges > 0) round(edge_weights, 3) else NA
  
  # Plot
  plot(g,
       layout = coords,
       vertex.size = vertex.size,
       vertex.color = vertex.color,
       vertex.label = igraph::V(g)$name,
       vertex.label.color = "black",
       vertex.label.cex = 0.9,
       edge.color = if (total_edges > 0) edge_colors else edge.color,
       edge.width = if (total_edges > 0) edge_widths else edge.width,
       edge.lty = if (total_edges > 0) edge_ltys else 1,
       edge.label = edge_labels,
       edge.label.cex = 0.7,
       main = main)
  
  # Add legend if weak edges present
  if (n_weak > 0) {
    legend("bottomright",
           legend = c("Strong edges", "Weak edges"),
           col = c(edge.color, weak.edge.color),
           lty = c(1, weak.edge.lty),
           lwd = c(edge.width, edge.width * 0.7),
           cex = 0.8,
           bg = "white")
  }
  
  invisible(g)
}


#' Plot Estimated Graph from Conditional Independence Testing
#'
#' Visualizes the graph structure estimated from data using conditional independence
#' tests. Edges are drawn between nodes that show significant conditional dependence
#' based on the mean log-odds ratio from \code{\link{ising_ci_test}}.
#'
#' @param samples Integer matrix of ±1 values (n_samples × n_nodes), typically from
#'                \code{result$samples}.
#' @param edges   Optional. Two-column matrix of true edges for comparison. If provided,
#'                edges will be colored to distinguish true positives from false positives.
#' @param weak_edges Optional. Two-column matrix of true weak edges for comparison.
#' @param threshold Numeric. Threshold for |mean_log_OR| above which an edge is drawn.
#'                  Higher values = more conservative (fewer edges). Default: 0.1.
#' @param min_stratum_size Minimum stratum count for CI testing. Passed to
#'                         \code{\link{ising_ci_test}}. Default: 10.
#' @param layout   Character or function specifying the layout algorithm.
#'                 Default: "circle".
#' @param vertex.size Numeric. Size of the vertices. Default: 25.
#' @param vertex.color Character. Color of the vertices. Default: "lightgreen".
#' @param edge.color.correct Character. Color for correctly identified edges (true positives).
#'                           Only used if \code{edges} is provided. Default: "darkgreen".
#' @param edge.color.weak Character. Color for correctly identified weak edges.
#'                        Only used if \code{weak_edges} is provided. Default: "orange".
#' @param edge.color.false Character. Color for incorrectly identified edges (false positives).
#'                         Only used if \code{edges} is provided. Default: "red".
#' @param edge.color.default Character. Color for edges when true edges not provided. Default: "gray50".
#' @param edge.width Numeric. Width of the edges. Default: 2.
#' @param main Character. Main title for the plot. Default: "Estimated Graph from CI Tests".
#' @param show.log_OR Logical. If TRUE, edge labels show the mean log-odds ratio values.
#'                    Default: FALSE.
#'
#' @return Invisibly returns a list with the igraph object and the CI test results.
#'
#' @examples
#' edges <- matrix(c(1,2, 2,3, 3,4, 4,5, 1,5), ncol = 2, byrow = TRUE)
#' result <- ising_generate(n_nodes = 5, edges = edges, n_samples = 1000, seed = 42)
#' 
#' # Plot estimated graph
#' plot_ising_estimated(result$samples)
#' 
#' # Compare with true edges
#' plot_ising_estimated(result$samples, edges = edges, threshold = 0.15)
#' 
#' # More conservative threshold
#' plot_ising_estimated(result$samples, edges = edges, threshold = 0.3)
#'
#' @export
plot_ising_estimated <- function(samples,
                                 edges = NULL,
                                 weak_edges = NULL,
                                 threshold = 0.1,
                                 min_stratum_size = 10,
                                 layout = "circle",
                                 vertex.size = 25,
                                 vertex.color = "lightgreen",
                                 edge.color.correct = "darkgreen",
                                 edge.color.weak = "orange",
                                 edge.color.false = "red",
                                 edge.color.default = "gray50",
                                 edge.width = 2,
                                 main = "Estimated Graph from CI Tests",
                                 show.log_OR = FALSE) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for plotting. Install it with: install.packages('igraph')")
  }
  
  # Run CI test on all pairs
  ci_results <- ising_ci_test(samples, edges = if (is.null(edges)) matrix(0, 0, 2) else edges, 
                              min_stratum_size = min_stratum_size)
  
  # Identify edges based on threshold
  ci_results$estimated_edge <- abs(ci_results$mean_log_OR) > threshold & !is.na(ci_results$mean_log_OR)
  
  # Extract node information
  n_nodes <- ncol(samples)
  
  # Create igraph object
  g <- igraph::graph.empty(n = n_nodes, directed = FALSE)
  igraph::V(g)$name <- paste0("A", seq_len(n_nodes))
  
  # Filter to estimated edges
  estimated <- ci_results[ci_results$estimated_edge, ]
  
  if (nrow(estimated) > 0) {
    # Extract node numbers from names like "A1", "A2"
    node_i <- as.integer(gsub("A", "", estimated$node_i))
    node_j <- as.integer(gsub("A", "", estimated$node_j))
    
    # Add edges
    edge_list <- as.vector(rbind(node_i, node_j))
    g <- igraph::add_edges(g, edge_list)
    
    # Determine edge colors if true edges provided
    if (!is.null(edges) && nrow(edges) > 0) {
      edge_set <- apply(edges, 1, function(r) paste(sort(r), collapse = "-"))
      weak_set <- if (!is.null(weak_edges) && nrow(weak_edges) > 0) {
        apply(weak_edges, 1, function(r) paste(sort(r), collapse = "-"))
      } else {
        character(0)
      }
      
      estimated_pairs <- paste(pmin(node_i, node_j), pmax(node_i, node_j), sep = "-")
      is_strong <- estimated_pairs %in% edge_set
      is_weak <- estimated_pairs %in% weak_set
      
      edge_colors <- ifelse(is_strong, edge.color.correct,
                     ifelse(is_weak, edge.color.weak, edge.color.false))
    } else {
      edge_colors <- rep(edge.color.default, nrow(estimated))
    }
    
    # Prepare edge labels
    edge_labels <- if (show.log_OR) round(estimated$mean_log_OR, 3) else NA
  } else {
    edge_colors <- character(0)
    edge_labels <- NA
  }
  
  # Set layout
  if (is.character(layout)) {
    layout_func <- switch(layout,
      "circle" = igraph::layout_in_circle,
      "kamada.kawai" = igraph::layout_with_kk,
      "fruchterman.reingold" = igraph::layout_with_fr,
      "grid" = igraph::layout_on_grid,
      igraph::layout_in_circle
    )
    coords <- layout_func(g)
  } else {
    coords <- layout(g)
  }
  
  # Plot
  plot(g,
       layout = coords,
       vertex.size = vertex.size,
       vertex.color = vertex.color,
       vertex.label = igraph::V(g)$name,
       vertex.label.color = "black",
       vertex.label.cex = 0.9,
       edge.color = edge_colors,
       edge.width = edge.width,
       edge.label = edge_labels,
       edge.label.cex = 0.7,
       main = main)
  
  # Add legend if comparing with true edges
  if (!is.null(edges) && nrow(edges) > 0 && nrow(estimated) > 0) {
    if (!is.null(weak_edges) && nrow(weak_edges) > 0) {
      legend("bottomright",
             legend = c("Strong edge (TP)", "Weak edge (TP)", "False Positive"),
             col = c(edge.color.correct, edge.color.weak, edge.color.false),
             lwd = edge.width,
             cex = 0.8,
             bg = "white")
    } else {
      legend("bottomright",
             legend = c("True Positive", "False Positive"),
             col = c(edge.color.correct, edge.color.false),
             lwd = edge.width,
             cex = 0.8,
             bg = "white")
    }
  }
  
  invisible(list(graph = g, ci_results = ci_results))
}
