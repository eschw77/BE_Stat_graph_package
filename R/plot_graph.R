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


#' Plot Ising Graph from Theta Matrix
#'
#' Visualizes an Ising graph directly from a coupling matrix \code{theta}.
#' An undirected edge is drawn between nodes \eqn{i} and \eqn{j} when
#' \code{|theta[i, j]| > threshold}. Positive and negative couplings are
#' colored separately, and edge width scales with \code{|theta[i, j]|}.
#'
#' @param theta Numeric square matrix of Ising couplings.
#' @param threshold Numeric. Only edges with \code{|theta[i, j]| > threshold}
#'   are plotted. Default: 0.
#' @param layout Character or function specifying the layout algorithm. Options
#'   include "circle", "kamada.kawai", "fruchterman.reingold", and "grid".
#'   Default: "circle".
#' @param vertex.size Numeric. Size of the vertices. Default: 25.
#' @param vertex.color Character. Color of the vertices. Default: "lightblue".
#' @param edge.color.positive Character. Color for positive couplings.
#'   Default: "steelblue4".
#' @param edge.color.negative Character. Color for negative couplings.
#'   Default: "firebrick3".
#' @param edge.width.range Numeric length-2 vector giving minimum and maximum
#'   edge width used for scaling by \code{|theta|}. Default: \code{c(1, 5)}.
#' @param main Character. Main title for the plot.
#'   Default: "Ising Graph from Theta".
#' @param show.weights Logical. If TRUE, edge labels show theta values.
#'   Default: FALSE.
#'
#' @return Invisibly returns the igraph object.
#'
#' @examples
#' theta <- matrix(0, 4, 4)
#' theta[1, 2] <- theta[2, 1] <- 0.7
#' theta[2, 3] <- theta[3, 2] <- -0.5
#' theta[3, 4] <- theta[4, 3] <- 0.3
#' plot_ising_theta(theta, threshold = 0.1, show.weights = TRUE)
#'
#' @export
plot_ising_theta <- function(theta,
                             threshold = 0,
                             layout = "circle",
                             vertex.size = 25,
                             vertex.color = "lightblue",
                             edge.color.positive = "steelblue4",
                             edge.color.negative = "firebrick3",
                             edge.width.range = c(1, 5),
                             main = "Ising Graph from Theta",
                             show.weights = FALSE) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for plotting. Install it with: install.packages('igraph')")
  }

  if (!is.matrix(theta) || !is.numeric(theta) || nrow(theta) != ncol(theta)) {
    stop("theta must be a numeric square matrix.")
  }
  if (anyNA(theta)) {
    stop("theta must not contain NA values.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a single non-negative number.")
  }
  if (!is.numeric(edge.width.range) || length(edge.width.range) != 2 ||
      edge.width.range[1] <= 0 || edge.width.range[2] < edge.width.range[1]) {
    stop("edge.width.range must be length-2 with 0 < min <= max.")
  }

  n_nodes <- nrow(theta)

  # Use the upper triangle so each undirected edge appears once.
  idx <- which(upper.tri(theta) & abs(theta) > threshold, arr.ind = TRUE)

  g <- igraph::graph.empty(n = n_nodes, directed = FALSE)
  igraph::V(g)$name <- paste0("A", seq_len(n_nodes))

  edge_colors <- character(0)
  edge_widths <- numeric(0)
  edge_labels <- NA

  if (nrow(idx) > 0) {
    edge_vals <- theta[idx]
    edge_list <- as.vector(t(idx))
    g <- igraph::add_edges(g, edge_list)

    edge_colors <- ifelse(edge_vals >= 0, edge.color.positive, edge.color.negative)

    abs_vals <- abs(edge_vals)
    if (diff(range(abs_vals)) < .Machine$double.eps) {
      edge_widths <- rep(mean(edge.width.range), length(abs_vals))
    } else {
      edge_widths <- edge.width.range[1] +
        (abs_vals - min(abs_vals)) / (max(abs_vals) - min(abs_vals)) *
        (edge.width.range[2] - edge.width.range[1])
    }

    edge_labels <- if (show.weights) round(edge_vals, 3) else NA
  }

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

  plot(g,
       layout = coords,
       vertex.size = vertex.size,
       vertex.color = vertex.color,
       vertex.label = igraph::V(g)$name,
       vertex.label.color = "black",
       vertex.label.cex = 0.9,
       edge.color = if (length(edge_colors)) edge_colors else "gray80",
       edge.width = if (length(edge_widths)) edge_widths else 1,
       edge.label = edge_labels,
       edge.label.cex = 0.7,
       main = main)

  if (nrow(idx) > 0) {
    legend("bottomright",
           legend = c("Positive coupling", "Negative coupling"),
           col = c(edge.color.positive, edge.color.negative),
           lty = 1,
           lwd = 2,
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




#' Function to plot detected conditional dependence graph
#' @param avg_results Data frame with columns node_i, node_j, reject_at_alpha
#' @param threshold Numeric. Threshold for reject_at_alpha to consider an edge detected. Default is 0.5. Usually should set to 1-alpha for a specific alpha level.
#' @param main Title for the plot. Default is "Detected Graph at Threshold".
#' @return Invisibly returns the igraph object of the detected graph.
#' @export 
plot_detected_graph <- function(avg_results, threshold = 0.5, main = "Detected Graph at Threshold") {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for plotting. Install it with: install.packages('igraph')")
  }
  
  # Filter edges where reject_at_alpha >= threshold
  detected_edges <- avg_results %>%
    dplyr::filter(reject_at_alpha >= threshold) %>%
    dplyr::select(node_i, node_j)
  
  # Get number of nodes from avg_results
  n_nodes <- max(c(avg_results$node_i, avg_results$node_j))
  
  if (nrow(detected_edges) == 0) {
    # No edges detected
    g <- igraph::make_empty_graph(n = n_nodes, directed = FALSE)
  } else {
    # Create edge list and graph
    edge_list <- as.matrix(detected_edges)
    g <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    # Ensure all nodes are included even if isolated
    g <- igraph::add_vertices(g, n_nodes - igraph::vcount(g))
  }
  
  # Set vertex names
  igraph::V(g)$name <- 1:n_nodes
  
  # Plot with consistent layout
  layout <- igraph::layout_in_circle(g)
  
  plot(g,
       layout = layout,
       vertex.size = 30,
       vertex.label = igraph::V(g)$name,
       vertex.label.cex = 1.2,
       edge.width = 2,
       edge.color = "black",
       main = main)
  
  # Also print summary
  cat("Threshold:", threshold, "\n")
  cat("Detected edges:\n")
  if (nrow(detected_edges) > 0) {
    print(detected_edges)
  } else {
    cat("  No edges detected at this threshold\n")
  }
  
  invisible(g)
}


#' Plot histograms of theta values (e.g. level of CI) for specified edges across simulations
#' @param theta_list List of theta matrices (p x p) from multiple simulations
#' @param edge_list Two-column matrix of node pairs (i, j) for which to plot histograms of theta values, should typically be the ising edges input 
#' @param main_prefix Prefix for the main title of each histogram, default is "Theta Hist"
#' @return Invisibly returns NULL after plotting histograms
#' @export
plot_theta_histograms <- function(theta_list, edge_list, main_prefix = "Theta Hist") {
  if (length(theta_list) == 0) {
    return(invisible(NULL))
  }

  pair_idx <- edge_list

  n_plots <- nrow(pair_idx)
  if (n_plots == 0) {
    return(invisible(NULL))
  }

  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)
  old_par <- par(mfrow = c(n_rows, n_cols))
  on.exit(par(old_par), add = TRUE)

  for (k in 1:n_plots) {
    i <- pair_idx[k, 1]
    j <- pair_idx[k, 2]
    vals <- vapply(theta_list, function(th) th[i, j], numeric(1))
    hist(
      vals,
      breaks = 20,
      col = "gray70",
      border = "white",
      main = paste0(main_prefix, ": ", i, "-", j),
      xlab = expression(theta),
      ylab = "Frequency"
    )
  }

  invisible(NULL)
}


#' Plot Ising multinomial distribution from Theta
#'
#' Uses \code{generate_ising_to_multinomial()} to compute the multinomial
#' distribution induced by an Ising coupling matrix, then displays the
#' probabilities in a table-style plot (or bar plot) with explicit state labels.
#'
#' @param theta Coupling matrix (p x p) defining the Ising model.
#' @param bias Optional bias vector (length p). If NULL, assumed to be zero.
#' @param sort_desc Logical. If TRUE, sort rows by decreasing probability.
#'   Default: TRUE.
#' @param max_rows Optional integer. If provided, only the top \code{max_rows}
#'   states are shown in the plot.
#' @param display Character. One of \code{"table"} (default) or
#'   \code{"barplot"}.
#' @param digits Integer. Number of digits for displayed probabilities.
#'   Default: 5.
#' @param main Plot title.
#'
#' @return Invisibly returns a data frame with columns
#'   \code{configuration} and \code{probability}.
#'
#' @examples
#' theta <- matrix(0, 3, 3)
#' theta[1, 2] <- theta[2, 1] <- 0.5
#' theta[2, 3] <- theta[3, 2] <- -0.4
#' plot_ising_multinomial(theta, max_rows = 8)
#'
#' @export
plot_ising_multinomial <- function(theta,
                                   bias = NULL,
                                   sort_desc = TRUE,
                                   max_rows = NULL,
                                   display = c("table", "barplot"),
                                   digits = 5,
                                   main = "Multinomial Distribution from Ising Theta") {
  display <- match.arg(display)

  if (!is.matrix(theta) || !is.numeric(theta) || nrow(theta) != ncol(theta)) {
    stop("theta must be a numeric square matrix.")
  }
  if (!is.null(bias) && (!is.numeric(bias) || length(bias) != nrow(theta))) {
    stop("bias must be numeric with length equal to nrow(theta), or NULL.")
  }
  if (!is.null(max_rows) && (!is.numeric(max_rows) || length(max_rows) != 1 || max_rows < 1)) {
    stop("max_rows must be NULL or a positive integer.")
  }

  dist_obj <- generate_ising_to_multinomial(theta, bias)
  states <- dist_obj$states
  probs <- dist_obj$probabilities
  p <- ncol(states)

  config_labels <- apply(states, 1, function(row_vals) {
    parts <- paste0("A", seq_len(p), "=", row_vals)
    paste(parts, collapse = ", ")
  })

  out <- data.frame(
    configuration = config_labels,
    probability = probs,
    stringsAsFactors = FALSE
  )

  if (sort_desc) {
    out <- out[order(out$probability, decreasing = TRUE), , drop = FALSE]
  }

  plot_df <- out
  if (!is.null(max_rows)) {
    plot_df <- head(plot_df, as.integer(max_rows))
  }

  if (display == "barplot") {
    graphics::barplot(
      height = plot_df$probability,
      names.arg = plot_df$configuration,
      las = 2,
      cex.names = 0.7,
      col = "gray70",
      border = "white",
      ylab = "Probability",
      main = main
    )
  } else {
    n_rows <- nrow(plot_df)
    graphics::plot.new()
    graphics::title(main = main)

    # Header row
    top_y <- 0.92
    row_height <- 0.8 / (n_rows + 1)
    x_conf <- 0.02
    x_prob <- 0.98

    graphics::text(x_conf, top_y, "Configuration", adj = c(0, 0.5), font = 2)
    graphics::text(x_prob, top_y, "Probability", adj = c(1, 0.5), font = 2)
    graphics::segments(0.02, top_y - row_height / 2, 0.98, top_y - row_height / 2)

    for (i in seq_len(n_rows)) {
      y <- top_y - i * row_height
      graphics::text(x_conf, y, plot_df$configuration[i], adj = c(0, 0.5), cex = 0.85)
      graphics::text(
        x_prob,
        y,
        format(round(plot_df$probability[i], digits = digits), nsmall = digits),
        adj = c(1, 0.5),
        cex = 0.85
      )
    }
  }

  invisible(out)
}

