# Graph construction, dimensionality reduction, node classification.
# Port of spqrp/spqrp/clustering_helpers.py.

# Build the kNN graph + 2D coords. Returns list(graph, coords_2d, sample_names).
# `mds_backend` = "cmdscale" | "smacof" (smacof requires the optional dep).
create_graph_based_on_reduction_method <- function(method,
                                                    dist_matrix,
                                                    sample_names,
                                                    n_umap_neighbors = 15L,
                                                    random_state = 42L,
                                                    precomputed_graph = NULL,
                                                    n_neighbors = 1L,
                                                    mds_backend = c("cmdscale", "smacof")) {
  method <- toupper(method)
  mds_backend <- match.arg(mds_backend)
  n_samples <- length(sample_names)
  d <- stats::as.dist(dist_matrix)

  coords_2d <- switch(
    method,
    "PCA" = {
      # Classical MDS via double-centred -d^2/2 -> top-2 eigen, equivalent
      # to the J %*% (-0.5 * d^2) %*% J construction in core.py.
      stats::cmdscale(d, k = 2L)
    },
    "UMAP" = {
      rlang::check_installed("uwot")
      # uwot's internal kNN routine (`dist_nn`) aborts with
      # "1:k argument of length 0" on tied or near-tied distances, which
      # happens often on small/synthetic cohorts. Python's `umap-learn`
      # uses sklearn's NearestNeighbors which handles ties; we replicate
      # that here by computing the kNN ourselves with base R's `order()`
      # (stable; falls back to index for ties) and passing it via
      # `nn_method = list(idx, dist)`. That skips uwot's brittle path
      # entirely.
      d_mat <- as.matrix(d)
      k_eff <- min(n_umap_neighbors, n_samples - 1L)
      # uwot's nn_method format wants self as the first neighbour, then
      # k_eff actual nearest neighbours, so total k_eff + 1 columns.
      idx  <- matrix(0L,  n_samples, k_eff + 1L)
      dnn  <- matrix(0.0, n_samples, k_eff + 1L)
      for (i in seq_len(n_samples)) {
        ord <- order(d_mat[i, ])
        # First entry is always self (distance 0); subsequent are NN.
        # If self isn't at position 1 due to ties, force it there.
        ord <- c(i, ord[ord != i])
        keep <- ord[seq_len(k_eff + 1L)]
        idx[i, ] <- keep
        dnn[i, ] <- d_mat[i, keep]
      }
      # Seed uwot's internal RNG locally via withr::with_seed instead of
      # calling set.seed() directly -- mutating the user's global RNG
      # state would violate CRAN policy.
      withr::with_seed(random_state, {
        uwot::umap(
          X = NULL,
          nn_method = list(idx = idx, dist = dnn),
          n_components = 2L
        )
      })
    },
    "MDS" = {
      if (mds_backend == "smacof") {
        rlang::check_installed("smacof")
        smacof::smacofSym(d, ndim = 2L, init = "torgerson")$conf
      } else {
        stats::cmdscale(d, k = 2L)
      }
    },
    cli::cli_abort("Unknown method {.val {method}}; use one of PCA/UMAP/MDS.")
  )

  G <- igraph::make_empty_graph(n = n_samples, directed = FALSE)
  G <- igraph::set_vertex_attr(G, "name", value = sample_names)

  if (!is.null(precomputed_graph)) {
    pg_edges <- igraph::as_edgelist(precomputed_graph)
    pg_edges <- pg_edges[pg_edges[, 1L] %in% sample_names &
                           pg_edges[, 2L] %in% sample_names, , drop = FALSE]
    if (nrow(pg_edges) > 0L) {
      weights <- vapply(seq_len(nrow(pg_edges)), function(k) {
        u <- pg_edges[k, 1L]; v <- pg_edges[k, 2L]
        dist_matrix[u, v]
      }, numeric(1L))
      edge_pairs <- as.vector(t(pg_edges))
      G <- igraph::add_edges(G, edge_pairs, attr = list(weight = weights))
    }
  } else {
    edges_acc <- character()
    weights_acc <- numeric()
    seen <- character()
    for (i in seq_len(n_samples)) {
      ord <- order(dist_matrix[i, ])
      take <- ord[2:(n_neighbors + 1L)]
      for (j in take) {
        pair_key <- paste(sort(c(sample_names[i], sample_names[j])), collapse = "||")
        if (pair_key %in% seen) next
        seen <- c(seen, pair_key)
        edges_acc <- c(edges_acc, sample_names[i], sample_names[j])
        weights_acc <- c(weights_acc, dist_matrix[i, j])
      }
    }
    if (length(edges_acc) > 0L) {
      G <- igraph::add_edges(G, edges_acc, attr = list(weight = weights_acc))
    }
  }

  list(G = G, coords_2d = coords_2d, sample_names = sample_names)
}

# Iteratively remove the largest-weight edge inside any component whose
# size exceeds `max_size`, until all components are <= max_size.
split_big_component_edges_by_weight <- function(G, max_size) {
  G <- G
  repeat {
    comps <- igraph::components(G)
    big <- which(comps$csize > max_size)
    if (length(big) == 0L) break
    cid <- big[1L]
    members <- which(comps$membership == cid)
    subg <- igraph::induced_subgraph(G, vids = members)
    eweights <- igraph::E(subg)$weight
    if (length(eweights) == 0L) break
    max_idx <- which.max(eweights)
    edge_endpoints <- igraph::ends(subg, igraph::E(subg)[max_idx])
    edge_id <- igraph::get_edge_ids(G, c(edge_endpoints[1L, 1L], edge_endpoints[1L, 2L]))
    if (length(edge_id) == 0L || edge_id == 0L) break
    G <- igraph::delete_edges(G, edge_id)
  }
  G
}

# Classify nodes: TP cluster members, nodes in FP edges (cross-patient),
# singletons (alone for their patient, isolated in graph), and "isolated"
# nodes (no edge to a same-patient sample even though such a sample exists).
identify_clusters_singletons <- function(G, sample_to_patient,
                                          samples_by_patient,
                                          drawn_pairs, sample_names) {
  nodes_in_fp <- character()
  for (k in seq_along(drawn_pairs)) {
    pair <- drawn_pairs[[k]]
    if (sample_to_patient[pair[1L]] != sample_to_patient[pair[2L]]) {
      nodes_in_fp <- c(nodes_in_fp, pair[1L], pair[2L])
    }
  }
  nodes_in_fp <- unique(nodes_in_fp)

  singleton_nodes <- vapply(sample_names, function(s) {
    length(samples_by_patient[[sample_to_patient[[s]]]]) == 1L
  }, logical(1L))
  singleton_nodes <- sample_names[singleton_nodes]
  degs <- igraph::degree(G)
  singleton_graph <- intersect(singleton_nodes,
                                names(degs)[degs == 0L])

  nodes_in_tp <- character()
  comps <- igraph::components(G)
  for (cid in seq_len(comps$no)) {
    members <- igraph::V(G)$name[comps$membership == cid]
    if (length(members) < 2L) next
    patients <- unique(sample_to_patient[members])
    if (length(patients) == 1L) {
      pid <- patients[1L]
      if (length(samples_by_patient[[pid]]) >= 2L) {
        nodes_in_tp <- c(nodes_in_tp, members)
      }
    }
  }
  nodes_in_tp <- unique(nodes_in_tp)

  isolated_nodes <- character()
  for (s in sample_names) {
    if (s %in% singleton_nodes ||
        s %in% nodes_in_fp ||
        s %in% nodes_in_tp) next
    pid <- sample_to_patient[[s]]
    mates <- setdiff(samples_by_patient[[pid]], s)
    if (length(mates) == 0L) next
    has_mate_edge <- any(vapply(mates, function(m) {
      m %in% igraph::V(G)$name &&
        length(igraph::get_edge_ids(G, c(s, m))) > 0L &&
        igraph::get_edge_ids(G, c(s, m)) != 0L
    }, logical(1L)))
    if (!has_mate_edge) isolated_nodes <- c(isolated_nodes, s)
  }

  list(
    nodes_in_tp_clusters     = nodes_in_tp,
    nodes_in_fp_cluster      = nodes_in_fp,
    singleton_nodes          = singleton_graph,
    isolated_nodes           = unique(isolated_nodes)
  )
}

# Per-patient convex hull or single-edge layer for ggplot.
tp_cluster_hull_data <- function(samples_by_patient,
                                  connected_samples,
                                  drawn_pairs,
                                  sample_index,
                                  coords_2d) {
  hulls <- list()
  edges <- list()
  for (pid in names(samples_by_patient)) {
    sample_list <- samples_by_patient[[pid]]
    filtered <- sample_list[
      sample_list %in% connected_samples &
        vapply(sample_list, function(s) {
          any(vapply(setdiff(sample_list, s), function(other) {
            any(vapply(drawn_pairs, function(p) {
              setequal(p, c(s, other))
            }, logical(1L)))
          }, logical(1L)))
        }, logical(1L))
    ]
    if (length(filtered) < 2L) next
    if (length(filtered) == 2L) {
      i <- sample_index[[filtered[1L]]]; j <- sample_index[[filtered[2L]]]
      edges[[length(edges) + 1L]] <- data.frame(
        x = coords_2d[i, 1L], y = coords_2d[i, 2L],
        xend = coords_2d[j, 1L], yend = coords_2d[j, 2L],
        patient = pid
      )
      next
    }
    coords_sub <- coords_2d[match(filtered, names(sample_index)), , drop = FALSE]
    hull_idx <- tryCatch(grDevices::chull(coords_sub[, 1L], coords_sub[, 2L]),
                          error = function(e) integer())
    if (length(hull_idx) < 3L) next
    hulls[[length(hulls) + 1L]] <- data.frame(
      x = coords_sub[hull_idx, 1L],
      y = coords_sub[hull_idx, 2L],
      patient = pid,
      stringsAsFactors = FALSE
    )
  }
  list(hulls = hulls, edges = edges)
}

# Hand-rolled ARI and NMI on two clusterings (each a named integer vector
# of cluster IDs by sample). Matches sklearn semantics.
adjusted_rand_index <- function(true_labels, pred_labels) {
  tab <- table(true_labels, pred_labels)
  n <- sum(tab)
  if (n < 2L) return(0)
  sum_comb_c <- sum(choose(rowSums(tab), 2L))
  sum_comb_k <- sum(choose(colSums(tab), 2L))
  sum_comb   <- sum(choose(tab,         2L))
  expected   <- sum_comb_c * sum_comb_k / choose(n, 2L)
  max_index  <- 0.5 * (sum_comb_c + sum_comb_k)
  if (max_index == expected) return(0)
  (sum_comb - expected) / (max_index - expected)
}

normalized_mutual_info <- function(true_labels, pred_labels) {
  tab <- table(true_labels, pred_labels)
  n <- sum(tab)
  if (n < 2L) return(0)
  pxy <- tab / n
  px  <- rowSums(pxy)
  py  <- colSums(pxy)
  mi <- 0
  for (i in seq_along(px)) {
    for (j in seq_along(py)) {
      if (pxy[i, j] > 0) {
        mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx + hy == 0) return(0)
  2 * mi / (hx + hy)
}

# Pairwise classification metrics computed from the transitive closure of
# `drawn_pairs` (i.e. samples in the same connected component of G are
# predicted to share a cluster).
transitive_performance <- function(sample_names, drawn_pairs, sample_to_patient,
                                    quiet = TRUE) {
  G <- igraph::make_empty_graph(n = length(sample_names), directed = FALSE)
  G <- igraph::set_vertex_attr(G, "name", value = sample_names)
  if (length(drawn_pairs) > 0L) {
    edge_vec <- unlist(drawn_pairs, use.names = FALSE)
    G <- igraph::add_edges(G, edge_vec)
  }

  comps <- igraph::components(G)
  pred_cluster <- comps$membership
  names(pred_cluster) <- igraph::V(G)$name
  true_cluster <- as.integer(factor(sample_to_patient[sample_names]))
  names(true_cluster) <- sample_names

  combos <- utils::combn(sample_names, 2L)
  TP <- FP <- FN <- TN <- 0L
  false_negatives <- list()
  for (k in seq_len(ncol(combos))) {
    s1 <- combos[1L, k]; s2 <- combos[2L, k]
    same_patient <- sample_to_patient[[s1]] == sample_to_patient[[s2]]
    same_cluster <- pred_cluster[[s1]] == pred_cluster[[s2]]
    if (same_patient && same_cluster) TP <- TP + 1L
    else if (same_patient && !same_cluster) {
      FN <- FN + 1L
      false_negatives[[length(false_negatives) + 1L]] <- c(s1, s2)
    }
    else if (!same_patient && same_cluster) FP <- FP + 1L
    else TN <- TN + 1L
  }
  total <- TP + FP + FN + TN

  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1 <- if ((precision + sensitivity) > 0)
    2 * precision * sensitivity / (precision + sensitivity) else 0
  accuracy <- if (total > 0) (TP + TN) / total else 0
  bACC <- if ((TP + FN) > 0 && (TN + FP) > 0)
    0.5 * (TP / (TP + FN) + TN / (TN + FP)) else 0

  ari <- adjusted_rand_index(true_cluster[sample_names], pred_cluster[sample_names])
  nmi <- normalized_mutual_info(true_cluster[sample_names], pred_cluster[sample_names])

  if (!quiet) {
    cli::cli_h2("Pairwise Clustering Performance (Transitive)")
    cli::cli_text("TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    cli::cli_text("Precision:         {sprintf('%.3f', precision)}")
    cli::cli_text("Sensitivity:       {sprintf('%.3f', sensitivity)}")
    cli::cli_text("F1 Score:          {sprintf('%.3f', f1)}")
    cli::cli_text("Accuracy:          {sprintf('%.3f', accuracy)}")
    cli::cli_text("Balanced Accuracy: {sprintf('%.3f', bACC)}")

    cli::cli_h2("Overall Clustering Agreement")
    cli::cli_text("Adjusted Rand Index (ARI):    {sprintf('%.3f', ari)}")
    cli::cli_text("Normalized Mutual Info (NMI): {sprintf('%.3f', nmi)}")

    if (length(false_negatives) > 0L) {
      cli::cli_h3(
        "False Negative (FN) pairs - same patient but not transitively connected: {length(false_negatives)}"
      )
      for (pair in false_negatives) {
        cli::cli_text("  - {pair[1L]} <-> {pair[2L]}")
      }
    } else {
      cli::cli_text("No False Negative (FN) pairs found.")
    }
  }

  list(
    precision = precision,
    sensitivity = sensitivity,
    f1 = f1,
    accuracy = accuracy,
    `balanced accuracy` = bACC,
    ari = ari,
    nmi = nmi,
    false_negatives = false_negatives
  )
}
