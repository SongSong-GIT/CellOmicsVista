.libPaths(c("D:/R/R-4.3.1/library/seuratV5", .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(StabMap)
  library(BiocNeighbors)
  library(harmony)
  library(uwot)
  library(igraph)
  library(ggplot2)
  library(scatterplot3d)
})
set.seed(2026)

args <- commandArgs(trailingOnly = TRUE)
destination_root <- if (length(args) >= 1L) args[[1L]] else
  "G:/multiomics_integration_data"
manifest_file <- if (length(args) >= 2L) args[[2L]] else
  file.path(destination_root, "disease_manifests", "batch_manifest.csv")
max_disease_cells <- as.integer(Sys.getenv("MAX_CELLS_PER_DISEASE", "20000"))
selected_diseases <- strsplit(Sys.getenv("DISEASES", ""), ",", fixed = TRUE)[[1L]]
selected_diseases <- trimws(selected_diseases[nzchar(selected_diseases)])
force_diseases <- strsplit(Sys.getenv("FORCE_DISEASES", ""), ",", fixed = TRUE)[[1L]]
force_diseases <- trimws(force_diseases[nzchar(force_diseases)])

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_argument)) sub("^--file=", "", file_argument[[1L]]) else
  sys.frame(1)$ofile
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "multiomics_io.R"), encoding = "UTF-8")

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE,
                     fileEncoding = "UTF-8-BOM")
if (!"dataset_key" %in% names(manifest))
  manifest$dataset_key <- paste(manifest$accession, manifest$omics, sep = "__")
manifest <- manifest[file.exists(manifest$output_file), , drop = FALSE]
if (length(selected_diseases))
  manifest <- manifest[manifest$disease %in% selected_diseases, , drop = FALSE]

slugify <- function(x) {
  value <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", value)
}

row_rms_scale <- function(x) {
  rms <- sqrt(pmax(Matrix::rowMeans(x^2), 1e-8))
  answer <- as_dgc(Matrix::Diagonal(x = 1 / rms) %*% x)
  dimnames(answer) <- dimnames(x)
  answer
}

make_link_matrix <- function(links, from_column, features, genes) {
  if (is.null(links) || !nrow(links))
    return(Matrix::sparseMatrix(i = integer(), j = integer(),
                                dims = c(length(genes), length(features)),
                                dimnames = list(genes, features)))
  i <- match(toupper(links$gene), genes)
  j <- match(links[[from_column]], features)
  weight <- if ("weight" %in% names(links)) links$weight else rep(1, nrow(links))
  keep <- !is.na(i) & !is.na(j) & is.finite(weight) & weight > 0
  answer <- Matrix::sparseMatrix(i = i[keep], j = j[keep], x = weight[keep],
                                 dims = c(length(genes), length(features)),
                                 dimnames = list(genes, features))
  answer <- Matrix::Diagonal(x = 1 / pmax(Matrix::rowSums(answer), 1)) %*% answer
  as_dgc(answer)
}

knn_metrics <- function(embedding, metadata, k = 30L) {
  neighbors <- BiocNeighbors::findKNN(
    embedding, k = min(k, nrow(embedding) - 1L), get.distance = FALSE,
    BNPARAM = BiocNeighbors::AnnoyParam(ntrees = 50, search.mult = 100,
                                       distance = "Euclidean"))$index
  batches <- as.character(metadata$dataset)
  neighbor_batches <- matrix(batches[as.vector(neighbors)], nrow = nrow(neighbors))
  observed <- mean(neighbor_batches != batches)
  expected <- 1 - sum(prop.table(table(batches))^2)
  labels <- as.character(metadata$broad_celltype)
  neighbor_labels <- matrix(labels[as.vector(neighbors)], nrow = nrow(neighbors))
  valid <- labels != "Unassigned"
  consistency <- if (any(valid))
    mean(rowMeans(neighbor_labels[valid, , drop = FALSE] == labels[valid])) else NA_real_
  list(index = neighbors, neighbor_batches = neighbor_batches, observed = observed,
       expected = expected, normalized = observed / expected,
       celltype_consistency = consistency)
}

save_png <- function(path, plot, width = 9, height = 6.5) {
  grDevices::png(path, width = width, height = height, units = "in", res = 300,
                 type = "cairo")
  print(plot)
  invisible(grDevices::dev.off())
}

save_umap_3d_plot <- function(embedding, groups, path, title, legend_title) {
  groups <- as.character(groups)
  groups[is.na(groups) | !nzchar(groups)] <- "Unknown"
  levels <- sort(unique(groups))
  palette <- setNames(grDevices::hcl.colors(length(levels), "Dark 3"), levels)
  point_colors <- unname(palette[groups])
  show_legend <- length(levels) <= 24L
  grDevices::png(path, width = 2100, height = 1600, res = 200)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 3.5, if (show_legend) 10 else 2), xpd = TRUE)
  scatterplot3d::scatterplot3d(
    embedding[, 1L], embedding[, 2L], embedding[, 3L],
    color = point_colors, pch = 16, cex.symbols = 0.45,
    angle = 55, grid = TRUE, box = FALSE,
    xlab = "UMAP 1", ylab = "UMAP 2", zlab = "UMAP 3",
    main = title
  )
  if (show_legend) {
    graphics::legend("topright", inset = c(0.01, 0.01), legend = levels,
                     col = unname(palette), pch = 16, cex = 0.65,
                     bty = "n", title = legend_title)
  }
  invisible(path)
}

subset_source <- function(source, cells) {
  source$gene <- source$gene[, cells, drop = FALSE]
  if (!is.null(source$protein)) source$protein <- source$protein[, cells, drop = FALSE]
  if (!is.null(source$peak)) source$peak <- source$peak[, cells, drop = FALSE]
  if (!is.null(source$coords)) source$coords <- source$coords[cells, , drop = FALSE]
  source$metadata <- source$metadata[cells, , drop = FALSE]
  source
}

results_root <- file.path(destination_root, "results")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
integration_status_file <- file.path(destination_root, "logs", "integration_status.csv")
integration_status <- if (file.exists(integration_status_file))
  read.csv(integration_status_file, stringsAsFactors = FALSE,
           fileEncoding = "UTF-8-BOM") else data.frame()

record_integration <- function(record) {
  if (nrow(integration_status))
    integration_status <<- integration_status[integration_status$disease != record$disease, ,
                                               drop = FALSE]
  integration_status <<- rbind(integration_status, record)
  write.csv(integration_status, integration_status_file, row.names = FALSE,
            fileEncoding = "UTF-8")
}

diseases <- unique(manifest$disease)
for (disease_index in seq_along(diseases)) {
  disease <- diseases[[disease_index]]
  rows <- manifest[manifest$disease == disease, , drop = FALSE]
  if (length(unique(rows$omics)) < 2L || nrow(rows) < 2L) next
  output_dir <- file.path(results_root, slugify(disease))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  complete_file <- file.path(output_dir, "complete_results.rds")
  if (file.exists(complete_file) && !disease %in% force_diseases) {
    log_message(sprintf("[%d/%d] %s already integrated", disease_index,
                        length(diseases), disease))
    next
  }
  started <- Sys.time()
  result <- tryCatch({
    log_message(sprintf("[%d/%d] START %s (%d datasets, %d omics)", disease_index,
                        length(diseases), disease, nrow(rows), length(unique(rows$omics))))
    sources <- lapply(rows$output_file, readRDS)
    names(sources) <- rows$dataset_key
    total_cells <- sum(vapply(sources, function(x) ncol(x$gene), integer(1)))
    if (total_cells > max_disease_cells) {
      quota <- max(100L, floor(max_disease_cells / length(sources)))
      sources <- lapply(sources, function(x) {
        metadata <- x$metadata
        cells <- sample_stratified(metadata$cell_id,
                                   paste(metadata$sample, metadata$broad_celltype), quota)
        subset_source(x, cells)
      })
    }

    common_genes <- Reduce(intersect, lapply(sources, function(x) rownames(x$gene)))
    if (length(common_genes) < 300L)
      stop("Only ", length(common_genes), " genes are shared by every dataset")
    # Rank shared genes by pooled nonzero variability, then fix the feature set before StabMap.
    candidate_genes <- common_genes
    pooled_score <- setNames(numeric(length(candidate_genes)), candidate_genes)
    for (source in sources) {
      x <- source$gene[candidate_genes, , drop = FALSE]
      mu <- Matrix::rowMeans(x)
      xx <- x
      xx@x <- xx@x^2
      pooled_score <- pooled_score + pmax(Matrix::rowMeans(xx) - mu^2, 0)
    }
    selected_genes <- names(head(sort(pooled_score, decreasing = TRUE), 1000L))
    total_n <- sum(vapply(sources, function(x) ncol(x$gene), integer(1)))
    gene_sum <- gene_sumsq <- setNames(numeric(length(selected_genes)), selected_genes)
    for (source in sources) {
      x <- source$gene[selected_genes, , drop = FALSE]
      gene_sum <- gene_sum + Matrix::rowSums(x)
      xx <- x
      xx@x <- xx@x^2
      gene_sumsq <- gene_sumsq + Matrix::rowSums(xx)
    }
    gene_sd <- sqrt(pmax((gene_sumsq - gene_sum^2 / total_n) / pmax(total_n - 1L, 1L),
                         1e-6))

    assays <- vector("list", length(sources))
    names(assays) <- names(sources)
    for (name in names(sources)) {
      source <- sources[[name]]
      gene <- as_dgc(Matrix::Diagonal(x = 1 / gene_sd) %*%
                       source$gene[selected_genes, , drop = FALSE])
      dimnames(gene) <- list(selected_genes, colnames(source$gene))
      extras <- list()
      if (!is.null(source$protein) && nrow(source$protein)) {
        protein <- row_rms_scale(source$protein)
        link <- make_link_matrix(source$protein_links, "protein", rownames(protein),
                                 selected_genes)
        gene <- gene + 0.08 * (link %*% protein)
        extras$protein <- 0.08 * protein
      }
      if (!is.null(source$peak) && nrow(source$peak)) {
        peak <- row_rms_scale(source$peak)
        link <- make_link_matrix(source$peak_links, "peak", rownames(peak), selected_genes)
        gene <- gene + 0.06 * (link %*% peak)
        extras$peak <- 0.06 * peak
      }
      assays[[name]] <- as_dgc(do.call(rbind, c(list(gene), extras)))
    }
    saveRDS(assays, file.path(output_dir, "01_StabMap_assay_list.rds"))
    topology <- StabMap::mosaicDataTopology(assays)
    if (!igraph::is_connected(topology)) stop("StabMap feature topology is disconnected")
    grDevices::png(file.path(output_dir, "01_mosaic_topology.png"), width = 7,
                   height = 7, units = "in", res = 300, type = "cairo")
    plot(topology, vertex.size = 30, vertex.label.cex = 0.65)
    invisible(grDevices::dev.off())

    richness <- vapply(sources, function(x)
      as.integer(!is.null(x$protein)) + as.integer(!is.null(x$peak)) + ncol(x$gene) / 1e6,
      numeric(1))
    references <- names(head(sort(richness, decreasing = TRUE), min(3L, length(sources))))
    log_message("Running actual StabMap; references:", paste(references, collapse = ", "))
    stab_raw <- StabMap::stabMap(
      assay_list = assays, reference_list = references,
      ncomponentsReference = 30, ncomponentsSubset = 30,
      suppressMessages = FALSE, projectAll = FALSE, restrictFeatures = FALSE,
      maxFeatures = max(vapply(assays, nrow, integer(1))), plot = FALSE,
      scale.center = TRUE, scale.scale = FALSE
    )
    stab_raw <- as.matrix(stab_raw)
    if (any(!is.finite(stab_raw))) stop("StabMap returned non-finite values")
    saveRDS(stab_raw, file.path(output_dir, "02_StabMap_embedding.rds"))

    metadata <- do.call(rbind, lapply(sources, `[[`, "metadata"))
    rownames(metadata) <- metadata$cell_id
    metadata <- metadata[rownames(stab_raw), , drop = FALSE]
    if (anyNA(metadata$dataset)) stop("StabMap cells are missing from metadata")
    raw_qc <- knn_metrics(stab_raw, metadata)

    log_message("Running Harmony on the StabMap embedding")
    stab_harmony <- harmony::RunHarmony(
      stab_raw, meta_data = metadata, vars_use = "dataset",
      theta = 2, lambda = 1, sigma = 0.1, max_iter = 20,
      block_size = 0.05, verbose = TRUE
    )
    stab_harmony <- as.matrix(stab_harmony)
    rownames(stab_harmony) <- rownames(stab_raw)
    saveRDS(stab_harmony, file.path(output_dir, "03_StabMap_Harmony_embedding.rds"))
    harmony_qc <- knn_metrics(stab_harmony, metadata)

    final_embedding <- stab_harmony
    spatial_datasets <- names(Filter(function(x) !is.null(x$coords), sources))
    for (name in spatial_datasets) {
      coords <- sources[[name]]$coords
      ids <- intersect(rownames(coords), rownames(final_embedding))
      coords <- coords[ids, , drop = FALSE]
      if (length(ids) < 7L) next
      index <- BiocNeighbors::findKNN(coords, k = 6, get.distance = FALSE)$index
      weights <- Matrix::sparseMatrix(i = rep(seq_len(nrow(coords)), each = ncol(index)),
                                      j = as.vector(t(index)), x = 1,
                                      dims = c(nrow(coords), nrow(coords)))
      weights <- weights + Matrix::t(weights)
      weights@x[] <- 1
      weights <- Matrix::Diagonal(x = 1 / pmax(Matrix::rowSums(weights), 1)) %*% weights
      final_embedding[ids, ] <- 0.90 * stab_harmony[ids, , drop = FALSE] +
        0.10 * as.matrix(weights %*% stab_harmony[ids, , drop = FALSE])
    }
    saveRDS(final_embedding, file.path(output_dir, "04_final_spatial_embedding.rds"))
    final_qc <- knn_metrics(final_embedding, metadata)

    metrics <- data.frame(
      stage = c("StabMap", "StabMap + Harmony", "StabMap + Harmony + spatial 6-NN"),
      observed_cross_dataset_neighbors = c(raw_qc$observed, harmony_qc$observed,
                                           final_qc$observed),
      expected_cross_dataset_neighbors = c(raw_qc$expected, harmony_qc$expected,
                                           final_qc$expected),
      normalized_mixing = c(raw_qc$normalized, harmony_qc$normalized,
                            final_qc$normalized),
      neighbor_celltype_consistency = c(raw_qc$celltype_consistency,
                                        harmony_qc$celltype_consistency,
                                        final_qc$celltype_consistency)
    )
    write.csv(metrics, file.path(output_dir, "stage_metrics.csv"), row.names = FALSE)

    neighbors <- final_qc$index
    edges <- unique(data.frame(
      from = pmin(rep(seq_len(nrow(final_embedding)), each = ncol(neighbors)),
                  as.vector(t(neighbors))),
      to = pmax(rep(seq_len(nrow(final_embedding)), each = ncol(neighbors)),
                as.vector(t(neighbors)))
    ))
    edges <- edges[edges$from != edges$to, ]
    graph <- igraph::graph_from_data_frame(edges, directed = FALSE,
      vertices = data.frame(name = as.character(seq_len(nrow(final_embedding)))))
    community <- igraph::cluster_louvain(graph)
    cluster <- factor(igraph::membership(community))
    umap <- uwot::umap(final_embedding, n_neighbors = 30, min_dist = 0.30,
                       metric = "cosine", n_components = 2, n_threads = 1,
                       verbose = TRUE)
    set.seed(2026)
    umap_3d <- uwot::umap(final_embedding, n_neighbors = 30, min_dist = 0.30,
                          metric = "cosine", n_components = 3, n_threads = 1,
                          verbose = TRUE)
    rownames(umap) <- rownames(final_embedding)
    rownames(umap_3d) <- rownames(final_embedding)
    metadata$UMAP1 <- umap[, 1L]
    metadata$UMAP2 <- umap[, 2L]
    metadata$UMAP3D_1 <- umap_3d[, 1L]
    metadata$UMAP3D_2 <- umap_3d[, 2L]
    metadata$UMAP3D_3 <- umap_3d[, 3L]
    metadata$cluster <- cluster
    write.csv(metadata, file.path(output_dir, "UMAP_metadata.csv"), row.names = FALSE)
    write.csv(metadata, file.path(output_dir, "UMAP_3D_plot_data.csv"), row.names = FALSE)

    theme <- theme_classic(base_size = 11) +
      theme(legend.key.height = grid::unit(0.35, "cm"),
            axis.text = element_text(colour = "black"))
    plot_dataset <- ggplot(metadata, aes(UMAP1, UMAP2, color = dataset)) +
      geom_point(size = 0.35, alpha = 0.7) + theme +
      labs(title = paste(disease, "StabMap + Harmony + spatial graph"), color = "Dataset")
    plot_omics <- ggplot(metadata, aes(UMAP1, UMAP2, color = omics)) +
      geom_point(size = 0.35, alpha = 0.7) + theme +
      labs(title = paste(disease, "by omics"), color = "Omics")
    plot_type <- ggplot(metadata, aes(UMAP1, UMAP2, color = broad_celltype)) +
      geom_point(size = 0.35, alpha = 0.7) + theme +
      labs(title = paste(disease, "broad cell types"), color = "Cell type")
    save_png(file.path(output_dir, "UMAP_by_dataset.png"), plot_dataset)
    save_png(file.path(output_dir, "UMAP_by_omics.png"), plot_omics)
    save_png(file.path(output_dir, "UMAP_by_celltype.png"), plot_type)
    save_umap_3d_plot(
      umap_3d, metadata$dataset, file.path(output_dir, "UMAP_3D_by_dataset.png"),
      paste(disease, "3D UMAP by dataset"), "Dataset"
    )
    save_umap_3d_plot(
      umap_3d, metadata$omics, file.path(output_dir, "UMAP_3D_by_omics.png"),
      paste(disease, "3D UMAP by omics"), "Omics"
    )
    save_umap_3d_plot(
      umap_3d, metadata$broad_celltype,
      file.path(output_dir, "UMAP_3D_by_celltype.png"),
      paste(disease, "3D UMAP by cell type"), "Cell type"
    )

    for (name in spatial_datasets) {
      coords <- sources[[name]]$coords
      ids <- intersect(rownames(coords), rownames(metadata))
      spatial_data <- cbind(metadata[ids, , drop = FALSE], coords[ids, , drop = FALSE])
      spatial_plot <- ggplot(spatial_data, aes(x, y, color = cluster)) +
        geom_point(size = 1) + scale_y_reverse() + coord_fixed() + theme +
        labs(title = paste(disease, name, "spatial clusters"), color = "Cluster")
      save_png(file.path(output_dir, paste0("Spatial_", name, "_clusters.png")),
               spatial_plot, 8, 7)
    }

    saveRDS(list(raw_StabMap = stab_raw, Harmony = stab_harmony,
                 final = final_embedding, umap = umap, umap_3d = umap_3d,
                 metadata = metadata,
                 metrics = metrics, references = references,
                 selected_genes = selected_genes), complete_file)
    writeLines(c(
      paste0("Disease: ", disease),
      "Core embedding: StabMap 0.1.8",
      paste0("References: ", paste(references, collapse = ", ")),
      "Batch correction: Harmony by dataset (theta=2, lambda=1)",
      "Spatial refinement: 10% propagation over each sample's true-coordinate 6-NN graph",
      "3D visualization: three-component UMAP rendered with scatterplot3d (T5.R method)",
      paste0("Datasets: ", paste(names(sources), collapse = ", "))
    ), file.path(output_dir, "README_method.txt"))
    data.frame(disease = disease, status = "success", datasets = length(sources),
               omics = length(unique(metadata$omics)), cells = nrow(metadata),
               elapsed_minutes = round(as.numeric(difftime(Sys.time(), started,
                                                            units = "mins")), 2),
               message = "", updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
               stringsAsFactors = FALSE)
  }, error = function(e) {
    log_message("FAILED", disease, conditionMessage(e))
    data.frame(disease = disease, status = "failed", datasets = nrow(rows),
               omics = length(unique(rows$omics)), cells = NA_integer_,
               elapsed_minutes = round(as.numeric(difftime(Sys.time(), started,
                                                            units = "mins")), 2),
               message = conditionMessage(e),
               updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
               stringsAsFactors = FALSE)
  })
  record_integration(result)
  rm(result); invisible(gc(full = TRUE))
}

cat("Disease-level integration pass complete.\n")
