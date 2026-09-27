gc()

# ============================================================
# Batch Module 3: Cell State Association
#
# Input:
#   E:/数据库/GSE*.rds
#
# Output:
#   E:/数据库/state_association_batch
#
# Main plots:
#   01_state_score_boxplot_by_celltype
#   02_state_score_heatmap
#   03_state_pseudotime_trend
#   04_state_expression_scatter_example
#   05_state_group_gene_boxplot_example
#
# Main web tables:
#   all_state_score_boxplot_raw.csv
#   all_state_score_boxplot_summary.csv
#   all_state_pseudotime_trend.csv
#   all_state_expression_scatter.csv
#   all_state_expression_correlation.csv
#   all_state_group_gene_boxplot_raw.csv
#   all_state_group_gene_boxplot_summary.csv
# ============================================================

# ============================================================
# 0. Load packages
# ============================================================
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

set.seed(123)

# ============================================================
# 1. Path settings
# ============================================================

rds_dir <- "E:/数据库"

out_root <- "E:/数据库/state_association_batch"

pseudotime_root <- "E:/数据库/pseudotime_web_data"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

max_dataset_n <- 10

skip_existing <- FALSE

save_example_plots <- TRUE

# ============================================================
# 2. Parameters
# ============================================================

assay_use <- "RNA"
slot_use <- "data"

preferred_celltype_col <- "celltype_label.main"

candidate_celltype_cols <- unique(c(
  preferred_celltype_col,
  "celltype_label.main",
  "celltype2",
  "celltype",
  "cell_type",
  "celltype_label.fine",
  "SingleR_annotation",
  "SingleR.labels",
  "annotation",
  "CellType",
  "celltype_label"
))

use_only_pseudotime_cells <- TRUE

pseudotime_bin_number <- 50

selected_cor_genes <- unique(c(
  "MKI67", "TOP2A", "PCNA",
  "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "GZMB", "PRF1",
  "LYZ", "S100A8", "S100A9", "IL1B", "TNF",
  "MS4A1", "CD79A", "CD79B",
  "SOX2", "NANOG", "POU5F1", "EPCAM", "LGR5",
  "VEGFA", "HIF1A", "CXCL8", "CCL2"
))

max_cor_gene_number <- 30

# ============================================================
# 3. Helper functions
# ============================================================

write_table <- function(df, file) {
  data.table::fwrite(df, file, na = "")
}

get_dataset_id <- function(file_path) {
  
  x <- basename(file_path)
  x <- gsub("\\.rds$", "", x, ignore.case = TRUE)
  
  gse <- sub(".*?(GSE[0-9]+).*", "\\1", x)
  
  if (grepl("^GSE[0-9]+$", gse)) {
    return(gse)
  }
  
  return(x)
}

get_assay_matrix <- function(obj, assay, layer_name = "data") {
  
  mat <- tryCatch(
    {
      GetAssayData(
        object = obj,
        assay = assay,
        layer = layer_name
      )
    },
    error = function(e1) {
      GetAssayData(
        object = obj,
        assay = assay,
        slot = layer_name
      )
    }
  )
  
  return(mat)
}

resolve_celltype_col <- function(meta_df) {
  
  available_cols <- candidate_celltype_cols[
    candidate_celltype_cols %in% colnames(meta_df)
  ]
  
  if (length(available_cols) == 0) {
    stop(
      "No cell type column found. Available columns: ",
      paste(colnames(meta_df), collapse = ", ")
    )
  }
  
  return(available_cols[1])
}

match_genes <- function(genes, features) {
  
  gene_map <- setNames(
    features,
    toupper(features)
  )
  
  genes_use <- unname(
    gene_map[toupper(unique(genes))]
  )
  
  genes_use <- unique(
    genes_use[!is.na(genes_use)]
  )
  
  return(genes_use)
}

score_one_gene_set <- function(expr_mat, genes) {
  
  genes_use <- match_genes(
    genes = genes,
    features = rownames(expr_mat)
  )
  
  if (length(genes_use) == 0) {
    return(rep(NA_real_, ncol(expr_mat)))
  }
  
  score <- Matrix::colMeans(
    expr_mat[genes_use, , drop = FALSE]
  )
  
  return(as.numeric(score))
}

scale_safe <- function(x) {
  
  if (length(x) == 0) {
    return(numeric())
  }
  
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  
  s <- sd(x, na.rm = TRUE)
  
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  
  return(as.numeric(scale(x)))
}

find_rds_files <- function(rds_dir) {
  
  files <- list.files(
    path = rds_dir,
    pattern = "^GSE[0-9]+\\.rds$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  
  files <- sort(files)
  
  return(files)
}

find_pseudotime_file <- function(dataset_id, pseudotime_root) {
  
  candidate_files <- c(
    file.path(pseudotime_root, dataset_id, paste0(dataset_id, "_pseudotime_trajectory_for_web.csv")),
    file.path(pseudotime_root, dataset_id, paste0(dataset_id, "_pseudotime_result.csv")),
    file.path(pseudotime_root, paste0(dataset_id, "_pseudotime_trajectory_for_web.csv")),
    file.path(pseudotime_root, paste0(dataset_id, "_pseudotime_result.csv"))
  )
  
  candidate_files <- candidate_files[file.exists(candidate_files)]
  
  if (length(candidate_files) > 0) {
    return(candidate_files[1])
  }
  
  all_file <- file.path(
    pseudotime_root,
    "all_pseudotime_trajectory_for_web.csv"
  )
  
  if (file.exists(all_file)) {
    return(all_file)
  }
  
  all_file2 <- file.path(
    pseudotime_root,
    "all_pseudotime_result.csv"
  )
  
  if (file.exists(all_file2)) {
    return(all_file2)
  }
  
  return(NA_character_)
}

standardize_pseudotime_table <- function(pseudotime_df, dataset_id) {
  
  if ("dataset_id" %in% colnames(pseudotime_df)) {
    pseudotime_df <- pseudotime_df %>%
      dplyr::filter(dataset_id == dataset_id)
  }
  
  if (!"cell_id" %in% colnames(pseudotime_df)) {
    stop("Pseudotime table does not contain cell_id.")
  }
  
  if (!"pseudotime" %in% colnames(pseudotime_df)) {
    
    if ("Pseudotime" %in% colnames(pseudotime_df)) {
      pseudotime_df$pseudotime <- pseudotime_df$Pseudotime
    } else {
      stop("Pseudotime table does not contain pseudotime.")
    }
  }
  
  if (!"state" %in% colnames(pseudotime_df)) {
    
    if ("State" %in% colnames(pseudotime_df)) {
      pseudotime_df$state <- pseudotime_df$State
    } else if ("trajectory_state" %in% colnames(pseudotime_df)) {
      pseudotime_df$state <- pseudotime_df$trajectory_state
    } else {
      pseudotime_df$state <- NA
    }
  }
  
  if (!"branch" %in% colnames(pseudotime_df)) {
    pseudotime_df$branch <- paste0("Branch_", pseudotime_df$state)
  }
  
  if (!"trajectory_cluster" %in% colnames(pseudotime_df)) {
    pseudotime_df$trajectory_cluster <- NA
  }
  
  pseudotime_use <- pseudotime_df %>%
    dplyr::filter(
      !is.na(cell_id),
      cell_id != ""
    ) %>%
    dplyr::distinct(cell_id, .keep_all = TRUE) %>%
    dplyr::select(
      cell_id,
      pseudotime,
      trajectory_state = state,
      branch,
      trajectory_cluster
    )
  
  return(pseudotime_use)
}

cor_one_pair <- function(df) {
  
  df <- df %>%
    dplyr::filter(
      !is.na(state_score),
      !is.na(log_expression)
    )
  
  if (nrow(df) < 5) {
    return(data.frame(
      n_cell = nrow(df),
      spearman_cor = NA_real_,
      spearman_p = NA_real_,
      pearson_cor = NA_real_,
      pearson_p = NA_real_
    ))
  }
  
  spearman_test <- suppressWarnings(
    cor.test(
      df$state_score,
      df$log_expression,
      method = "spearman"
    )
  )
  
  pearson_test <- suppressWarnings(
    cor.test(
      df$state_score,
      df$log_expression,
      method = "pearson"
    )
  )
  
  data.frame(
    n_cell = nrow(df),
    spearman_cor = as.numeric(spearman_test$estimate),
    spearman_p = spearman_test$p.value,
    pearson_cor = as.numeric(pearson_test$estimate),
    pearson_p = pearson_test$p.value
  )
}

theme_example <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

save_plot <- function(p, plot_dir, filename, width = 8, height = 6) {
  
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = p,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}

# ============================================================
# 4. State gene sets
# ============================================================

cc_genes_obj <- tryCatch(
  {
    get("cc.genes.updated.2019", envir = asNamespace("Seurat"))
  },
  error = function(e) {
    get("cc.genes", envir = asNamespace("Seurat"))
  }
)

s_genes <- cc_genes_obj$s.genes
g2m_genes <- cc_genes_obj$g2m.genes

state_gene_sets <- list(
  
  CellCycle_S = s_genes,
  
  CellCycle_G2M = g2m_genes,
  
  Proliferation = c(
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4",
    "MCM5", "MCM6", "MCM7", "UBE2C", "BIRC5", "CCNB1", "CCNB2"
  ),
  
  Stemness = c(
    "SOX2", "NANOG", "POU5F1", "KLF4", "MYC",
    "LIN28A", "PROM1", "ALDH1A1", "LGR5", "EPCAM"
  ),
  
  Immune_Activity = c(
    "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "GZMB",
    "PRF1", "IFNG", "CXCL9", "CXCL10", "HLA-DRA",
    "HLA-DRB1", "CD74"
  ),
  
  Inflammation = c(
    "IL1B", "IL6", "TNF", "NFKBIA", "CXCL8",
    "CCL2", "CCL3", "CCL4", "JUN", "FOS"
  ),
  
  Hypoxia = c(
    "HIF1A", "VEGFA", "LDHA", "ENO1", "PGK1",
    "SLC2A1", "CA9", "BNIP3", "NDRG1"
  )
)

# ============================================================
# 5. Single dataset function
# ============================================================

process_one_dataset <- function(rds_file) {
  
  dataset_id_value <- get_dataset_id(rds_file)
  
  message("================================================")
  message("Start dataset: ", dataset_id_value)
  message("Input RDS: ", rds_file)
  
  outdir <- file.path(out_root, dataset_id_value)
  plot_dir <- file.path(outdir, "example_plots")
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  da_v4 <- readRDS(rds_file)
  
  if (!assay_use %in% Assays(da_v4)) {
    stop("Assay not found: ", assay_use)
  }
  
  DefaultAssay(da_v4) <- assay_use
  
  celltype_col <- resolve_celltype_col(da_v4@meta.data)
  
  message("Cell type column used: ", celltype_col)
  
  expr_data <- get_assay_matrix(
    da_v4,
    assay = assay_use,
    layer_name = slot_use
  )
  
  if ((nrow(expr_data) == 0 || ncol(expr_data) == 0) &&
      slot_use == "data") {
    
    message("Data slot is empty. Running NormalizeData.")
    
    da_v4 <- NormalizeData(
      da_v4,
      assay = assay_use,
      verbose = FALSE
    )
    
    expr_data <- get_assay_matrix(
      da_v4,
      assay = assay_use,
      layer_name = slot_use
    )
  }
  
  expr_data <- as(expr_data, "sparseMatrix")
  
  meta_df <- da_v4@meta.data
  meta_df$cell_id <- rownames(meta_df)
  meta_df$cell_type <- as.character(meta_df[[celltype_col]])
  
  if ("orig.ident" %in% colnames(meta_df)) {
    meta_df$sample_name <- as.character(meta_df$orig.ident)
  } else {
    meta_df$sample_name <- dataset_id_value
  }
  
  if ("seurat_clusters" %in% colnames(meta_df)) {
    meta_df$cluster_id <- as.character(meta_df$seurat_clusters)
  } else {
    meta_df$cluster_id <- NA
  }
  
  rownames(meta_df) <- meta_df$cell_id
  
  # ============================================================
  # Merge pseudotime if available
  # ============================================================
  
  pseudotime_file <- find_pseudotime_file(
    dataset_id = dataset_id_value,
    pseudotime_root = pseudotime_root
  )
  
  if (!is.na(pseudotime_file)) {
    
    message("Reading pseudotime table: ", pseudotime_file)
    
    pseudotime_df <- data.table::fread(
      pseudotime_file,
      data.table = FALSE,
      check.names = FALSE
    )
    
    pseudotime_use <- standardize_pseudotime_table(
      pseudotime_df = pseudotime_df,
      dataset_id = dataset_id_value
    )
    
    if (use_only_pseudotime_cells) {
      
      matched_cells <- pseudotime_use$cell_id[
        pseudotime_use$cell_id %in% colnames(expr_data)
      ]
      
      message("Matched pseudotime cells: ", length(matched_cells))
      
      if (length(matched_cells) >= 50) {
        
        expr_data <- expr_data[, matched_cells, drop = FALSE]
        meta_df <- meta_df[matched_cells, , drop = FALSE]
        
        pseudotime_use <- pseudotime_use[
          match(matched_cells, pseudotime_use$cell_id),
          ,
          drop = FALSE
        ]
      }
    }
    
    meta_df$cell_id <- rownames(meta_df)
    
    meta_df <- meta_df %>%
      dplyr::left_join(
        pseudotime_use,
        by = "cell_id"
      )
    
    rownames(meta_df) <- meta_df$cell_id
    
  } else {
    
    message("No pseudotime table found. Continuing without pseudotime.")
    
    meta_df$pseudotime <- NA
    meta_df$trajectory_state <- NA
    meta_df$branch <- NA
    meta_df$trajectory_cluster <- NA
  }
  
  meta_df <- meta_df[colnames(expr_data), , drop = FALSE]
  
  stopifnot(
    identical(
      rownames(meta_df),
      colnames(expr_data)
    )
  )
  
  # ============================================================
  # Compute state scores
  # ============================================================
  
  state_score_mat <- data.frame(
    cell_id = colnames(expr_data),
    stringsAsFactors = FALSE
  )
  
  for (score_name in names(state_gene_sets)) {
    
    message("Computing state score: ", score_name)
    
    state_score_mat[[score_name]] <- score_one_gene_set(
      expr_mat = expr_data,
      genes = state_gene_sets[[score_name]]
    )
  }
  
  if (all(c("CellCycle_S", "CellCycle_G2M") %in% colnames(state_score_mat))) {
    state_score_mat$CellCycle_Diff_G2M_minus_S <-
      state_score_mat$CellCycle_G2M - state_score_mat$CellCycle_S
  }
  
  score_cols <- setdiff(
    colnames(state_score_mat),
    "cell_id"
  )
  
  state_score_cell <- meta_df %>%
    dplyr::left_join(
      state_score_mat,
      by = "cell_id"
    ) %>%
    dplyr::mutate(
      dataset_id = dataset_id_value
    ) %>%
    dplyr::select(
      dataset_id,
      cell_id,
      cell_type,
      sample_name,
      cluster_id,
      pseudotime,
      trajectory_state,
      branch,
      trajectory_cluster,
      dplyr::all_of(score_cols)
    )
  
  write_table(
    state_score_cell,
    file.path(outdir, paste0(dataset_id_value, "_state_score_cell.csv"))
  )
  
  state_score_long <- state_score_cell %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(score_cols),
      names_to = "state_name",
      values_to = "state_score"
    ) %>%
    dplyr::filter(!is.na(state_score)) %>%
    dplyr::select(
      dataset_id,
      cell_id,
      cell_type,
      sample_name,
      cluster_id,
      pseudotime,
      trajectory_state,
      branch,
      trajectory_cluster,
      state_name,
      state_score
    )
  
  write_table(
    state_score_long,
    file.path(outdir, paste0(dataset_id_value, "_state_score_long.csv"))
  )
  
  # ============================================================
  # Boxplot raw / summary
  # ============================================================
  
  state_score_boxplot_raw <- state_score_long %>%
    dplyr::filter(
      !is.na(cell_type),
      cell_type != "",
      !is.na(state_score)
    ) %>%
    dplyr::mutate(
      id = dplyr::row_number()
    ) %>%
    dplyr::select(
      id,
      dataset_id,
      cell_id,
      cell_type,
      sample_name,
      cluster_id,
      state_name,
      state_score
    )
  
  write_table(
    state_score_boxplot_raw,
    file.path(outdir, paste0(dataset_id_value, "_state_score_boxplot_raw.csv"))
  )
  
  state_score_boxplot_summary <- state_score_boxplot_raw %>%
    dplyr::group_by(dataset_id, state_name, cell_type) %>%
    dplyr::summarise(
      cell_count = dplyr::n(),
      min_score = min(state_score, na.rm = TRUE),
      q1_score = quantile(state_score, 0.25, na.rm = TRUE),
      median_score = median(state_score, na.rm = TRUE),
      q3_score = quantile(state_score, 0.75, na.rm = TRUE),
      max_score = max(state_score, na.rm = TRUE),
      mean_score = mean(state_score, na.rm = TRUE),
      sd_score = sd(state_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(state_name) %>%
    dplyr::mutate(
      z_score = scale_safe(mean_score)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      id = dplyr::row_number()
    ) %>%
    dplyr::select(
      id,
      dataset_id,
      state_name,
      cell_type,
      cell_count,
      min_score,
      q1_score,
      median_score,
      q3_score,
      max_score,
      mean_score,
      sd_score,
      z_score
    )
  
  write_table(
    state_score_boxplot_summary,
    file.path(outdir, paste0(dataset_id_value, "_state_score_boxplot_summary.csv"))
  )
  
  # ============================================================
  # Pseudotime trend
  # ============================================================
  
  if (sum(!is.na(state_score_cell$pseudotime)) > 0) {
    
    ptime_min <- min(state_score_cell$pseudotime, na.rm = TRUE)
    ptime_max <- max(state_score_cell$pseudotime, na.rm = TRUE)
    
    if (is.finite(ptime_min) &&
        is.finite(ptime_max) &&
        ptime_max > ptime_min) {
      
      ptime_breaks <- seq(
        ptime_min,
        ptime_max,
        length.out = pseudotime_bin_number + 1
      )
      
      state_score_long$pseudotime_bin <- cut(
        state_score_long$pseudotime,
        breaks = ptime_breaks,
        labels = FALSE,
        include.lowest = TRUE
      )
      
      state_pseudotime_trend <- state_score_long %>%
        dplyr::filter(
          !is.na(pseudotime),
          !is.na(pseudotime_bin)
        ) %>%
        dplyr::group_by(
          dataset_id,
          state_name,
          pseudotime_bin
        ) %>%
        dplyr::summarise(
          pseudotime_min = min(pseudotime, na.rm = TRUE),
          pseudotime_max = max(pseudotime, na.rm = TRUE),
          pseudotime_mid = mean(pseudotime, na.rm = TRUE),
          mean_score = mean(state_score, na.rm = TRUE),
          median_score = median(state_score, na.rm = TRUE),
          sd_score = sd(state_score, na.rm = TRUE),
          cell_count = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::group_by(state_name) %>%
        dplyr::mutate(
          scaled_mean_score = scale_safe(mean_score)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          id = dplyr::row_number()
        ) %>%
        dplyr::select(
          id,
          dataset_id,
          state_name,
          pseudotime_bin,
          pseudotime_min,
          pseudotime_max,
          pseudotime_mid,
          mean_score,
          median_score,
          sd_score,
          scaled_mean_score,
          cell_count
        )
    } else {
      state_pseudotime_trend <- data.frame()
    }
    
  } else {
    state_pseudotime_trend <- data.frame()
  }
  
  write_table(
    state_pseudotime_trend,
    file.path(outdir, paste0(dataset_id_value, "_state_pseudotime_trend.csv"))
  )
  
  # ============================================================
  # State-expression scatter / correlation
  # ============================================================
  
  selected_cor_genes_use <- match_genes(
    genes = selected_cor_genes,
    features = rownames(expr_data)
  )
  
  selected_cor_genes_use <- selected_cor_genes_use[
    1:min(max_cor_gene_number, length(selected_cor_genes_use))
  ]
  
  if (length(selected_cor_genes_use) == 0) {
    
    state_expression_scatter <- data.frame()
    state_expression_correlation <- data.frame()
    state_group_gene_boxplot_raw <- data.frame()
    state_group_gene_boxplot_summary <- data.frame()
    
  } else {
    
    expr_selected <- expr_data[
      selected_cor_genes_use,
      state_score_cell$cell_id,
      drop = FALSE
    ]
    
    expr_selected_dense <- as.matrix(expr_selected)
    
    expr_selected_long <- as.data.frame(
      as.table(expr_selected_dense)
    )
    
    colnames(expr_selected_long) <- c(
      "gene_symbol",
      "cell_id",
      "expression"
    )
    
    expr_selected_long$gene_symbol <- as.character(expr_selected_long$gene_symbol)
    expr_selected_long$cell_id <- as.character(expr_selected_long$cell_id)
    expr_selected_long$expression <- as.numeric(expr_selected_long$expression)
    expr_selected_long$log_expression <- log1p(expr_selected_long$expression)
    
    state_score_for_join <- state_score_long %>%
      dplyr::select(
        dataset_id,
        cell_id,
        cell_type,
        sample_name,
        cluster_id,
        pseudotime,
        trajectory_state,
        state_name,
        state_score
      )
    
    state_expression_scatter <- expr_selected_long %>%
      dplyr::left_join(
        state_score_for_join,
        by = "cell_id"
      ) %>%
      dplyr::filter(!is.na(state_score)) %>%
      dplyr::mutate(
        id = dplyr::row_number()
      ) %>%
      dplyr::select(
        id,
        dataset_id,
        cell_id,
        cell_type,
        sample_name,
        cluster_id,
        pseudotime,
        trajectory_state,
        state_name,
        state_score,
        gene_symbol,
        expression,
        log_expression
      )
    
    write_table(
      state_expression_scatter,
      file.path(outdir, paste0(dataset_id_value, "_state_expression_scatter.csv"))
    )
    
    state_expression_correlation <- state_expression_scatter %>%
      dplyr::group_by(dataset_id, state_name, gene_symbol) %>%
      dplyr::group_modify(~ cor_one_pair(.x)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(dataset_id, state_name) %>%
      dplyr::mutate(
        spearman_p_adj = p.adjust(spearman_p, method = "BH"),
        pearson_p_adj = p.adjust(pearson_p, method = "BH")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(
        state_name,
        spearman_p_adj,
        dplyr::desc(abs(spearman_cor))
      ) %>%
      dplyr::mutate(
        id = dplyr::row_number()
      ) %>%
      dplyr::select(
        id,
        dataset_id,
        state_name,
        gene_symbol,
        n_cell,
        spearman_cor,
        spearman_p,
        spearman_p_adj,
        pearson_cor,
        pearson_p,
        pearson_p_adj
      )
    
    write_table(
      state_expression_correlation,
      file.path(outdir, paste0(dataset_id_value, "_state_expression_correlation.csv"))
    )
    
    state_score_group <- state_score_long %>%
      dplyr::group_by(state_name) %>%
      dplyr::mutate(
        q_low = quantile(state_score, 1 / 3, na.rm = TRUE),
        q_high = quantile(state_score, 2 / 3, na.rm = TRUE),
        state_score_group = dplyr::case_when(
          state_score <= q_low ~ "Low",
          state_score >= q_high ~ "High",
          TRUE ~ "Mid"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        dataset_id,
        cell_id,
        cell_type,
        sample_name,
        cluster_id,
        state_name,
        state_score,
        state_score_group
      )
    
    state_group_gene_boxplot_raw <- expr_selected_long %>%
      dplyr::left_join(
        state_score_group,
        by = "cell_id"
      ) %>%
      dplyr::filter(!is.na(state_score_group)) %>%
      dplyr::mutate(
        id = dplyr::row_number(),
        log_expression = log1p(expression)
      ) %>%
      dplyr::select(
        id,
        dataset_id,
        cell_id,
        cell_type,
        sample_name,
        cluster_id,
        state_name,
        state_score_group,
        state_score,
        gene_symbol,
        expression,
        log_expression
      )
    
    write_table(
      state_group_gene_boxplot_raw,
      file.path(outdir, paste0(dataset_id_value, "_state_group_gene_boxplot_raw.csv"))
    )
    
    state_group_gene_boxplot_summary <- state_group_gene_boxplot_raw %>%
      dplyr::group_by(
        dataset_id,
        state_name,
        state_score_group,
        gene_symbol
      ) %>%
      dplyr::summarise(
        cell_count = dplyr::n(),
        min_expression = min(log_expression, na.rm = TRUE),
        q1_expression = quantile(log_expression, 0.25, na.rm = TRUE),
        median_expression = median(log_expression, na.rm = TRUE),
        q3_expression = quantile(log_expression, 0.75, na.rm = TRUE),
        max_expression = max(log_expression, na.rm = TRUE),
        mean_expression = mean(log_expression, na.rm = TRUE),
        sd_expression = sd(log_expression, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        id = dplyr::row_number()
      )
    
    write_table(
      state_group_gene_boxplot_summary,
      file.path(outdir, paste0(dataset_id_value, "_state_group_gene_boxplot_summary.csv"))
    )
  }
  
  # ============================================================
  # Example plots
  # ============================================================
  
  if (save_example_plots) {
    
    p_state_box <- ggplot(
      state_score_boxplot_raw,
      aes(
        x = cell_type,
        y = state_score,
        fill = cell_type
      )
    ) +
      geom_boxplot(outlier.size = 0.2) +
      facet_wrap(~ state_name, scales = "free_y") +
      coord_flip() +
      theme_example +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 8)
      ) +
      labs(
        x = "Cell type",
        y = "State score",
        title = paste0(dataset_id_value, " cell state scores by cell type")
      )
    
    save_plot(
      p_state_box,
      plot_dir,
      paste0(dataset_id_value, "_01_state_score_boxplot_by_celltype"),
      width = 12,
      height = 8
    )
    
    heatmap_df <- state_score_boxplot_summary %>%
      dplyr::filter(is.finite(z_score))
    
    if (nrow(heatmap_df) > 0) {
      
      p_heatmap <- ggplot(
        heatmap_df,
        aes(
          x = cell_type,
          y = state_name,
          fill = z_score
        )
      ) +
        geom_tile(color = "white", linewidth = 0.2) +
        theme_example +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(
          x = "Cell type",
          y = "Cell state",
          fill = "Z-score",
          title = paste0(dataset_id_value, " state activity heatmap")
        )
      
      save_plot(
        p_heatmap,
        plot_dir,
        paste0(dataset_id_value, "_02_state_score_heatmap"),
        width = 9,
        height = 5.5
      )
    }
    
    if (nrow(state_pseudotime_trend) > 0) {
      
      p_trend <- ggplot(
        state_pseudotime_trend,
        aes(
          x = pseudotime_mid,
          y = scaled_mean_score,
          color = state_name,
          group = state_name
        )
      ) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 0.8) +
        facet_wrap(~ state_name, scales = "free_y", ncol = 2) +
        theme_example +
        theme(
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 8)
        ) +
        labs(
          x = "Pseudotime",
          y = "Scaled mean state score",
          title = paste0(dataset_id_value, " state scores along pseudotime")
        )
      
      save_plot(
        p_trend,
        plot_dir,
        paste0(dataset_id_value, "_03_state_pseudotime_trend"),
        width = 9,
        height = 8
      )
    }
    
    if (nrow(state_expression_scatter) > 0 &&
        nrow(state_expression_correlation) > 0) {
      
      plot_state <- state_expression_correlation$state_name[1]
      plot_gene <- state_expression_correlation$gene_symbol[1]
      
      scatter_df <- state_expression_scatter %>%
        dplyr::filter(
          state_name == plot_state,
          gene_symbol == plot_gene
        )
      
      if (nrow(scatter_df) > 0) {
        
        p_scatter <- ggplot(
          scatter_df,
          aes(
            x = state_score,
            y = log_expression,
            color = cell_type
          )
        ) +
          geom_point(size = 0.5, alpha = 0.5) +
          geom_smooth(method = "lm", se = FALSE, color = "black") +
          theme_example +
          labs(
            x = paste0(plot_state, " score"),
            y = paste0(plot_gene, " log1p expression"),
            color = "Cell type",
            title = paste0(dataset_id_value, " ", plot_gene, " vs ", plot_state)
          )
        
        save_plot(
          p_scatter,
          plot_dir,
          paste0(dataset_id_value, "_04_state_expression_scatter_example"),
          width = 7,
          height = 5.5
        )
      }
    }
    
    if (nrow(state_group_gene_boxplot_raw) > 0) {
      
      plot_state2 <- unique(state_group_gene_boxplot_raw$state_name)[1]
      
      genes_show <- unique(state_group_gene_boxplot_raw$gene_symbol)
      genes_show <- genes_show[1:min(6, length(genes_show))]
      
      group_box_df <- state_group_gene_boxplot_raw %>%
        dplyr::filter(
          state_name == plot_state2,
          gene_symbol %in% genes_show
        )
      
      p_group_box <- ggplot(
        group_box_df,
        aes(
          x = state_score_group,
          y = log_expression,
          fill = state_score_group
        )
      ) +
        geom_boxplot(outlier.size = 0.2) +
        facet_wrap(~ gene_symbol, scales = "free_y") +
        theme_example +
        theme(
          legend.position = "none",
          strip.text = element_text(face = "bold")
        ) +
        labs(
          x = paste0(plot_state2, " group"),
          y = "log1p expression",
          title = paste0(dataset_id_value, " gene expression across ", plot_state2, " groups")
        )
      
      save_plot(
        p_group_box,
        plot_dir,
        paste0(dataset_id_value, "_05_state_group_gene_boxplot_example"),
        width = 10,
        height = 6
      )
    }
  }
  
  message("Finished dataset: ", dataset_id_value)
  
  result <- list(
    state_score_cell = state_score_cell,
    state_score_long = state_score_long,
    state_score_boxplot_raw = state_score_boxplot_raw,
    state_score_boxplot_summary = state_score_boxplot_summary,
    state_pseudotime_trend = state_pseudotime_trend,
    state_expression_scatter = state_expression_scatter,
    state_expression_correlation = state_expression_correlation,
    state_group_gene_boxplot_raw = state_group_gene_boxplot_raw,
    state_group_gene_boxplot_summary = state_group_gene_boxplot_summary,
    cell_number = ncol(expr_data)
  )
  
  rm(
    da_v4,
    expr_data,
    meta_df,
    state_score_mat,
    state_score_cell,
    state_score_long
  )
  
  gc()
  
  return(result)
}

# ============================================================
# 6. Batch loop
# ============================================================

rds_files <- find_rds_files(rds_dir)

if (length(rds_files) == 0) {
  stop("No GSE*.rds files were found in: ", rds_dir)
}

rds_files <- rds_files[
  1:min(max_dataset_n, length(rds_files))
]

message("Datasets to process:")
print(basename(rds_files))

run_log <- data.frame(
  dataset_id = character(),
  status = character(),
  message = character(),
  input_file = character(),
  cell_number = numeric(),
  output_dir = character(),
  stringsAsFactors = FALSE
)

all_state_score_cell_list <- list()
all_state_score_long_list <- list()
all_state_score_boxplot_raw_list <- list()
all_state_score_boxplot_summary_list <- list()
all_state_pseudotime_trend_list <- list()
all_state_expression_scatter_list <- list()
all_state_expression_correlation_list <- list()
all_state_group_gene_boxplot_raw_list <- list()
all_state_group_gene_boxplot_summary_list <- list()

for (rds_file in rds_files) {
  
  dataset_id_value <- get_dataset_id(rds_file)
  outdir_i <- file.path(out_root, dataset_id_value)
  
  tryCatch({
    
    res <- process_one_dataset(rds_file)
    
    all_state_score_cell_list[[dataset_id_value]] <- res$state_score_cell
    all_state_score_long_list[[dataset_id_value]] <- res$state_score_long
    all_state_score_boxplot_raw_list[[dataset_id_value]] <- res$state_score_boxplot_raw
    all_state_score_boxplot_summary_list[[dataset_id_value]] <- res$state_score_boxplot_summary
    all_state_pseudotime_trend_list[[dataset_id_value]] <- res$state_pseudotime_trend
    all_state_expression_scatter_list[[dataset_id_value]] <- res$state_expression_scatter
    all_state_expression_correlation_list[[dataset_id_value]] <- res$state_expression_correlation
    all_state_group_gene_boxplot_raw_list[[dataset_id_value]] <- res$state_group_gene_boxplot_raw
    all_state_group_gene_boxplot_summary_list[[dataset_id_value]] <- res$state_group_gene_boxplot_summary
    
    run_log <- rbind(
      run_log,
      data.frame(
        dataset_id = dataset_id_value,
        status = "success",
        message = "Finished",
        input_file = rds_file,
        cell_number = res$cell_number,
        output_dir = outdir_i,
        stringsAsFactors = FALSE
      )
    )
    
  }, error = function(e) {
    
    message("Dataset failed: ", dataset_id_value)
    message("Error message: ", e$message)
    
    run_log <- rbind(
      run_log,
      data.frame(
        dataset_id = dataset_id_value,
        status = "failed",
        message = e$message,
        input_file = rds_file,
        cell_number = NA,
        output_dir = outdir_i,
        stringsAsFactors = FALSE
      )
    )
    
    gc()
  })
}

# ============================================================
# 7. Save merged tables
# ============================================================

save_merged <- function(list_obj, filename) {
  
  list_obj <- Filter(
    function(x) {
      is.data.frame(x) && nrow(x) > 0
    },
    list_obj
  )
  
  if (length(list_obj) == 0) {
    message("No data to merge for: ", filename)
    return(NULL)
  }
  
  merged_df <- data.table::rbindlist(
    list_obj,
    use.names = TRUE,
    fill = TRUE
  )
  
  merged_df$id <- seq_len(nrow(merged_df))
  
  write_table(
    merged_df,
    file.path(out_root, filename)
  )
  
  message("Merged table saved: ", filename, " rows: ", nrow(merged_df))
  
  return(merged_df)
}

save_merged(
  all_state_score_cell_list,
  "all_state_score_cell.csv"
)

save_merged(
  all_state_score_long_list,
  "all_state_score_long.csv"
)

save_merged(
  all_state_score_boxplot_raw_list,
  "all_state_score_boxplot_raw.csv"
)

save_merged(
  all_state_score_boxplot_summary_list,
  "all_state_score_boxplot_summary.csv"
)

save_merged(
  all_state_pseudotime_trend_list,
  "all_state_pseudotime_trend.csv"
)

save_merged(
  all_state_expression_scatter_list,
  "all_state_expression_scatter.csv"
)

save_merged(
  all_state_expression_correlation_list,
  "all_state_expression_correlation.csv"
)

save_merged(
  all_state_group_gene_boxplot_raw_list,
  "all_state_group_gene_boxplot_raw.csv"
)

save_merged(
  all_state_group_gene_boxplot_summary_list,
  "all_state_group_gene_boxplot_summary.csv"
)

# ============================================================
# 8. Save run log
# ============================================================

log_file <- file.path(
  out_root,
  "batch_state_association_run_log.csv"
)

write_table(
  run_log,
  log_file
)

message("================================================")
message("All state association jobs finished.")
message("Output root: ", out_root)
message("Run log: ", log_file)
