############################################################
## Batch CITE-seq-specific analysis after WNN + SingleR
##
## Outputs for every annotated Seurat RDS:
##   1. Annotation validation: UMAP, RNA/ADT marker DotPlots,
##      WNN-cluster × annotation tables and heatmap
##   2. ADT markers for each annotated cell type
##   3. Matched RNA–protein correlations within cell types
##   4. WNN RNA/ADT modality-weight summaries and plots
##   5. RNA differential genes plus GO/KEGG enrichment
##   6. Web-ready plotting tables and a figure-to-table manifest
##
## Designed for:
##   I:/CITE-seq 7.16/mouse 7.16/GSE*/
##   GSE*_WNN_SingleR_annotated.rds
##
## Plot titles, labels, legends and axes are in English.
############################################################


## =========================================================
## 0. First-time package installation
## Run only when a package is missing, then restart R/RStudio.
## =========================================================

# install.packages(
#   c(
#     "Seurat",
#     "SeuratObject",
#     "ggplot2",
#     "patchwork",
#     "Matrix"
#   ),
#   dependencies = TRUE
# )
#
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
#
# BiocManager::install(
#   c(
#     "clusterProfiler",
#     "AnnotationDbi",
#     "org.Mm.eg.db",
#     "org.Hs.eg.db"
#   ),
#   ask = FALSE,
#   update = FALSE
# )


## =========================================================
## 1. Load packages
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(clusterProfiler)
  library(AnnotationDbi)
})

set.seed(1234)

options(
  stringsAsFactors = FALSE,
  future.globals.maxSize = 20 * 1024^3
)

## Keep processing sequential to reduce peak memory.
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

cat(
  "Seurat version: ",
  as.character(packageVersion("Seurat")),
  "\n",
  sep = ""
)

cat(
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject")),
  "\n",
  sep = ""
)


## =========================================================
## 2. User settings
## Only this section normally needs editing.
## =========================================================

## Root folder containing GSE folders.
root_dir <- "I:\\CITE-seq 7.19\\mouse 7.19\\Seurat_RDS1"

## Annotated Seurat RDS filename pattern.
input_rds_pattern <- "_WNN_SingleR_annotated\\.rds$"

## Results are saved inside each GSE folder:
## GSEXXXXX/CITE_specific_analysis_GO_KEGG_webdata/
output_folder_name <- "CITE_specific_analysis_GO_KEGG_webdata"

## Species controls the built-in RNA marker panel and
## the RNA–protein matching dictionary.
species <- "mouse"
## Allowed values: "mouse" or "human"

## Resume settings.
overwrite_completed <- FALSE

## First test one RDS, then set Inf for the full batch.
test_first_n <- Inf

## Save an updated object containing:
## celltype_for_analysis and any derived ADT.weight column.
## FALSE saves disk space and avoids another large RDS.
save_analysis_ready_rds <- FALSE

## Preferred annotation columns, in priority order.
## The first usable column is selected. Missing values are then
## filled from lower-priority columns.
preferred_celltype_columns <- c(
  "celltype_final",
  "celltype_label.main_pruned",
  "celltype_label.main",
  "celltype_label.fine_pruned",
  "celltype_label.fine"
)

## Preferred cluster columns.
preferred_cluster_columns <- c(
  "wnn_clusters",
  "seurat_clusters",
  "rna_clusters"
)

## Preferred 2D reduction.
preferred_reductions <- c(
  "wnn.umap",
  "umap",
  "rna.umap"
)

## Minimum group sizes.
min_cells_for_marker_group <- 10L
min_cells_for_correlation <- 20L

## Memory-saving downsampling.
max_cells_per_group_for_markers <- 3000L
max_cells_per_group_for_correlation <- 3000L
max_cells_per_group_for_weight_plot <- 2000L
max_cells_for_scatter <- 8000L

## ADT marker thresholds.
adt_marker_logfc_threshold <- 0.10
adt_marker_min_pct <- 0.10
adt_top_markers_per_celltype <- 10L
adt_top_features_for_dotplot <- 40L

## Optional RNA marker discovery.
## Known RNA marker DotPlots are always produced when markers exist.
## FindAllMarkers on RNA is much heavier, so it is disabled by default.
run_rna_findallmarkers <- FALSE
rna_marker_logfc_threshold <- 0.25
rna_marker_min_pct <- 0.10
rna_top_markers_per_celltype <- 10L
rna_top_features_for_dotplot <- 40L

## RNA DEG and enrichment settings.
## Enrichment requires positive RNA markers for each cell type
## versus all remaining cells.
run_rna_deg_enrichment <- TRUE
rna_deg_padj_threshold <- 0.05
rna_deg_min_log2fc_for_enrichment <- 0.25
rna_deg_min_pct_for_enrichment <- 0.10
min_genes_for_enrichment <- 10L
max_genes_per_celltype_for_enrichment <- 1000L

go_ontologies <- c("BP", "CC", "MF")
enrichment_pvalue_cutoff <- 0.05
enrichment_qvalue_cutoff <- 0.20
enrichment_top_terms_per_database <- 10L

## Export one or more compact CSV tables that reproduce every plot.
export_web_plot_data <- TRUE
web_plot_data_folder_name <- "06_web_plot_data"

## RNA–protein correlation settings.
include_all_cells_correlation <- TRUE
max_correlation_pairs_in_heatmap <- 30L
number_of_scatter_pairs <- 6L

## Plot format.
plot_dpi <- 300


## =========================================================
## 3. Basic checks
## =========================================================

if (!dir.exists(root_dir)) {
  stop(
    "Root folder does not exist:\n",
    root_dir
  )
}

species <- tolower(species)

if (!species %in% c("mouse", "human")) {
  stop("species must be either 'mouse' or 'human'.")
}

orgdb_package <- if (species == "mouse") {
  "org.Mm.eg.db"
} else {
  "org.Hs.eg.db"
}

kegg_organism <- if (species == "mouse") {
  "mmu"
} else {
  "hsa"
}

if (!requireNamespace(orgdb_package, quietly = TRUE)) {
  stop(
    "Missing annotation package: ",
    orgdb_package,
    "\nInstall it with BiocManager::install('",
    orgdb_package,
    "')."
  )
}

org_db <- get(
  orgdb_package,
  envir = asNamespace(orgdb_package)
)


## =========================================================
## 4. General helper functions
## =========================================================

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}


safe_filename <- function(x) {
  
  x <- as.character(x)
  x <- gsub("[\\\\/:*?\"<>|]+", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  
  ifelse(
    nchar(x) == 0,
    "unnamed",
    x
  )
}


extract_gse_id <- function(path) {
  
  hit <- regmatches(
    path,
    regexpr(
      "GSE[0-9]+",
      path,
      ignore.case = TRUE
    )
  )
  
  if (length(hit) == 0 || is.na(hit) || hit == "") {
    return(
      tools::file_path_sans_ext(
        basename(path)
      )
    )
  }
  
  toupper(hit)
}


write_csv_safe <- function(x, file, row.names = FALSE) {
  
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(
    x,
    file = file,
    row.names = row.names,
    na = ""
  )
}


save_plot_safe <- function(
    plot_object,
    filename,
    width,
    height
) {
  
  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    dpi = plot_dpi,
    limitsize = FALSE
  )
}


write_text_safe <- function(text, file) {
  
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  writeLines(
    text = as.character(text),
    con = file,
    useBytes = TRUE
  )
}



normalize_path_for_table <- function(x) {
  gsub("\\\\", "/", as.character(x))
}


relative_path_safe <- function(path, root) {
  
  path_normalized <- normalize_path_for_table(
    normalizePath(
      path,
      winslash = "/",
      mustWork = FALSE
    )
  )
  
  root_normalized <- normalize_path_for_table(
    normalizePath(
      root,
      winslash = "/",
      mustWork = FALSE
    )
  )
  
  root_prefix <- paste0(
    sub("/+$", "", root_normalized),
    "/"
  )
  
  if (startsWith(path_normalized, root_prefix)) {
    substring(
      path_normalized,
      nchar(root_prefix) + 1L
    )
  } else {
    path_normalized
  }
}


get_assay_names_safe <- function(object) {
  as.character(
    base::names(object@assays)
  )
}


get_reduction_names_safe <- function(object) {
  as.character(
    base::names(object@reductions)
  )
}


get_graph_names_safe <- function(object) {
  as.character(
    base::names(object@graphs)
  )
}


get_neighbor_names_safe <- function(object) {
  as.character(
    base::names(object@neighbors)
  )
}


has_name_safe <- function(
    available_names,
    required_name
) {
  any(
    as.character(available_names) ==
      as.character(required_name)
  )
}


has_all_names_safe <- function(
    available_names,
    required_names
) {
  
  required_names <- as.character(
    required_names
  )
  
  all(
    vapply(
      required_names,
      function(current_name) {
        has_name_safe(
          available_names,
          current_name
        )
      },
      logical(1)
    )
  )
}


as_numeric_safe <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


factor_order_index <- function(x) {
  
  values <- as.character(x)
  
  levels_use <- if (is.factor(x)) {
    levels(x)
  } else {
    unique(values)
  }
  
  match(
    values,
    levels_use
  )
}


rbind_fill_data_frames <- function(data_list) {
  
  data_list <- Filter(
    function(x) {
      is.data.frame(x) &&
        nrow(x) > 0
    },
    data_list
  )
  
  if (length(data_list) == 0) {
    return(data.frame())
  }
  
  all_columns <- unique(
    unlist(
      lapply(
        data_list,
        colnames
      ),
      use.names = FALSE
    )
  )
  
  aligned <- lapply(
    data_list,
    function(current_data) {
      
      missing_columns <- setdiff(
        all_columns,
        colnames(current_data)
      )
      
      if (length(missing_columns) > 0) {
        for (current_column in missing_columns) {
          current_data[[current_column]] <- NA
        }
      }
      
      current_data <- current_data[
        ,
        all_columns,
        drop = FALSE
      ]
      
      for (current_column in colnames(current_data)) {
        if (is.factor(current_data[[current_column]])) {
          current_data[[current_column]] <-
            as.character(
              current_data[[current_column]]
            )
        }
      }
      
      current_data
    }
  )
  
  output <- do.call(
    rbind,
    aligned
  )
  
  rownames(output) <- NULL
  
  output
}


standardize_dotplot_data <- function(
    plot_object,
    dataset_id,
    plot_id,
    plot_file,
    modality
) {
  
  plot_data <- tryCatch(
    {
      as.data.frame(
        plot_object$data
      )
    },
    error = function(e) {
      data.frame()
    }
  )
  
  if (nrow(plot_data) == 0) {
    return(data.frame())
  }
  
  feature_column <- intersect(
    c(
      "features.plot",
      "feature",
      "features"
    ),
    colnames(plot_data)
  )
  
  group_column <- intersect(
    c(
      "id",
      "celltype",
      "ident"
    ),
    colnames(plot_data)
  )
  
  if (
    length(feature_column) == 0 ||
    length(group_column) == 0
  ) {
    return(data.frame())
  }
  
  get_column_or_na <- function(
    column_name,
    default = NA_real_
  ) {
    
    if (column_name %in% colnames(plot_data)) {
      return(plot_data[[column_name]])
    }
    
    rep(
      default,
      nrow(plot_data)
    )
  }
  
  output <- data.frame(
    dataset_id = dataset_id,
    plot_id = plot_id,
    plot_file = normalize_path_for_table(
      plot_file
    ),
    modality = modality,
    cell_type = as.character(plot_data[[group_column[1]]]),
    cell_type_order = factor_order_index(
      plot_data[[group_column[1] ]]
    ),
    feature = as.character(
      plot_data[[feature_column[1]] ]
    ),
    feature_order = factor_order_index(
      plot_data[[feature_column[1]]]
    ),
    avg_expression = as_numeric_safe(
      get_column_or_na("avg.exp")
    ),
    avg_expression_scaled = as_numeric_safe(
      get_column_or_na(
        "avg.exp.scaled"
      )
    ),
    percent_expressing = as_numeric_safe(
      get_column_or_na("pct.exp")
    ),
    stringsAsFactors = FALSE
  )
  
  output
}


select_first_existing_metadata_column <- function(
    object,
    candidates
) {
  
  candidates <- as.character(
    candidates
  )
  
  hit <- candidates[
    candidates %in%
      colnames(object@meta.data)
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}


build_cell_embedding_plot_data <- function(
    object,
    dataset_id,
    reduction_name,
    cluster_column,
    annotation_column,
    rna_weight_column = NA_character_,
    adt_weight_column = NA_character_
) {
  
  embedding <- Embeddings(
    object,
    reduction_name
  )
  
  if (ncol(embedding) < 2) {
    stop(
      "Reduction ",
      reduction_name,
      " has fewer than two dimensions."
    )
  }
  
  embedding <- embedding[
    ,
    1:2,
    drop = FALSE
  ]
  
  cell_ids <- rownames(
    embedding
  )
  
  metadata <- object@meta.data[
    cell_ids,
    ,
    drop = FALSE
  ]
  
  sample_column <-
    select_first_existing_metadata_column(
      object,
      c(
        "sample",
        "sample_id",
        "Sample",
        "orig.ident",
        "donor",
        "donor_id"
      )
    )
  
  condition_column <-
    select_first_existing_metadata_column(
      object,
      c(
        "condition",
        "group",
        "treatment",
        "status",
        "timepoint"
      )
    )
  
  get_metadata_values <- function(
    column_name,
    fallback
  ) {
    
    if (
      is.na(column_name) ||
      !column_name %in%
      colnames(metadata)
    ) {
      return(
        rep(
          fallback,
          nrow(metadata)
        )
      )
    }
    
    as.character(
      metadata[[column_name]]
    )
  }
  
  get_weight_values <- function(
    column_name
  ) {
    
    if (
      is.na(column_name) ||
      !column_name %in%
      colnames(metadata)
    ) {
      return(
        rep(
          NA_real_,
          nrow(metadata)
        )
      )
    }
    
    as_numeric_safe(
      metadata[[column_name]]
    )
  }
  
  output <- data.frame(
    cell_id = cell_ids,
    UMAP_1 = as.numeric(
      embedding[, 1]
    ),
    UMAP_2 = as.numeric(
      embedding[, 2]
    ),
    cell_type = as.character(
      metadata$celltype_for_analysis
    ),
    cluster = as.character(
      metadata[[cluster_column]]
    ),
    sample = get_metadata_values(
      sample_column,
      dataset_id
    ),
    condition = get_metadata_values(
      condition_column,
      NA_character_
    ),
    dataset_id = dataset_id,
    annotation_source_column =
      annotation_column,
    cluster_source_column =
      cluster_column,
    reduction = reduction_name,
    RNA_weight = get_weight_values(
      rna_weight_column
    ),
    ADT_weight = get_weight_values(
      adt_weight_column
    ),
    stringsAsFactors = FALSE
  )
  
  output
}


build_plot_data_manifest <- function(
    output_dir,
    web_data_dir
) {
  
  figure_files <- list.files(
    output_dir,
    pattern = "\\.png$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(figure_files) == 0) {
    return(data.frame())
  }
  
  figure_relative <- vapply(
    figure_files,
    relative_path_safe,
    character(1),
    root = output_dir
  )
  
  web_relative <- relative_path_safe(
    web_data_dir,
    output_dir
  )
  
  rows <- vector(
    "list",
    length(figure_relative)
  )
  
  for (i in seq_along(figure_relative)) {
    
    current_figure <- normalize_path_for_table(
      figure_relative[i]
    )
    
    current_base <- basename(
      current_figure
    )
    
    data_file <- NA_character_
    filter_rule <- ""
    plot_type <- ""
    x_column <- ""
    y_column <- ""
    aesthetic_columns <- ""
    description <- ""
    
    if (
      current_base %in% c(
        "UMAP_annotated_celltypes.png",
        "UMAP_clusters.png",
        "WNN_modality_weights_on_UMAP.png"
      )
    ) {
      
      data_file <- file.path(
        web_relative,
        "cell_embedding_metadata_plot_data.csv"
      )
      
      x_column <- "UMAP_1"
      y_column <- "UMAP_2"
      
      if (
        current_base ==
        "UMAP_annotated_celltypes.png"
      ) {
        plot_type <- "UMAP scatter"
        aesthetic_columns <- "cell_type"
        description <-
          "Use cell_type for point color or grouping."
      } else if (
        current_base ==
        "UMAP_clusters.png"
      ) {
        plot_type <- "UMAP scatter"
        aesthetic_columns <- "cluster"
        description <-
          "Use cluster for point color or grouping."
      } else {
        plot_type <-
          "UMAP continuous feature panels"
        aesthetic_columns <-
          "RNA_weight; ADT_weight"
        description <-
          "Draw one panel for RNA_weight and one for ADT_weight."
      }
      
    } else if (
      current_base ==
      "cluster_by_celltype_heatmap.png"
    ) {
      
      data_file <- file.path(
        web_relative,
        "cluster_celltype_heatmap_plot_data.csv"
      )
      plot_type <- "Heatmap"
      x_column <- "CellType"
      y_column <- "Cluster"
      aesthetic_columns <- "Proportion"
      description <-
        "Fill tiles with within-cluster cell-type proportion."
      
    } else if (
      current_base %in% c(
        "RNA_marker_validation_dotplot.png",
        "ADT_marker_validation_dotplot.png",
        "ADT_top_markers_dotplot.png",
        "RNA_top_markers_dotplot.png"
      )
    ) {
      
      data_file <- file.path(
        web_relative,
        "marker_dotplot_all_plot_data.csv"
      )
      plot_type <- "Dot plot"
      filter_rule <- paste0(
        "plot_file == '",
        current_figure,
        "'"
      )
      x_column <- "feature"
      y_column <- "cell_type"
      aesthetic_columns <-
        "percent_expressing; avg_expression_scaled"
      description <-
        "Use percent_expressing as point size and avg_expression_scaled as point color."
      
    } else if (
      current_base ==
      "RNA_ADT_correlation_heatmap.png"
    ) {
      
      data_file <- file.path(
        web_relative,
        "RNA_ADT_correlation_heatmap_plot_data.csv"
      )
      plot_type <- "Heatmap"
      x_column <- "pair_label"
      y_column <- "celltype"
      aesthetic_columns <- "spearman_rho"
      description <-
        "Fill tiles with Spearman correlation coefficient."
      
    } else if (
      grepl(
        "^scatter_[0-9]+_",
        current_base
      )
    ) {
      
      data_file <- file.path(
        web_relative,
        "RNA_ADT_scatter_all_plot_data.csv"
      )
      plot_type <- "Scatter plot"
      filter_rule <- paste0(
        "plot_file == '",
        current_figure,
        "'"
      )
      x_column <- "RNA"
      y_column <- "ADT"
      aesthetic_columns <- "CellType"
      description <-
        "Use CellType for point color; each figure is selected with plot_file."
      
    } else if (
      current_base ==
      "WNN_modality_weight_violin.png"
    ) {
      
      data_file <- file.path(
        web_relative,
        "WNN_modality_weight_violin_plot_data.csv"
      )
      plot_type <- "Violin and box plot"
      x_column <- "CellType"
      y_column <- "Weight"
      aesthetic_columns <- "Modality"
      description <-
        "Facet or group by Modality."
      
    } else if (
      grepl(
        "enrichment",
        current_base,
        ignore.case = TRUE
      )
    ) {
      
      data_file <- file.path(
        web_relative,
        "RNA_DEG_GO_KEGG_enrichment_plot_data.csv"
      )
      plot_type <- "Enrichment dot plot"
      filter_rule <- paste0(
        "plot_file == '",
        current_figure,
        "'"
      )
      x_column <- "GeneRatio_numeric"
      y_column <- "term_label"
      aesthetic_columns <-
        "Count; minus_log10_p_adjust"
      description <-
        "Use Count as point size and minus_log10_p_adjust as point color."
      
    } else {
      
      plot_type <- "Unmapped figure"
      description <-
        "No explicit web plotting-data mapping was defined."
    }
    
    rows[[i]] <- data.frame(
      figure_file = current_figure,
      plot_type = plot_type,
      data_file = normalize_path_for_table(
        data_file
      ),
      filter_rule = filter_rule,
      x_column = x_column,
      y_column = y_column,
      aesthetic_columns =
        aesthetic_columns,
      description = description,
      stringsAsFactors = FALSE
    )
  }
  
  output <- do.call(
    rbind,
    rows
  )
  
  output <- output[
    order(output$figure_file),
    ,
    drop = FALSE
  ]
  
  rownames(output) <- NULL
  
  output
}


write_plot_data_readme <- function(
    file,
    manifest_file
) {
  
  write_text_safe(
    c(
      "WEB PLOTTING DATA",
      "",
      "1. cell_embedding_metadata_plot_data.csv",
      "   One row per cell. It combines UMAP coordinates, annotated cell type,",
      "   cluster, sample/dataset information and RNA/ADT modality weights.",
      "   It reproduces cell-type UMAP, cluster UMAP and modality-weight UMAP panels.",
      "",
      "2. marker_dotplot_all_plot_data.csv",
      "   Combined plotting data for RNA marker validation, ADT marker validation,",
      "   top differential ADT markers and top differential RNA markers.",
      "   Filter with plot_file before plotting.",
      "",
      "3. cluster_celltype_heatmap_plot_data.csv",
      "   Long-form cluster-by-cell-type row proportions.",
      "",
      "4. RNA_ADT_correlation_heatmap_plot_data.csv",
      "   Exact rows displayed in the RNA-ADT correlation heatmap.",
      "",
      "5. RNA_ADT_scatter_all_plot_data.csv",
      "   Combined cell-level data for every RNA-ADT scatter plot.",
      "   Filter with plot_file before plotting.",
      "",
      "6. WNN_modality_weight_violin_plot_data.csv",
      "   Long-form cell-level RNA and ADT weights for the violin plot.",
      "",
      "7. RNA_DEG_GO_KEGG_enrichment_plot_data.csv",
      "   Combined plotting rows for all GO and KEGG figures.",
      "   Filter with plot_file before plotting.",
      "",
      paste0(
        "The complete figure-to-table mapping is in: ",
        basename(manifest_file)
      )
    ),
    file
  )
}


is_valid_label <- function(x) {
  
  x <- trimws(as.character(x))
  
  !is.na(x) &
    x != "" &
    !toupper(x) %in% c(
      "NA",
      "N/A",
      "NULL",
      "UNASSIGNED",
      "UNKNOWN"
    )
}


downsample_cells_by_group <- function(
    object,
    group_column,
    maximum_per_group,
    seed = 1234
) {
  
  group_values <- as.character(
    object@meta.data[[group_column]]
  )
  
  names(group_values) <- rownames(
    object@meta.data
  )
  
  split_cells <- split(
    names(group_values),
    group_values
  )
  
  set.seed(seed)
  
  selected_cells <- unlist(
    lapply(
      split_cells,
      function(current_cells) {
        
        if (length(current_cells) <= maximum_per_group) {
          return(current_cells)
        }
        
        sample(
          current_cells,
          maximum_per_group
        )
      }
    ),
    use.names = FALSE
  )
  
  selected_cells
}


get_numeric_summary <- function(x) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(
      c(
        n = 0,
        mean = NA_real_,
        sd = NA_real_,
        median = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_,
        min = NA_real_,
        max = NA_real_
      )
    )
  }
  
  c(
    n = length(x),
    mean = mean(x),
    sd = if (length(x) > 1) sd(x) else NA_real_,
    median = median(x),
    q25 = unname(quantile(x, 0.25)),
    q75 = unname(quantile(x, 0.75)),
    min = min(x),
    max = max(x)
  )
}


## =========================================================
## 5. Seurat v4/v5 assay-layer helpers
## =========================================================

get_layer_names_safe <- function(
    object,
    assay_name
) {
  
  tryCatch(
    {
      Layers(object[[assay_name]])
    },
    error = function(e) {
      character(0)
    }
  )
}


join_assay_layers_if_needed <- function(
    object,
    assay_name
) {
  
  layer_names <- get_layer_names_safe(
    object,
    assay_name
  )
  
  if (length(layer_names) == 0) {
    return(object)
  }
  
  count_layers <- grep(
    "^counts($|\\.)",
    layer_names,
    value = TRUE
  )
  
  data_layers <- grep(
    "^data($|\\.)",
    layer_names,
    value = TRUE
  )
  
  needs_join <- (
    length(count_layers) > 1 ||
      length(data_layers) > 1 ||
      (
        length(data_layers) == 1 &&
          data_layers != "data"
      )
  )
  
  if (!needs_join) {
    return(object)
  }
  
  if (!exists("JoinLayers", mode = "function")) {
    
    warning(
      "Multiple layers were found in assay ",
      assay_name,
      ", but JoinLayers() is unavailable. ",
      "Some downstream functions may fail."
    )
    
    return(object)
  }
  
  message(
    "Joining layers for assay: ",
    assay_name
  )
  
  tryCatch(
    {
      JoinLayers(
        object = object,
        assay = assay_name
      )
    },
    error = function(e) {
      
      warning(
        "JoinLayers failed for assay ",
        assay_name,
        ": ",
        conditionMessage(e)
      )
      
      object
    }
  )
}


get_assay_matrix <- function(
    object,
    assay_name,
    layer_name = "data"
) {
  
  ## Preferred SeuratObject v5 accessor.
  result <- tryCatch(
    {
      LayerData(
        object = object,
        assay = assay_name,
        layer = layer_name
      )
    },
    error = function(e) {
      NULL
    }
  )
  
  if (!is.null(result)) {
    return(result)
  }
  
  ## Seurat v4/v5 compatibility fallback.
  result <- tryCatch(
    {
      GetAssayData(
        object = object,
        assay = assay_name,
        layer = layer_name
      )
    },
    error = function(e) {
      NULL
    }
  )
  
  if (!is.null(result)) {
    return(result)
  }
  
  result <- tryCatch(
    {
      GetAssayData(
        object = object,
        assay = assay_name,
        slot = layer_name
      )
    },
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(result)) {
    stop(
      "Unable to retrieve assay ",
      assay_name,
      " layer/slot ",
      layer_name,
      "."
    )
  }
  
  result
}


has_nonempty_data_layer <- function(
    object,
    assay_name
) {
  
  x <- tryCatch(
    {
      get_assay_matrix(
        object,
        assay_name,
        "data"
      )
    },
    error = function(e) {
      NULL
    }
  )
  
  !is.null(x) &&
    nrow(x) > 0 &&
    ncol(x) > 0
}


ensure_normalized_assays <- function(object) {
  
  if (!has_name_safe(get_assay_names_safe(object), "RNA")) {
    stop("RNA assay is missing.")
  }
  
  if (!has_name_safe(get_assay_names_safe(object), "ADT")) {
    stop("ADT assay is missing.")
  }
  
  object <- join_assay_layers_if_needed(
    object,
    "RNA"
  )
  
  object <- join_assay_layers_if_needed(
    object,
    "ADT"
  )
  
  if (!has_nonempty_data_layer(
    object,
    "RNA"
  )) {
    
    message(
      "RNA normalized data are missing; running NormalizeData."
    )
    
    object <- NormalizeData(
      object = object,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
  
  if (!has_nonempty_data_layer(
    object,
    "ADT"
  )) {
    
    message(
      "ADT normalized data are missing; running CLR normalization."
    )
    
    object <- NormalizeData(
      object = object,
      assay = "ADT",
      normalization.method = "CLR",
      margin = 2,
      verbose = FALSE
    )
  }
  
  object
}


## =========================================================
## 6. Annotation, cluster and reduction selection
## =========================================================

resolve_annotation_column <- function(object) {
  
  metadata_columns <- colnames(
    object@meta.data
  )
  
  usable_columns <- preferred_celltype_columns[
    preferred_celltype_columns %in%
      metadata_columns
  ]
  
  if (length(usable_columns) == 0) {
    stop(
      "No supported annotation column was found.\nAvailable columns:\n",
      paste(metadata_columns, collapse = ", ")
    )
  }
  
  selected_column <- NULL
  final_labels <- rep(
    NA_character_,
    ncol(object)
  )
  
  names(final_labels) <- colnames(object)
  
  for (current_column in usable_columns) {
    
    current_values <- as.character(
      object@meta.data[
        colnames(object),
        current_column
      ]
    )
    
    valid_current <- is_valid_label(
      current_values
    )
    
    if (
      is.null(selected_column) &&
      any(valid_current)
    ) {
      selected_column <- current_column
    }
    
    fill_positions <- (
      !is_valid_label(final_labels) &
        valid_current
    )
    
    final_labels[fill_positions] <-
      current_values[fill_positions]
  }
  
  if (is.null(selected_column)) {
    stop(
      "Annotation columns exist, but all labels are missing."
    )
  }
  
  final_labels[
    !is_valid_label(final_labels)
  ] <- "Unassigned"
  
  label_frequency <- sort(
    table(final_labels),
    decreasing = TRUE
  )
  
  object$celltype_for_analysis <- factor(
    final_labels[colnames(object)],
    levels = names(label_frequency)
  )
  
  list(
    object = object,
    selected_column = selected_column,
    source_columns = usable_columns
  )
}


resolve_cluster_column <- function(object) {
  
  metadata_columns <- colnames(
    object@meta.data
  )
  
  hit <- preferred_cluster_columns[
    preferred_cluster_columns %in%
      metadata_columns
  ]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  ## Fall back to current identities.
  object$cluster_for_analysis <- factor(
    as.character(Idents(object))
  )
  
  "cluster_for_analysis"
}


resolve_reduction <- function(object) {
  
  available_reductions <- get_reduction_names_safe(
    object
  )
  
  hit <- preferred_reductions[
    preferred_reductions %in%
      available_reductions
  ]
  
  if (length(hit) == 0) {
    stop(
      "No supported 2D reduction was found.\nAvailable reductions:\n",
      paste(available_reductions, collapse = ", ")
    )
  }
  
  hit[1]
}


## =========================================================
## 7. Built-in RNA marker panels
## =========================================================

get_rna_marker_panel <- function(species_name) {
  
  if (species_name == "mouse") {
    
    return(
      list(
        "T cell" = c(
          "Cd3d",
          "Cd3e",
          "Trac"
        ),
        "CD4 T" = c(
          "Cd4",
          "Il7r",
          "Ltb"
        ),
        "CD8 T" = c(
          "Cd8a",
          "Cd8b1",
          "Ccl5"
        ),
        "B cell" = c(
          "Cd79a",
          "Ms4a1",
          "Cd74",
          "Cd19"
        ),
        "Plasma cell" = c(
          "Jchain",
          "Mzb1",
          "Sdc1"
        ),
        "NK cell" = c(
          "Nkg7",
          "Klrd1",
          "Ncr1",
          "Gzmb"
        ),
        "Monocyte" = c(
          "Lyz2",
          "Ccr2",
          "Ly6c2",
          "S100a8",
          "S100a9"
        ),
        "Macrophage" = c(
          "Adgre1",
          "C1qa",
          "C1qb",
          "C1qc",
          "Fcgr1"
        ),
        "Neutrophil" = c(
          "Ly6g",
          "Csf3r",
          "S100a8",
          "S100a9"
        ),
        "Dendritic cell" = c(
          "Itgax",
          "Flt3",
          "H2-Aa",
          "H2-Ab1",
          "Xcr1"
        ),
        "Mast cell" = c(
          "Kit",
          "Fcer1a",
          "Cpa3",
          "Mcpt8"
        ),
        "Erythroid" = c(
          "Hbb-bs",
          "Hba-a1",
          "Gypa"
        ),
        "Platelet" = c(
          "Pf4",
          "Ppbp",
          "Itga2b"
        )
      )
    )
  }
  
  list(
    "T cell" = c(
      "CD3D",
      "CD3E",
      "TRAC"
    ),
    "CD4 T" = c(
      "CD4",
      "IL7R",
      "LTB"
    ),
    "CD8 T" = c(
      "CD8A",
      "CD8B",
      "CCL5"
    ),
    "B cell" = c(
      "CD79A",
      "MS4A1",
      "CD74",
      "CD19"
    ),
    "Plasma cell" = c(
      "JCHAIN",
      "MZB1",
      "SDC1"
    ),
    "NK cell" = c(
      "NKG7",
      "KLRD1",
      "NCR1",
      "GNLY"
    ),
    "Monocyte" = c(
      "LYZ",
      "CCR2",
      "S100A8",
      "S100A9",
      "FCN1"
    ),
    "Macrophage" = c(
      "C1QA",
      "C1QB",
      "C1QC",
      "FCGR1A"
    ),
    "Neutrophil" = c(
      "FCGR3B",
      "CSF3R",
      "S100A8",
      "S100A9"
    ),
    "Dendritic cell" = c(
      "ITGAX",
      "FCER1A",
      "HLA-DRA",
      "CLEC10A",
      "XCR1"
    ),
    "Mast cell" = c(
      "KIT",
      "TPSAB1",
      "CPA3"
    ),
    "Erythroid" = c(
      "HBB",
      "HBA1",
      "GYPA"
    ),
    "Platelet" = c(
      "PF4",
      "PPBP",
      "ITGA2B"
    )
  )
}


match_features_case_insensitive <- function(
    requested_features,
    available_features
) {
  
  feature_lookup <- setNames(
    available_features,
    toupper(available_features)
  )
  
  matched <- unname(
    feature_lookup[
      toupper(requested_features)
    ]
  )
  
  unique(
    matched[
      !is.na(matched)
    ]
  )
}


## =========================================================
## 8. ADT feature selection for validation DotPlot
## =========================================================

normalize_adt_label <- function(x) {
  
  x <- toupper(as.character(x))
  
  x <- gsub(
    "ANTI[-_ ]*(MOUSE|HUMAN)",
    " ",
    x
  )
  
  x <- gsub(
    "TOTALSEQ[-_ ]*[ABC]?",
    " ",
    x
  )
  
  x <- gsub(
    "ANTIBODY",
    " ",
    x
  )
  
  x <- gsub(
    "[^A-Z0-9]+",
    " ",
    x
  )
  
  trimws(
    gsub(
      "[[:space:]]+",
      " ",
      x
    )
  )
}


select_adt_features_for_dotplot <- function(
    adt_features,
    maximum_features
) {
  
  normalized <- normalize_adt_label(
    adt_features
  )
  
  priority_patterns <- c(
    "(^| )CD3($| )",
    "(^| )CD4($| )",
    "(^| )CD8A?($| )",
    "(^| )CD19($| )",
    "(^| )CD20($| )",
    "(^| )B220($| )",
    "(^| )NK1 1($| )",
    "(^| )NKP46($| )",
    "(^| )CD11B($| )",
    "(^| )CD11C($| )",
    "(^| )LY6C($| )",
    "(^| )LY6G($| )",
    "(^| )F4 80($| )",
    "(^| )CD64($| )",
    "(^| )MHC II($| )",
    "(^| )I A I E($| )",
    "(^| )CD44($| )",
    "(^| )CD62L($| )",
    "(^| )CD69($| )",
    "(^| )PD 1($| )",
    "(^| )PD1($| )",
    "(^| )CD25($| )",
    "(^| )CD127($| )",
    "(^| )KLRG1($| )",
    "(^| )CD86($| )",
    "(^| )CD80($| )",
    "(^| )CD40($| )",
    "(^| )CD24($| )",
    "(^| )CD38($| )",
    "(^| )CD138($| )",
    "(^| )CD103($| )",
    "(^| )CD117($| )",
    "(^| )SCA 1($| )",
    "(^| )TER119($| )"
  )
  
  selected_indices <- integer(0)
  
  for (current_pattern in priority_patterns) {
    
    current_hits <- grep(
      current_pattern,
      normalized,
      perl = TRUE
    )
    
    selected_indices <- unique(
      c(
        selected_indices,
        current_hits
      )
    )
    
    if (length(selected_indices) >= maximum_features) {
      break
    }
  }
  
  if (length(selected_indices) < maximum_features) {
    
    remaining <- setdiff(
      seq_along(adt_features),
      selected_indices
    )
    
    selected_indices <- c(
      selected_indices,
      head(
        remaining,
        maximum_features -
          length(selected_indices)
      )
    )
  }
  
  head(
    adt_features[selected_indices],
    maximum_features
  )
}


## =========================================================
## 9. RNA–protein matching dictionaries
## =========================================================

get_rna_adt_dictionary <- function(species_name) {
  
  if (species_name == "mouse") {
    
    return(
      data.frame(
        marker = c(
          "CD16/32",
          "B220",
          "MHC-II",
          "I-A/I-E",
          "F4/80",
          "NK1.1",
          "PD-L1",
          "PD-1",
          "TIM-3",
          "c-Kit",
          "Sca-1",
          "Siglec-H",
          "SIRP-alpha",
          "CD3",
          "CD4",
          "CD8",
          "CD19",
          "CD20",
          "CD45",
          "CD11b",
          "CD11c",
          "CD64",
          "Ly6C",
          "Ly6G",
          "CD44",
          "CD62L",
          "CD69",
          "CTLA-4",
          "LAG-3",
          "TIGIT",
          "ICOS",
          "CD25",
          "CD127",
          "KLRG1",
          "CD49b",
          "NKp46",
          "CD86",
          "CD80",
          "CD40",
          "CD24",
          "CD38",
          "CD138",
          "CD103",
          "CD90",
          "TER119",
          "CD117",
          "CXCR3",
          "CXCR4",
          "CCR2",
          "CCR7",
          "CD95",
          "CD28",
          "CD27",
          "CD43",
          "CD71",
          "CD150",
          "CD48",
          "CD34",
          "CD31",
          "CD41",
          "CD61",
          "CD14",
          "CD115",
          "CD206",
          "MARCO",
          "XCR1",
          "FCER1A"
        ),
        regex = c(
          "(^| )CD16 32($| )|(^| )CD1632($| )",
          "(^| )B220($| )|(^| )CD45R($| )",
          "(^| )MHC II($| )|(^| )MHCII($| )",
          "(^| )I A I E($| )|(^| )IAIE($| )",
          "(^| )F4 80($| )|(^| )F480($| )",
          "(^| )NK1 1($| )|(^| )NK11($| )",
          "(^| )PD L1($| )|(^| )PDL1($| )|(^| )CD274($| )",
          "(^| )PD 1($| )|(^| )PD1($| )|(^| )CD279($| )",
          "(^| )TIM 3($| )|(^| )TIM3($| )",
          "(^| )C KIT($| )|(^| )CKIT($| )",
          "(^| )SCA 1($| )|(^| )SCA1($| )",
          "(^| )SIGLEC H($| )|(^| )SIGLECH($| )",
          "(^| )SIRP A($| )|(^| )SIRPA($| )|(^| )CD172A($| )",
          "(^| )CD3($| )",
          "(^| )CD4($| )",
          "(^| )CD8($| )|(^| )CD8A($| )",
          "(^| )CD19($| )",
          "(^| )CD20($| )",
          "(^| )CD45($| )",
          "(^| )CD11B($| )",
          "(^| )CD11C($| )",
          "(^| )CD64($| )",
          "(^| )LY6C($| )",
          "(^| )LY6G($| )",
          "(^| )CD44($| )",
          "(^| )CD62L($| )",
          "(^| )CD69($| )",
          "(^| )CTLA 4($| )|(^| )CTLA4($| )",
          "(^| )LAG 3($| )|(^| )LAG3($| )",
          "(^| )TIGIT($| )",
          "(^| )ICOS($| )",
          "(^| )CD25($| )",
          "(^| )CD127($| )",
          "(^| )KLRG1($| )",
          "(^| )CD49B($| )",
          "(^| )NKP46($| )",
          "(^| )CD86($| )",
          "(^| )CD80($| )",
          "(^| )CD40($| )",
          "(^| )CD24($| )",
          "(^| )CD38($| )",
          "(^| )CD138($| )",
          "(^| )CD103($| )",
          "(^| )CD90($| )|(^| )THY1($| )",
          "(^| )TER119($| )",
          "(^| )CD117($| )",
          "(^| )CXCR3($| )",
          "(^| )CXCR4($| )",
          "(^| )CCR2($| )",
          "(^| )CCR7($| )",
          "(^| )CD95($| )|(^| )FAS($| )",
          "(^| )CD28($| )",
          "(^| )CD27($| )",
          "(^| )CD43($| )",
          "(^| )CD71($| )",
          "(^| )CD150($| )",
          "(^| )CD48($| )",
          "(^| )CD34($| )",
          "(^| )CD31($| )",
          "(^| )CD41($| )",
          "(^| )CD61($| )",
          "(^| )CD14($| )",
          "(^| )CD115($| )",
          "(^| )CD206($| )",
          "(^| )MARCO($| )",
          "(^| )XCR1($| )",
          "(^| )FCER1A($| )"
        ),
        genes = c(
          "Fcgr3;Fcgr2b",
          "Ptprc",
          "H2-Aa;H2-Ab1;H2-Eb1",
          "H2-Aa;H2-Ab1;H2-Eb1",
          "Adgre1",
          "Klrb1c",
          "Cd274",
          "Pdcd1",
          "Havcr2",
          "Kit",
          "Ly6a",
          "Siglech",
          "Sirpa",
          "Cd3d;Cd3e;Cd3g",
          "Cd4",
          "Cd8a;Cd8b1",
          "Cd19",
          "Ms4a1",
          "Ptprc",
          "Itgam",
          "Itgax",
          "Fcgr1",
          "Ly6c2",
          "Ly6g",
          "Cd44",
          "Sell",
          "Cd69",
          "Ctla4",
          "Lag3",
          "Tigit",
          "Icos",
          "Il2ra",
          "Il7r",
          "Klrg1",
          "Itga2",
          "Ncr1",
          "Cd86",
          "Cd80",
          "Cd40",
          "Cd24a",
          "Cd38",
          "Sdc1",
          "Itgae",
          "Thy1",
          "Ly76",
          "Kit",
          "Cxcr3",
          "Cxcr4",
          "Ccr2",
          "Ccr7",
          "Fas",
          "Cd28",
          "Cd27",
          "Spn",
          "Tfrc",
          "Slamf1",
          "Cd48",
          "Cd34",
          "Pecam1",
          "Itga2b",
          "Itgb3",
          "Cd14",
          "Csf1r",
          "Mrc1",
          "Marco",
          "Xcr1",
          "Fcer1a"
        ),
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    marker = c(
      "CD16/32",
      "MHC-II",
      "PD-L1",
      "PD-1",
      "TIM-3",
      "CD3",
      "CD4",
      "CD8",
      "CD19",
      "CD20",
      "CD45",
      "CD11b",
      "CD11c",
      "CD64",
      "CD14",
      "CD44",
      "CD62L",
      "CD69",
      "CTLA-4",
      "LAG-3",
      "TIGIT",
      "ICOS",
      "CD25",
      "CD127",
      "KLRG1",
      "CD49b",
      "NKp46",
      "CD86",
      "CD80",
      "CD40",
      "CD24",
      "CD38",
      "CD138",
      "CD103",
      "CD90",
      "CD117",
      "CXCR3",
      "CXCR4",
      "CCR2",
      "CCR7",
      "CD95",
      "CD28",
      "CD27",
      "CD43",
      "CD71",
      "CD150",
      "CD48",
      "CD34",
      "CD31",
      "CD41",
      "CD61",
      "CD115",
      "CD206",
      "MARCO",
      "XCR1",
      "FCER1A"
    ),
    regex = c(
      "(^| )CD16 32($| )|(^| )CD1632($| )",
      "(^| )MHC II($| )|(^| )MHCII($| )|(^| )HLA DR($| )",
      "(^| )PD L1($| )|(^| )PDL1($| )|(^| )CD274($| )",
      "(^| )PD 1($| )|(^| )PD1($| )|(^| )CD279($| )",
      "(^| )TIM 3($| )|(^| )TIM3($| )",
      "(^| )CD3($| )",
      "(^| )CD4($| )",
      "(^| )CD8($| )|(^| )CD8A($| )",
      "(^| )CD19($| )",
      "(^| )CD20($| )",
      "(^| )CD45($| )",
      "(^| )CD11B($| )",
      "(^| )CD11C($| )",
      "(^| )CD64($| )",
      "(^| )CD14($| )",
      "(^| )CD44($| )",
      "(^| )CD62L($| )",
      "(^| )CD69($| )",
      "(^| )CTLA 4($| )|(^| )CTLA4($| )",
      "(^| )LAG 3($| )|(^| )LAG3($| )",
      "(^| )TIGIT($| )",
      "(^| )ICOS($| )",
      "(^| )CD25($| )",
      "(^| )CD127($| )",
      "(^| )KLRG1($| )",
      "(^| )CD49B($| )",
      "(^| )NKP46($| )",
      "(^| )CD86($| )",
      "(^| )CD80($| )",
      "(^| )CD40($| )",
      "(^| )CD24($| )",
      "(^| )CD38($| )",
      "(^| )CD138($| )",
      "(^| )CD103($| )",
      "(^| )CD90($| )|(^| )THY1($| )",
      "(^| )CD117($| )|(^| )C KIT($| )",
      "(^| )CXCR3($| )",
      "(^| )CXCR4($| )",
      "(^| )CCR2($| )",
      "(^| )CCR7($| )",
      "(^| )CD95($| )|(^| )FAS($| )",
      "(^| )CD28($| )",
      "(^| )CD27($| )",
      "(^| )CD43($| )",
      "(^| )CD71($| )",
      "(^| )CD150($| )",
      "(^| )CD48($| )",
      "(^| )CD34($| )",
      "(^| )CD31($| )",
      "(^| )CD41($| )",
      "(^| )CD61($| )",
      "(^| )CD115($| )",
      "(^| )CD206($| )",
      "(^| )MARCO($| )",
      "(^| )XCR1($| )",
      "(^| )FCER1A($| )"
    ),
    genes = c(
      "FCGR3A;FCGR2A;FCGR2B",
      "HLA-DRA;HLA-DRB1;HLA-DPA1;HLA-DPB1",
      "CD274",
      "PDCD1",
      "HAVCR2",
      "CD3D;CD3E;CD3G",
      "CD4",
      "CD8A;CD8B",
      "CD19",
      "MS4A1",
      "PTPRC",
      "ITGAM",
      "ITGAX",
      "FCGR1A",
      "CD14",
      "CD44",
      "SELL",
      "CD69",
      "CTLA4",
      "LAG3",
      "TIGIT",
      "ICOS",
      "IL2RA",
      "IL7R",
      "KLRG1",
      "ITGA2",
      "NCR1",
      "CD86",
      "CD80",
      "CD40",
      "CD24",
      "CD38",
      "SDC1",
      "ITGAE",
      "THY1",
      "KIT",
      "CXCR3",
      "CXCR4",
      "CCR2",
      "CCR7",
      "FAS",
      "CD28",
      "CD27",
      "SPN",
      "TFRC",
      "SLAMF1",
      "CD48",
      "CD34",
      "PECAM1",
      "ITGA2B",
      "ITGB3",
      "CSF1R",
      "MRC1",
      "MARCO",
      "XCR1",
      "FCER1A"
    ),
    stringsAsFactors = FALSE
  )
}


build_rna_adt_pairs <- function(
    object,
    species_name
) {
  
  adt_features <- rownames(
    object[["ADT"]]
  )
  
  rna_features <- rownames(
    object[["RNA"]]
  )
  
  normalized_adt <- normalize_adt_label(
    adt_features
  )
  
  dictionary <- get_rna_adt_dictionary(
    species_name
  )
  
  rna_lookup <- setNames(
    rna_features,
    toupper(rna_features)
  )
  
  result_list <- list()
  result_index <- 0L
  
  for (feature_index in seq_along(adt_features)) {
    
    current_adt <- adt_features[feature_index]
    current_normalized <- normalized_adt[feature_index]
    
    matched_row <- NA_integer_
    
    for (dictionary_index in seq_len(nrow(dictionary))) {
      
      if (
        grepl(
          dictionary$regex[dictionary_index],
          current_normalized,
          perl = TRUE
        )
      ) {
        matched_row <- dictionary_index
        break
      }
    }
    
    marker_name <- NA_character_
    gene_candidates <- character(0)
    
    if (!is.na(matched_row)) {
      
      marker_name <- dictionary$marker[
        matched_row
      ]
      
      gene_candidates <- strsplit(
        dictionary$genes[matched_row],
        ";",
        fixed = TRUE
      )[[1]]
      
    } else {
      
      ## Fallback: direct case-insensitive match between a cleaned
      ## single-token ADT label and an RNA feature.
      token_candidates <- unique(
        unlist(
          strsplit(
            current_normalized,
            " ",
            fixed = TRUE
          )
        )
      )
      
      direct_hits <- unname(
        rna_lookup[
          token_candidates
        ]
      )
      
      direct_hits <- direct_hits[
        !is.na(direct_hits)
      ]
      
      if (length(direct_hits) > 0) {
        marker_name <- current_normalized
        gene_candidates <- direct_hits
      }
    }
    
    matched_genes <- unname(
      rna_lookup[
        toupper(gene_candidates)
      ]
    )
    
    matched_genes <- unique(
      matched_genes[
        !is.na(matched_genes)
      ]
    )
    
    if (length(matched_genes) == 0) {
      next
    }
    
    result_index <- result_index + 1L
    
    result_list[[result_index]] <- data.frame(
      marker = marker_name,
      adt_feature = current_adt,
      rna_genes = paste(
        matched_genes,
        collapse = ";"
      ),
      pair_label = paste0(
        marker_name,
        ": ",
        paste(
          matched_genes,
          collapse = "+"
        ),
        " vs ",
        current_adt
      ),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(result_list) == 0) {
    
    return(
      data.frame(
        marker = character(0),
        adt_feature = character(0),
        rna_genes = character(0),
        pair_label = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  pairs <- do.call(
    rbind,
    result_list
  )
  
  pairs <- pairs[
    !duplicated(
      paste(
        pairs$adt_feature,
        pairs$rna_genes,
        sep = "||"
      )
    ),
    ,
    drop = FALSE
  ]
  
  rownames(pairs) <- NULL
  
  pairs
}


## =========================================================
## 10. Marker-table helpers
## =========================================================

identify_marker_columns <- function(marker_table) {
  
  cluster_column <- intersect(
    c(
      "cluster",
      "ident",
      "group"
    ),
    colnames(marker_table)
  )
  
  fold_change_column <- intersect(
    c(
      "avg_log2FC",
      "avg_logFC",
      "avg_diff"
    ),
    colnames(marker_table)
  )
  
  if (length(cluster_column) == 0) {
    stop("No cluster column was found in the marker table.")
  }
  
  if (length(fold_change_column) == 0) {
    stop("No fold-change column was found in the marker table.")
  }
  
  list(
    cluster = cluster_column[1],
    fold_change = fold_change_column[1]
  )
}


top_markers_per_group <- function(
    marker_table,
    number_per_group
) {
  
  if (nrow(marker_table) == 0) {
    return(marker_table)
  }
  
  if (!"gene" %in% colnames(marker_table)) {
    marker_table$gene <- rownames(marker_table)
  }
  
  marker_columns <- identify_marker_columns(
    marker_table
  )
  
  split_markers <- split(
    marker_table,
    marker_table[[marker_columns$cluster]]
  )
  
  top_list <- lapply(
    split_markers,
    function(current_table) {
      
      ordering <- order(
        current_table[[marker_columns$fold_change]],
        decreasing = TRUE,
        na.last = TRUE
      )
      
      head(
        current_table[
          ordering,
          ,
          drop = FALSE
        ],
        number_per_group
      )
    }
  )
  
  output <- do.call(
    rbind,
    top_list
  )
  
  rownames(output) <- NULL
  
  output
}



## =========================================================
## 11. RNA DEG, GO and KEGG enrichment functions
## =========================================================

infer_annotation_keytype <- function(
    genes
) {
  
  genes <- unique(
    as.character(genes)
  )
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  if (length(genes) == 0) {
    return("SYMBOL")
  }
  
  ensembl_fraction <- mean(
    grepl(
      "^ENS[A-Z]*G[0-9]+",
      toupper(genes)
    )
  )
  
  if (is.finite(ensembl_fraction) &&
      ensembl_fraction >= 0.5) {
    "ENSEMBL"
  } else {
    "SYMBOL"
  }
}


map_genes_to_entrez <- function(
    genes,
    OrgDb,
    keytype
) {
  
  genes <- unique(
    as.character(genes)
  )
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  if (length(genes) == 0) {
    return(
      data.frame(
        input_gene = character(0),
        SYMBOL = character(0),
        ENTREZID = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  selected <- suppressWarnings(
    suppressMessages(
      AnnotationDbi::select(
        x = OrgDb,
        keys = genes,
        keytype = keytype,
        columns = unique(
          c(
            "SYMBOL",
            "ENTREZID"
          )
        )
      )
    )
  )
  
  if (nrow(selected) == 0) {
    return(
      data.frame(
        input_gene = character(0),
        SYMBOL = character(0),
        ENTREZID = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  selected$input_gene <- as.character(
    selected[[ keytype]]
  )
  
  if (!"SYMBOL" %in% colnames(selected)) {
    selected$SYMBOL <- NA_character_
  }
  
  if (!"ENTREZID" %in% colnames(selected)) {
    selected$ENTREZID <- NA_character_
  }
  
  selected <- selected[
    !is.na(selected$ENTREZID) &
      selected$ENTREZID != "",
    c(
      "input_gene",
      "SYMBOL",
      "ENTREZID"
    ),
    drop = FALSE
  ]
  
  selected$input_gene <- as.character(
    selected$input_gene
  )
  
  selected$SYMBOL <- as.character(
    selected$SYMBOL
  )
  
  selected$ENTREZID <- as.character(
    selected$ENTREZID
  )
  
  selected <- unique(
    selected
  )
  
  rownames(selected) <- NULL
  
  selected
}


parse_ratio_numeric <- function(x) {
  
  x <- as.character(x)
  
  vapply(
    x,
    function(current_value) {
      
      if (
        is.na(current_value) ||
        current_value == ""
      ) {
        return(NA_real_)
      }
      
      pieces <- strsplit(
        current_value,
        "/",
        fixed = TRUE
      )[[1]]
      
      if (length(pieces) != 2) {
        return(NA_real_)
      }
      
      numerator <- suppressWarnings(
        as.numeric(pieces[1])
      )
      
      denominator <- suppressWarnings(
        as.numeric(pieces[2])
      )
      
      if (
        !is.finite(numerator) ||
        !is.finite(denominator) ||
        denominator == 0
      ) {
        return(NA_real_)
      }
      
      numerator / denominator
    },
    numeric(1)
  )
}


standardize_enrichment_table <- function(
    enrichment_object,
    dataset_id,
    celltype,
    database,
    ontology,
    input_gene_count,
    mapped_gene_count
) {
  
  if (is.null(enrichment_object)) {
    return(data.frame())
  }
  
  enrichment_table <- tryCatch(
    {
      as.data.frame(
        enrichment_object
      )
    },
    error = function(e) {
      data.frame()
    }
  )
  
  if (nrow(enrichment_table) == 0) {
    return(data.frame())
  }
  
  get_enrichment_column <- function(
    column_name,
    default = NA
  ) {
    
    if (
      column_name %in%
      colnames(enrichment_table)
    ) {
      return(
        enrichment_table[[column_name]]
      )
    }
    
    rep(
      default,
      nrow(enrichment_table)
    )
  }
  
  output <- data.frame(
    dataset_id = dataset_id,
    celltype = celltype,
    database = database,
    ontology = ontology,
    ID = as.character(
      get_enrichment_column(
        "ID",
        NA_character_
      )
    ),
    Description = as.character(
      get_enrichment_column(
        "Description",
        NA_character_
      )
    ),
    GeneRatio = as.character(
      get_enrichment_column(
        "GeneRatio",
        NA_character_
      )
    ),
    GeneRatio_numeric = parse_ratio_numeric(
      get_enrichment_column(
        "GeneRatio",
        NA_character_
      )
    ),
    BgRatio = as.character(
      get_enrichment_column(
        "BgRatio",
        NA_character_
      )
    ),
    BgRatio_numeric = parse_ratio_numeric(
      get_enrichment_column(
        "BgRatio",
        NA_character_
      )
    ),
    pvalue = as_numeric_safe(
      get_enrichment_column("pvalue")
    ),
    p_adjust = as_numeric_safe(
      get_enrichment_column("p.adjust")
    ),
    qvalue = as_numeric_safe(
      get_enrichment_column("qvalue")
    ),
    Count = as.integer(
      as_numeric_safe(
        get_enrichment_column("Count")
      )
    ),
    geneID = as.character(
      get_enrichment_column(
        "geneID",
        NA_character_
      )
    ),
    input_gene_count =
      as.integer(input_gene_count),
    mapped_gene_count =
      as.integer(mapped_gene_count),
    stringsAsFactors = FALSE
  )
  
  output$minus_log10_p_adjust <-
    -log10(
      pmax(
        output$p_adjust,
        .Machine$double.xmin
      )
    )
  
  output
}


run_celltype_deg_enrichment <- function(
    marker_table,
    object,
    dataset_id,
    OrgDb,
    kegg_organism,
    go_ontologies,
    p_adjust_threshold,
    minimum_log2fc,
    minimum_pct,
    minimum_genes,
    maximum_genes,
    enrichment_pvalue_cutoff,
    enrichment_qvalue_cutoff
) {
  
  if (nrow(marker_table) == 0) {
    
    return(
      list(
        enrichment = data.frame(),
        gene_mapping = data.frame(),
        input_summary = data.frame()
      )
    )
  }
  
  if (!"gene" %in% colnames(marker_table)) {
    marker_table$gene <- rownames(
      marker_table
    )
  }
  
  marker_columns <- identify_marker_columns(
    marker_table
  )
  
  group_column <-
    marker_columns$cluster
  
  fold_change_column <-
    marker_columns$fold_change
  
  marker_table$gene <- as.character(
    marker_table$gene
  )
  
  marker_table$.fold_change <-
    as_numeric_safe(
      marker_table[[fold_change_column]]
    )
  
  keep <- (
    !is.na(marker_table$gene) &
      marker_table$gene != "" &
      is.finite(marker_table$.fold_change) &
      marker_table$.fold_change >=
      minimum_log2fc
  )
  
  if (
    "p_val_adj" %in%
    colnames(marker_table)
  ) {
    
    adjusted_values <- as_numeric_safe(
      marker_table$p_val_adj
    )
    
    keep <- keep &
      is.finite(adjusted_values) &
      adjusted_values <=
      p_adjust_threshold
  }
  
  if (
    "pct.1" %in%
    colnames(marker_table)
  ) {
    
    pct_values <- as_numeric_safe(
      marker_table$pct.1
    )
    
    keep <- keep &
      is.finite(pct_values) &
      pct_values >= minimum_pct
  }
  
  filtered_markers <- marker_table[
    keep,
    ,
    drop = FALSE
  ]
  
  if (nrow(filtered_markers) == 0) {
    
    return(
      list(
        enrichment = data.frame(),
        gene_mapping = data.frame(),
        input_summary = data.frame(
          dataset_id = dataset_id,
          celltype = character(0),
          input_marker_genes = integer(0),
          mapped_entrez_genes = integer(0),
          GO_terms = integer(0),
          KEGG_terms = integer(0),
          status = character(0),
          stringsAsFactors = FALSE
        )
      )
    )
  }
  
  filtered_markers[[group_column]] <- as.character(
    filtered_markers[[group_column]]
  )
  
  marker_split <- split(
    filtered_markers,
    filtered_markers[[group_column]]
  )
  
  selected_marker_tables <- lapply(
    marker_split,
    function(current_table) {
      
      ordering <- order(
        current_table$.fold_change,
        decreasing = TRUE,
        na.last = NA
      )
      
      current_table <- current_table[
        ordering,
        ,
        drop = FALSE
      ]
      
      current_table <- current_table[
        !duplicated(
          current_table$gene
        ),
        ,
        drop = FALSE
      ]
      
      head(
        current_table,
        maximum_genes
      )
    }
  )
  
  selected_marker_table <- do.call(
    rbind,
    selected_marker_tables
  )
  
  rownames(selected_marker_table) <- NULL
  
  all_input_genes <- unique(
    selected_marker_table$gene
  )
  
  keytype <- infer_annotation_keytype(
    rownames(object[["RNA"]])
  )
  
  marker_mapping <- map_genes_to_entrez(
    all_input_genes,
    OrgDb,
    keytype
  )
  
  universe_genes <- unique(
    rownames(object[["RNA"]])
  )
  
  universe_mapping <- map_genes_to_entrez(
    universe_genes,
    OrgDb,
    keytype
  )
  
  universe_entrez <- unique(
    universe_mapping$ENTREZID
  )
  
  enrichment_list <- list()
  mapping_list <- list()
  summary_list <- list()
  
  enrichment_index <- 0L
  mapping_index <- 0L
  summary_index <- 0L
  
  for (
    celltype_name in
    names(selected_marker_tables)
  ) {
    
    current_table <-
      selected_marker_tables[[celltype_name]]
    
    current_genes <- unique(
      current_table$gene
    )
    
    current_mapping <- marker_mapping[
      marker_mapping$input_gene %in%
        current_genes,
      ,
      drop = FALSE
    ]
    
    current_entrez <- unique(
      current_mapping$ENTREZID
    )
    
    mapping_index <- mapping_index + 1L
    
    if (nrow(current_mapping) > 0) {
      
      current_mapping$dataset_id <-
        dataset_id
      
      current_mapping$celltype <-
        celltype_name
      
      current_mapping$keytype <-
        keytype
      
      mapping_list[[mapping_index]] <-
        current_mapping[
          ,
          c(
            "dataset_id",
            "celltype",
            "keytype",
            "input_gene",
            "SYMBOL",
            "ENTREZID"
          ),
          drop = FALSE
        ]
      
    } else {
      
      mapping_list[[mapping_index]] <-
        data.frame()
    }
    
    go_term_count <- 0L
    kegg_term_count <- 0L
    current_status <- "COMPLETED"
    
    if (
      length(current_entrez) <
      minimum_genes
    ) {
      
      current_status <- paste0(
        "SKIPPED: only ",
        length(current_entrez),
        " mapped Entrez genes"
      )
      
    } else {
      
      for (
        current_ontology in
        go_ontologies
      ) {
        
        go_result <- tryCatch(
          {
            clusterProfiler::enrichGO(
              gene = current_entrez,
              universe = universe_entrez,
              OrgDb = OrgDb,
              keyType = "ENTREZID",
              ont = current_ontology,
              pAdjustMethod = "BH",
              pvalueCutoff =
                enrichment_pvalue_cutoff,
              qvalueCutoff =
                enrichment_qvalue_cutoff,
              readable = TRUE
            )
          },
          error = function(e) {
            
            warning(
              "GO enrichment failed for ",
              celltype_name,
              " (",
              current_ontology,
              "): ",
              conditionMessage(e)
            )
            
            NULL
          }
        )
        
        go_table <-
          standardize_enrichment_table(
            enrichment_object =
              go_result,
            dataset_id = dataset_id,
            celltype =
              celltype_name,
            database = paste0(
              "GO_",
              current_ontology
            ),
            ontology =
              current_ontology,
            input_gene_count =
              length(current_genes),
            mapped_gene_count =
              length(current_entrez)
          )
        
        if (nrow(go_table) > 0) {
          
          enrichment_index <-
            enrichment_index + 1L
          
          enrichment_list[[enrichment_index]] <- go_table
          
          go_term_count <- go_term_count +
            nrow(go_table)
        }
      }
      
      kegg_result <- tryCatch(
        {
          clusterProfiler::enrichKEGG(
            gene = current_entrez,
            universe = universe_entrez,
            organism =
              kegg_organism,
            keyType =
              "ncbi-geneid",
            pvalueCutoff =
              enrichment_pvalue_cutoff,
            pAdjustMethod = "BH",
            qvalueCutoff =
              enrichment_qvalue_cutoff,
            use_internal_data = FALSE
          )
        },
        error = function(e) {
          
          warning(
            "KEGG enrichment failed for ",
            celltype_name,
            ": ",
            conditionMessage(e),
            ". KEGG may require internet access."
          )
          
          NULL
        }
      )
      
      if (!is.null(kegg_result)) {
        
        kegg_result <- tryCatch(
          {
            clusterProfiler::setReadable(
              kegg_result,
              OrgDb = OrgDb,
              keyType = "ENTREZID"
            )
          },
          error = function(e) {
            kegg_result
          }
        )
      }
      
      kegg_table <-
        standardize_enrichment_table(
          enrichment_object =
            kegg_result,
          dataset_id = dataset_id,
          celltype =
            celltype_name,
          database = "KEGG",
          ontology = "KEGG",
          input_gene_count =
            length(current_genes),
          mapped_gene_count =
            length(current_entrez)
        )
      
      if (nrow(kegg_table) > 0) {
        
        enrichment_index <-
          enrichment_index + 1L
        
        enrichment_list[[enrichment_index]] <- kegg_table
        
        kegg_term_count <- nrow(
          kegg_table
        )
      }
    }
    
    summary_index <- summary_index + 1L
    
    summary_list[[summary_index]] <-
      data.frame(
        dataset_id = dataset_id,
        celltype = celltype_name,
        input_marker_genes =
          length(current_genes),
        mapped_entrez_genes =
          length(current_entrez),
        GO_terms =
          go_term_count,
        KEGG_terms =
          kegg_term_count,
        status =
          current_status,
        stringsAsFactors = FALSE
      )
  }
  
  list(
    enrichment =
      rbind_fill_data_frames(
        enrichment_list
      ),
    gene_mapping =
      rbind_fill_data_frames(
        mapping_list
      ),
    input_summary =
      rbind_fill_data_frames(
        summary_list
      )
  )
}


select_top_enrichment_terms <- function(
    enrichment_table,
    number_per_database
) {
  
  if (nrow(enrichment_table) == 0) {
    return(data.frame())
  }
  
  usable <- enrichment_table[
    is.finite(
      enrichment_table$p_adjust
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(usable) == 0) {
    return(data.frame())
  }
  
  split_key <- interaction(
    usable$celltype,
    usable$database,
    drop = TRUE,
    lex.order = TRUE
  )
  
  split_table <- split(
    usable,
    split_key
  )
  
  top_list <- lapply(
    split_table,
    function(current_table) {
      
      ordering <- order(
        current_table$p_adjust,
        -current_table$Count,
        na.last = TRUE
      )
      
      head(
        current_table[
          ordering,
          ,
          drop = FALSE
        ],
        number_per_database
      )
    }
  )
  
  output <- do.call(
    rbind,
    top_list
  )
  
  rownames(output) <- NULL
  
  output
}


prepare_enrichment_plot_data <- function(
    enrichment_table,
    database_selector,
    number_per_database,
    plot_file
) {
  
  if (nrow(enrichment_table) == 0) {
    return(data.frame())
  }
  
  selected <- if (
    identical(
      database_selector,
      "GO"
    )
  ) {
    
    enrichment_table[
      grepl(
        "^GO_",
        enrichment_table$database
      ),
      ,
      drop = FALSE
    ]
    
  } else if (
    identical(
      database_selector,
      "KEGG"
    )
  ) {
    
    enrichment_table[
      enrichment_table$database ==
        "KEGG",
      ,
      drop = FALSE
    ]
    
  } else {
    
    enrichment_table[
      enrichment_table$celltype ==
        database_selector,
      ,
      drop = FALSE
    ]
  }
  
  selected <- select_top_enrichment_terms(
    selected,
    number_per_database
  )
  
  if (nrow(selected) == 0) {
    return(data.frame())
  }
  
  selected$term_label <- paste0(
    "[",
    selected$database,
    "] ",
    selected$Description
  )
  
  selected$plot_file <-
    normalize_path_for_table(
      plot_file
    )
  
  selected$term_order <- ave(
    selected$p_adjust,
    selected$celltype,
    selected$database,
    FUN = function(x) {
      rank(
        x,
        ties.method = "first",
        na.last = "keep"
      )
    }
  )
  
  selected
}


order_enrichment_term_factor <- function(
    plot_data
) {
  
  if (nrow(plot_data) == 0) {
    return(plot_data)
  }
  
  plot_data$term_key <- paste(
    plot_data$term_label,
    plot_data$celltype,
    sep = "___CELLTYPE___"
  )
  
  split_data <- split(
    plot_data,
    plot_data$celltype
  )
  
  ordered_keys <- unlist(
    lapply(
      split_data,
      function(current_table) {
        
        current_table <- current_table[
          order(
            current_table$p_adjust,
            decreasing = TRUE,
            na.last = TRUE
          ),
          ,
          drop = FALSE
        ]
        
        current_table$term_key
      }
    ),
    use.names = FALSE
  )
  
  plot_data$term_key <- factor(
    plot_data$term_key,
    levels = unique(
      ordered_keys
    )
  )
  
  plot_data
}


make_enrichment_overview_plot <- function(
    plot_data,
    dataset_id,
    title_suffix,
    facet_columns = 3L
) {
  
  plot_data <- order_enrichment_term_factor(
    plot_data
  )
  
  ggplot(
    plot_data,
    aes(
      x = GeneRatio_numeric,
      y = term_key,
      size = Count,
      color = minus_log10_p_adjust
    )
  ) +
    geom_point(
      alpha = 0.85
    ) +
    facet_wrap(
      ~ celltype,
      scales = "free_y",
      ncol = facet_columns
    ) +
    scale_y_discrete(
      labels = function(x) {
        sub(
          "___CELLTYPE___.*$",
          "",
          x
        )
      }
    ) +
    labs(
      title = paste0(
        dataset_id,
        " — ",
        title_suffix
      ),
      x = "Gene ratio",
      y = "Enriched term",
      size = "Gene count",
      color = "-log10(adjusted P)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 7
      ),
      strip.text = element_text(
        face = "bold"
      )
    )
}


## =========================================================
## 12. RNA–protein correlation functions
## =========================================================

approximate_correlation_pvalue <- function(
    rho,
    n
) {
  
  if (
    is.na(rho) ||
    n < 3 ||
    abs(rho) >= 1
  ) {
    
    if (!is.na(rho) && abs(rho) == 1) {
      return(0)
    }
    
    return(NA_real_)
  }
  
  test_statistic <- rho *
    sqrt(
      (n - 2) /
        (1 - rho^2)
    )
  
  2 * pt(
    -abs(test_statistic),
    df = n - 2
  )
}


compute_rna_adt_correlations <- function(
    object,
    pair_table,
    group_column,
    minimum_cells,
    maximum_cells_per_group,
    include_all_cells = TRUE
) {
  
  if (nrow(pair_table) == 0) {
    
    return(
      data.frame(
        celltype = character(0),
        marker = character(0),
        adt_feature = character(0),
        rna_genes = character(0),
        pair_label = character(0),
        n_total = integer(0),
        n_used = integer(0),
        spearman_rho = numeric(0),
        p_value = numeric(0),
        p_adjusted = numeric(0),
        mean_rna = numeric(0),
        mean_adt = numeric(0),
        rna_positive_fraction = numeric(0),
        adt_positive_fraction = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  required_rna_genes <- unique(
    unlist(
      strsplit(
        pair_table$rna_genes,
        ";",
        fixed = TRUE
      )
    )
  )
  
  required_adt_features <- unique(
    pair_table$adt_feature
  )
  
  rna_matrix <- get_assay_matrix(
    object,
    "RNA",
    "data"
  )[
    required_rna_genes,
    ,
    drop = FALSE
  ]
  
  adt_matrix <- get_assay_matrix(
    object,
    "ADT",
    "data"
  )[
    required_adt_features,
    ,
    drop = FALSE
  ]
  
  group_values <- as.character(
    object@meta.data[
      colnames(object),
      group_column
    ]
  )
  
  names(group_values) <- colnames(object)
  
  group_cells <- split(
    names(group_values),
    group_values
  )
  
  if (include_all_cells) {
    group_cells <- c(
      list("All cells" = colnames(object)),
      group_cells
    )
  }
  
  result_list <- list()
  result_index <- 0L
  
  set.seed(1234)
  
  for (group_name in names(group_cells)) {
    
    current_cells_all <- group_cells[[group_name]]
    
    n_total <- length(
      current_cells_all
    )
    
    if (n_total < minimum_cells) {
      next
    }
    
    current_cells <- current_cells_all
    
    if (
      length(current_cells) >
      maximum_cells_per_group
    ) {
      
      current_cells <- sample(
        current_cells,
        maximum_cells_per_group
      )
    }
    
    for (pair_index in seq_len(nrow(pair_table))) {
      
      current_genes <- strsplit(
        pair_table$rna_genes[pair_index],
        ";",
        fixed = TRUE
      )[[1]]
      
      current_adt <- pair_table$adt_feature[
        pair_index
      ]
      
      if (length(current_genes) == 1) {
        
        rna_values <- as.numeric(
          rna_matrix[
            current_genes,
            current_cells,
            drop = TRUE
          ]
        )
        
      } else {
        
        rna_values <- as.numeric(
          Matrix::colMeans(
            rna_matrix[
              current_genes,
              current_cells,
              drop = FALSE
            ]
          )
        )
      }
      
      adt_values <- as.numeric(
        adt_matrix[
          current_adt,
          current_cells,
          drop = TRUE
        ]
      )
      
      valid_positions <- (
        is.finite(rna_values) &
          is.finite(adt_values)
      )
      
      rna_values <- rna_values[
        valid_positions
      ]
      
      adt_values <- adt_values[
        valid_positions
      ]
      
      n_used <- length(rna_values)
      
      if (
        n_used < minimum_cells ||
        length(unique(rna_values)) < 2 ||
        length(unique(adt_values)) < 2
      ) {
        rho <- NA_real_
        p_value <- NA_real_
      } else {
        
        rho <- suppressWarnings(
          cor(
            rna_values,
            adt_values,
            method = "spearman",
            use = "complete.obs"
          )
        )
        
        p_value <- approximate_correlation_pvalue(
          rho,
          n_used
        )
      }
      
      result_index <- result_index + 1L
      
      result_list[[result_index]] <- data.frame(
        celltype = group_name,
        marker = pair_table$marker[pair_index],
        adt_feature = current_adt,
        rna_genes = pair_table$rna_genes[pair_index],
        pair_label = pair_table$pair_label[pair_index],
        n_total = n_total,
        n_used = n_used,
        spearman_rho = rho,
        p_value = p_value,
        p_adjusted = NA_real_,
        mean_rna = if (n_used > 0) {
          mean(rna_values)
        } else {
          NA_real_
        },
        mean_adt = if (n_used > 0) {
          mean(adt_values)
        } else {
          NA_real_
        },
        rna_positive_fraction = if (n_used > 0) {
          mean(rna_values > 0)
        } else {
          NA_real_
        },
        adt_positive_fraction = if (n_used > 0) {
          mean(adt_values > 0)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(result_list) == 0) {
    return(data.frame())
  }
  
  result <- do.call(
    rbind,
    result_list
  )
  
  result$p_adjusted <- p.adjust(
    result$p_value,
    method = "BH"
  )
  
  rownames(result) <- NULL
  
  result
}


rank_correlation_pairs <- function(
    correlation_table
) {
  
  usable <- correlation_table[
    correlation_table$celltype != "All cells" &
      is.finite(
        correlation_table$spearman_rho
      ),
    ,
    drop = FALSE
  ]
  
  if (nrow(usable) == 0) {
    
    usable <- correlation_table[
      is.finite(
        correlation_table$spearman_rho
      ),
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(usable) == 0) {
    return(character(0))
  }
  
  pair_split <- split(
    abs(usable$spearman_rho),
    usable$pair_label
  )
  
  pair_scores <- vapply(
    pair_split,
    median,
    numeric(1),
    na.rm = TRUE
  )
  
  names(
    sort(
      pair_scores,
      decreasing = TRUE
    )
  )
}


make_pair_scatter_plot <- function(
    object,
    pair_row,
    group_column,
    maximum_cells
) {
  
  all_cells <- colnames(object)
  
  set.seed(1234)
  
  plot_cells <- if (
    length(all_cells) > maximum_cells
  ) {
    sample(
      all_cells,
      maximum_cells
    )
  } else {
    all_cells
  }
  
  current_genes <- strsplit(
    pair_row$rna_genes,
    ";",
    fixed = TRUE
  )[[1]]
  
  rna_matrix <- get_assay_matrix(
    object,
    "RNA",
    "data"
  )[
    current_genes,
    plot_cells,
    drop = FALSE
  ]
  
  adt_matrix <- get_assay_matrix(
    object,
    "ADT",
    "data"
  )[
    pair_row$adt_feature,
    plot_cells,
    drop = FALSE
  ]
  
  rna_values <- if (
    length(current_genes) == 1
  ) {
    
    as.numeric(
      rna_matrix[
        1,
        ,
        drop = TRUE
      ]
    )
    
  } else {
    
    as.numeric(
      Matrix::colMeans(
        rna_matrix
      )
    )
  }
  
  adt_values <- as.numeric(
    adt_matrix[
      1,
      ,
      drop = TRUE
    ]
  )
  
  plot_data <- data.frame(
    cell_id = plot_cells,
    RNA = rna_values,
    ADT = adt_values,
    CellType = as.character(
      object@meta.data[
        plot_cells,
        group_column
      ]
    ),
    stringsAsFactors = FALSE
  )
  
  valid <- (
    is.finite(plot_data$RNA) &
      is.finite(plot_data$ADT)
  )
  
  plot_data <- plot_data[
    valid,
    ,
    drop = FALSE
  ]
  
  rho <- suppressWarnings(
    cor(
      plot_data$RNA,
      plot_data$ADT,
      method = "spearman"
    )
  )
  
  ggplot(
    plot_data,
    aes(
      x = RNA,
      y = ADT,
      color = CellType
    )
  ) +
    geom_point(
      size = 0.55,
      alpha = 0.45
    ) +
    labs(
      title = paste0(
        pair_row$marker,
        " RNA–protein relationship"
      ),
      subtitle = paste0(
        "Spearman rho = ",
        format(
          round(rho, 3),
          nsmall = 3
        ),
        "; cells plotted = ",
        nrow(plot_data)
      ),
      x = paste0(
        "RNA normalized expression: ",
        gsub(
          ";",
          " + ",
          pair_row$rna_genes,
          fixed = TRUE
        )
      ),
      y = paste0(
        "ADT CLR expression: ",
        pair_row$adt_feature
      ),
      color = "Cell type"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "right"
    )
}


## =========================================================
## 13. WNN modality-weight analysis
## =========================================================

detect_modality_weight_columns <- function(object) {
  
  metadata_columns <- colnames(
    object@meta.data
  )
  
  rna_candidates <- metadata_columns[
    grepl(
      "RNA.*weight|weight.*RNA",
      metadata_columns,
      ignore.case = TRUE
    )
  ]
  
  adt_candidates <- metadata_columns[
    grepl(
      "ADT.*weight|weight.*ADT|protein.*weight|weight.*protein",
      metadata_columns,
      ignore.case = TRUE
    )
  ]
  
  rna_column <- if (
    "RNA.weight" %in% metadata_columns
  ) {
    "RNA.weight"
  } else if (
    length(rna_candidates) > 0
  ) {
    rna_candidates[1]
  } else {
    NA_character_
  }
  
  adt_column <- if (
    "ADT.weight" %in% metadata_columns
  ) {
    "ADT.weight"
  } else if (
    length(adt_candidates) > 0
  ) {
    adt_candidates[1]
  } else {
    NA_character_
  }
  
  if (
    !is.na(rna_column) &&
    is.na(adt_column)
  ) {
    
    rna_values <- suppressWarnings(
      as.numeric(
        as.character(
          object@meta.data[[rna_column]]
        )
      )
    )
    
    finite_values <- rna_values[
      is.finite(rna_values)
    ]
    
    if (
      length(finite_values) > 0 &&
      min(finite_values) >= -1e-6 &&
      max(finite_values) <= 1 + 1e-6
    ) {
      
      object$ADT.weight.derived <-
        1 - rna_values
      
      adt_column <- "ADT.weight.derived"
    }
  }
  
  list(
    object = object,
    rna_column = rna_column,
    adt_column = adt_column
  )
}


summarize_modality_weights <- function(
    object,
    group_column,
    rna_weight_column,
    adt_weight_column
) {
  
  available <- c(
    rna_weight_column,
    adt_weight_column
  )
  
  available <- available[
    !is.na(available)
  ]
  
  if (length(available) == 0) {
    return(data.frame())
  }
  
  modality_names <- c(
    RNA = rna_weight_column,
    ADT = adt_weight_column
  )
  
  modality_names <- modality_names[
    !is.na(modality_names)
  ]
  
  group_values <- as.character(
    object@meta.data[[group_column]]
  )
  
  result_list <- list()
  result_index <- 0L
  
  for (modality_name in names(modality_names)) {
    
    current_column <- modality_names[[modality_name]]
    
    current_values <- suppressWarnings(
      as.numeric(
        as.character(
          object@meta.data[[current_column]]
        )
      )
    )
    
    split_values <- split(
      current_values,
      group_values
    )
    
    for (group_name in names(split_values)) {
      
      summary_values <- get_numeric_summary(
        split_values[[group_name]]
      )
      
      result_index <- result_index + 1L
      
      result_list[[result_index]] <- data.frame(
        celltype = group_name,
        modality = modality_name,
        source_column = current_column,
        n = unname(summary_values["n"]),
        mean = unname(summary_values["mean"]),
        sd = unname(summary_values["sd"]),
        median = unname(summary_values["median"]),
        q25 = unname(summary_values["q25"]),
        q75 = unname(summary_values["q75"]),
        min = unname(summary_values["min"]),
        max = unname(summary_values["max"]),
        stringsAsFactors = FALSE
      )
    }
  }
  
  output <- do.call(
    rbind,
    result_list
  )
  
  rownames(output) <- NULL
  
  output
}


make_modality_weight_long_data <- function(
    object,
    group_column,
    rna_weight_column,
    adt_weight_column,
    maximum_per_group
) {
  
  selected_cells <- downsample_cells_by_group(
    object,
    group_column,
    maximum_per_group
  )
  
  result_list <- list()
  
  if (!is.na(rna_weight_column)) {
    
    result_list[["RNA"]] <- data.frame(
      Cell = selected_cells,
      CellType = as.character(
        object@meta.data[
          selected_cells,
          group_column
        ]
      ),
      Modality = "RNA",
      Weight = suppressWarnings(
        as.numeric(
          as.character(
            object@meta.data[
              selected_cells,
              rna_weight_column
            ]
          )
        )
      ),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.na(adt_weight_column)) {
    
    result_list[["ADT"]] <- data.frame(
      Cell = selected_cells,
      CellType = as.character(
        object@meta.data[
          selected_cells,
          group_column
        ]
      ),
      Modality = "ADT",
      Weight = suppressWarnings(
        as.numeric(
          as.character(
            object@meta.data[
              selected_cells,
              adt_weight_column
            ]
          )
        )
      ),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(result_list) == 0) {
    return(data.frame())
  }
  
  output <- do.call(
    rbind,
    result_list
  )
  
  output <- output[
    is.finite(output$Weight),
    ,
    drop = FALSE
  ]
  
  rownames(output) <- NULL
  
  output
}


## =========================================================
## 14. Main analysis for one annotated RDS
## =========================================================

analyze_one_annotated_rds <- function(
    input_rds_file
) {
  
  gse_id <- extract_gse_id(
    input_rds_file
  )
  
  gse_folder <- dirname(
    input_rds_file
  )
  
  output_dir <- file.path(
    gse_folder,
    output_folder_name
  )
  
  completion_flag <- file.path(
    output_dir,
    "analysis_complete.flag"
  )
  
  if (
    file.exists(completion_flag) &&
    !overwrite_completed
  ) {
    
    message(
      "[SKIP] ",
      gse_id,
      " already completed."
    )
    
    return(
      data.frame(
        GSE = gse_id,
        status = "SKIPPED_COMPLETED",
        input_rds = input_rds_file,
        output_dir = output_dir,
        cells = NA_integer_,
        RNA_features = NA_integer_,
        ADT_features = NA_integer_,
        annotation_column = NA_character_,
        cluster_column = NA_character_,
        reduction = NA_character_,
        matched_RNA_ADT_pairs = NA_integer_,
        message = "Completion flag already exists.",
        stringsAsFactors = FALSE
      )
    )
  }
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  annotation_dir <- file.path(
    output_dir,
    "01_annotation_validation"
  )
  
  adt_marker_dir <- file.path(
    output_dir,
    "02_ADT_markers"
  )
  
  correlation_dir <- file.path(
    output_dir,
    "03_RNA_ADT_correlation"
  )
  
  weight_dir <- file.path(
    output_dir,
    "04_WNN_modality_weights"
  )
  
  enrichment_dir <- file.path(
    output_dir,
    "05_RNA_DEG_GO_KEGG"
  )
  
  web_data_dir <- file.path(
    output_dir,
    web_plot_data_folder_name
  )
  
  dir.create(
    annotation_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    adt_marker_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    correlation_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    weight_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    enrichment_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    web_data_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dotplot_data_list <- list()
  scatter_plot_data_list <- list()
  enrichment_plot_data_list <- list()
  
  log_file <- file.path(
    output_dir,
    "analysis_log.txt"
  )
  
  write_text_safe(
    c(
      paste0("GSE: ", gse_id),
      paste0("Input RDS: ", input_rds_file),
      paste0("Started: ", timestamp_now()),
      ""
    ),
    log_file
  )
  
  result <- tryCatch(
    {
      
      message(
        "\n",
        paste(rep("=", 72), collapse = ""),
        "\n[START] ",
        gse_id,
        "\nInput: ",
        input_rds_file
      )
      
      obj <- readRDS(
        input_rds_file
      )
      
      if (!inherits(obj, "Seurat")) {
        stop(
          "The RDS does not contain a Seurat object. Class: ",
          paste(class(obj), collapse = ", ")
        )
      }
      
      assay_names <- get_assay_names_safe(
        obj
      )
      
      if (!has_all_names_safe(
        assay_names,
        c("RNA", "ADT")
      )) {
        stop(
          "Both RNA and ADT assays are required. Available assays: ",
          paste(assay_names, collapse = ", ")
        )
      }
      
      obj <- ensure_normalized_assays(
        obj
      )
      
      annotation_result <- resolve_annotation_column(
        obj
      )
      
      obj <- annotation_result$object
      annotation_column <- annotation_result$selected_column
      
      cluster_column <- resolve_cluster_column(
        obj
      )
      
      if (
        cluster_column == "cluster_for_analysis"
      ) {
        
        obj$cluster_for_analysis <- factor(
          as.character(Idents(obj))
        )
      }
      
      reduction_name <- resolve_reduction(
        obj
      )
      
      celltype_counts <- sort(
        table(
          obj$celltype_for_analysis
        ),
        decreasing = TRUE
      )
      
      cluster_counts <- sort(
        table(
          obj@meta.data[[cluster_column]]
        ),
        decreasing = TRUE
      )
      
      object_summary <- data.frame(
        GSE = gse_id,
        input_rds = input_rds_file,
        cells = ncol(obj),
        RNA_features = nrow(obj[["RNA"]]),
        ADT_features = nrow(obj[["ADT"]]),
        selected_annotation_column = annotation_column,
        annotation_source_columns = paste(
          annotation_result$source_columns,
          collapse = ";"
        ),
        cluster_column = cluster_column,
        reduction = reduction_name,
        stringsAsFactors = FALSE
      )
      
      write_csv_safe(
        object_summary,
        file.path(
          output_dir,
          "object_summary.csv"
        )
      )
      
      write_csv_safe(
        data.frame(
          celltype = names(celltype_counts),
          cells = as.integer(celltype_counts),
          fraction = as.numeric(
            celltype_counts /
              sum(celltype_counts)
          ),
          stringsAsFactors = FALSE
        ),
        file.path(
          annotation_dir,
          "celltype_counts.csv"
        )
      )
      
      write_csv_safe(
        data.frame(
          cluster = names(cluster_counts),
          cells = as.integer(cluster_counts),
          stringsAsFactors = FALSE
        ),
        file.path(
          annotation_dir,
          "cluster_counts.csv"
        )
      )
      
      
      ## ---------------------------------------------------
      ## 13.1 Annotation UMAP
      ## ---------------------------------------------------
      
      p_annotation <- DimPlot(
        object = obj,
        reduction = reduction_name,
        group.by = "celltype_for_analysis",
        label = TRUE,
        repel = TRUE,
        raster = NULL
      ) +
        ggtitle(
          paste0(
            gse_id,
            " — annotated cell types"
          )
        ) +
        theme(
          plot.title = element_text(
            face = "bold"
          )
        )
      
      save_plot_safe(
        p_annotation,
        file.path(
          annotation_dir,
          "UMAP_annotated_celltypes.png"
        ),
        width = 11,
        height = 8
      )
      
      p_clusters <- DimPlot(
        object = obj,
        reduction = reduction_name,
        group.by = cluster_column,
        label = TRUE,
        repel = TRUE,
        raster = NULL
      ) +
        ggtitle(
          paste0(
            gse_id,
            " — WNN clusters"
          )
        ) +
        theme(
          plot.title = element_text(
            face = "bold"
          )
        )
      
      save_plot_safe(
        p_clusters,
        file.path(
          annotation_dir,
          "UMAP_clusters.png"
        ),
        width = 10,
        height = 8
      )
      
      
      ## ---------------------------------------------------
      ## 13.2 Cluster × annotation table and heatmap
      ## ---------------------------------------------------
      
      cluster_values <- as.character(
        obj@meta.data[[cluster_column]]
      )
      
      annotation_values <- as.character(
        obj$celltype_for_analysis
      )
      
      cluster_annotation_table <- table(
        cluster = cluster_values,
        celltype = annotation_values
      )
      
      cluster_annotation_proportion <- prop.table(
        cluster_annotation_table,
        margin = 1
      )
      
      write_csv_safe(
        as.data.frame.matrix(
          cluster_annotation_table
        ),
        file.path(
          annotation_dir,
          "cluster_by_celltype_counts.csv"
        ),
        row.names = TRUE
      )
      
      write_csv_safe(
        as.data.frame.matrix(
          cluster_annotation_proportion
        ),
        file.path(
          annotation_dir,
          "cluster_by_celltype_row_proportions.csv"
        ),
        row.names = TRUE
      )
      
      heatmap_data <- as.data.frame(
        cluster_annotation_proportion
      )
      
      colnames(heatmap_data) <- c(
        "Cluster",
        "CellType",
        "Proportion"
      )
      
      p_cluster_heatmap <- ggplot(
        heatmap_data,
        aes(
          x = CellType,
          y = Cluster,
          fill = Proportion
        )
      ) +
        geom_tile() +
        scale_fill_gradient(
          low = "white",
          high = "steelblue"
        ) +
        labs(
          title = paste0(
            gse_id,
            " — cluster-to-annotation agreement"
          ),
          x = "Annotated cell type",
          y = "WNN cluster",
          fill = "Row\nproportion"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(
            face = "bold"
          ),
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          )
        )
      
      save_plot_safe(
        p_cluster_heatmap,
        file.path(
          annotation_dir,
          "cluster_by_celltype_heatmap.png"
        ),
        width = max(
          9,
          0.45 *
            nlevels(
              obj$celltype_for_analysis
            ) +
            5
        ),
        height = max(
          6,
          0.35 *
            length(
              unique(cluster_values)
            ) +
            3
        )
      )
      
      if (isTRUE(export_web_plot_data)) {
        
        heatmap_plot_data <- heatmap_data
        heatmap_plot_data$dataset_id <- gse_id
        heatmap_plot_data$Cluster_order <-
          factor_order_index(
            heatmap_plot_data$Cluster
          )
        heatmap_plot_data$CellType_order <-
          factor_order_index(
            heatmap_plot_data$CellType
          )
        heatmap_plot_data$plot_file <-
          normalize_path_for_table(
            file.path(
              "01_annotation_validation",
              "cluster_by_celltype_heatmap.png"
            )
          )
        
        heatmap_plot_data <- heatmap_plot_data[
          ,
          c(
            "dataset_id",
            "plot_file",
            "Cluster",
            "Cluster_order",
            "CellType",
            "CellType_order",
            "Proportion"
          ),
          drop = FALSE
        ]
        
        write_csv_safe(
          heatmap_plot_data,
          file.path(
            web_data_dir,
            "cluster_celltype_heatmap_plot_data.csv"
          )
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.3 Known RNA marker DotPlot
      ## ---------------------------------------------------
      
      rna_marker_panel <- get_rna_marker_panel(
        species
      )
      
      all_requested_rna_markers <- unique(
        unlist(
          rna_marker_panel,
          use.names = FALSE
        )
      )
      
      available_rna_markers <-
        match_features_case_insensitive(
          all_requested_rna_markers,
          rownames(obj[["RNA"]])
        )
      
      if (length(available_rna_markers) > 0) {
        
        DefaultAssay(obj) <- "RNA"
        
        p_rna_dot <- DotPlot(
          object = obj,
          features = available_rna_markers,
          group.by = "celltype_for_analysis",
          assay = "RNA",
          scale = TRUE
        ) +
          RotatedAxis() +
          labs(
            title = paste0(
              gse_id,
              " — RNA marker validation"
            ),
            x = "RNA marker",
            y = "Annotated cell type"
          ) +
          theme(
            plot.title = element_text(
              face = "bold"
            )
          )
        
        if (isTRUE(export_web_plot_data)) {
          
          dotplot_data_list[[length(dotplot_data_list) + 1L]] <- standardize_dotplot_data(
            plot_object = p_rna_dot,
            dataset_id = gse_id,
            plot_id =
              "RNA_marker_validation",
            plot_file = file.path(
              "01_annotation_validation",
              "RNA_marker_validation_dotplot.png"
            ),
            modality = "RNA"
          )
        }
        
        save_plot_safe(
          p_rna_dot,
          file.path(
            annotation_dir,
            "RNA_marker_validation_dotplot.png"
          ),
          width = max(
            12,
            0.35 *
              length(
                available_rna_markers
              ) +
              5
          ),
          height = max(
            7,
            0.35 *
              nlevels(
                obj$celltype_for_analysis
              ) +
              3
          )
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.4 ADT validation DotPlot
      ## ---------------------------------------------------
      
      adt_dot_features <-
        select_adt_features_for_dotplot(
          rownames(obj[["ADT"]]),
          adt_top_features_for_dotplot
        )
      
      if (length(adt_dot_features) > 0) {
        
        DefaultAssay(obj) <- "ADT"
        
        p_adt_dot <- DotPlot(
          object = obj,
          features = adt_dot_features,
          group.by = "celltype_for_analysis",
          assay = "ADT",
          scale = TRUE
        ) +
          RotatedAxis() +
          labs(
            title = paste0(
              gse_id,
              " — ADT marker validation"
            ),
            x = "ADT marker",
            y = "Annotated cell type"
          ) +
          theme(
            plot.title = element_text(
              face = "bold"
            )
          )
        
        if (isTRUE(export_web_plot_data)) {
          
          dotplot_data_list[[length(dotplot_data_list) + 1L]] <- standardize_dotplot_data(
            plot_object = p_adt_dot,
            dataset_id = gse_id,
            plot_id =
              "ADT_marker_validation",
            plot_file = file.path(
              "01_annotation_validation",
              "ADT_marker_validation_dotplot.png"
            ),
            modality = "ADT"
          )
        }
        
        save_plot_safe(
          p_adt_dot,
          file.path(
            annotation_dir,
            "ADT_marker_validation_dotplot.png"
          ),
          width = max(
            12,
            0.38 *
              length(
                adt_dot_features
              ) +
              5
          ),
          height = max(
            7,
            0.35 *
              nlevels(
                obj$celltype_for_analysis
              ) +
              3
          )
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.5 ADT markers by annotated cell type
      ## ---------------------------------------------------
      
      group_sizes <- table(
        obj$celltype_for_analysis
      )
      
      eligible_groups <- names(
        group_sizes[group_sizes >=min_cells_for_marker_group]
      )
      
      marker_cells <- colnames(obj)[
        as.character(
          obj$celltype_for_analysis
        ) %in%
          eligible_groups
      ]
      
      adt_markers <- data.frame()
      
      if (
        length(eligible_groups) >= 2 &&
        length(marker_cells) > 0
      ) {
        
        marker_object <- subset(
          obj,
          cells = marker_cells
        )
        
        Idents(marker_object) <-
          "celltype_for_analysis"
        
        DefaultAssay(marker_object) <- "ADT"
        
        adt_markers <- FindAllMarkers(
          object = marker_object,
          assay = "ADT",
          slot = "data",
          only.pos = TRUE,
          test.use = "wilcox",
          logfc.threshold =
            adt_marker_logfc_threshold,
          min.pct = adt_marker_min_pct,
          max.cells.per.ident =
            max_cells_per_group_for_markers,
          random.seed = 1234,
          densify = FALSE,
          verbose = FALSE
        )
        
        if (nrow(adt_markers) > 0) {
          
          if (!"gene" %in% colnames(adt_markers)) {
            adt_markers$gene <- rownames(
              adt_markers
            )
          }
          
          write_csv_safe(
            adt_markers,
            file.path(
              adt_marker_dir,
              "ADT_all_positive_markers.csv"
            )
          )
          
          top_adt_markers <-
            top_markers_per_group(
              adt_markers,
              adt_top_markers_per_celltype
            )
          
          write_csv_safe(
            top_adt_markers,
            file.path(
              adt_marker_dir,
              "ADT_top_markers_per_celltype.csv"
            )
          )
          
          top_adt_features <- unique(
            top_adt_markers$gene
          )
          
          top_adt_features <- head(
            top_adt_features,
            adt_top_features_for_dotplot
          )
          
          if (length(top_adt_features) > 0) {
            
            DefaultAssay(marker_object) <-
              "ADT"
            
            p_top_adt <- DotPlot(
              object = marker_object,
              features = top_adt_features,
              group.by =
                "celltype_for_analysis",
              assay = "ADT",
              scale = TRUE
            ) +
              RotatedAxis() +
              labs(
                title = paste0(
                  gse_id,
                  " — top differential ADT markers"
                ),
                x = "ADT marker",
                y = "Annotated cell type"
              ) +
              theme(
                plot.title = element_text(
                  face = "bold"
                )
              )
            
            if (isTRUE(export_web_plot_data)) {
              
              dotplot_data_list[[length(dotplot_data_list) + 1L]] <- standardize_dotplot_data(
                plot_object = p_top_adt,
                dataset_id = gse_id,
                plot_id =
                  "ADT_top_differential_markers",
                plot_file = file.path(
                  "02_ADT_markers",
                  "ADT_top_markers_dotplot.png"
                ),
                modality = "ADT"
              )
            }
            
            save_plot_safe(
              p_top_adt,
              file.path(
                adt_marker_dir,
                "ADT_top_markers_dotplot.png"
              ),
              width = max(
                12,
                0.38 *
                  length(top_adt_features) +
                  5
              ),
              height = max(
                7,
                0.35 *
                  length(eligible_groups) +
                  3
              )
            )
          }
        }
        
        rm(marker_object)
        gc()
      }
      
      if (nrow(adt_markers) == 0) {
        
        write_text_safe(
          c(
            "No ADT marker table was produced.",
            paste0(
              "Eligible groups with at least ",
              min_cells_for_marker_group,
              " cells: ",
              length(eligible_groups)
            )
          ),
          file.path(
            adt_marker_dir,
            "ADT_markers_not_available.txt"
          )
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.6 RNA differential genes and GO/KEGG enrichment
      ## ---------------------------------------------------
      
      rna_markers <- data.frame()
      top_rna_markers <- data.frame()
      enrichment_table <- data.frame()
      enrichment_gene_mapping <- data.frame()
      enrichment_input_summary <- data.frame()
      
      run_rna_marker_step <- (
        isTRUE(run_rna_findallmarkers) ||
          isTRUE(run_rna_deg_enrichment)
      )
      
      if (
        run_rna_marker_step &&
        length(eligible_groups) >= 2 &&
        length(marker_cells) > 0
      ) {
        
        rna_marker_object <- subset(
          obj,
          cells = marker_cells
        )
        
        Idents(rna_marker_object) <-
          "celltype_for_analysis"
        
        DefaultAssay(rna_marker_object) <-
          "RNA"
        
        rna_markers <- FindAllMarkers(
          object = rna_marker_object,
          assay = "RNA",
          slot = "data",
          only.pos = TRUE,
          test.use = "wilcox",
          logfc.threshold =
            rna_marker_logfc_threshold,
          min.pct = rna_marker_min_pct,
          max.cells.per.ident =
            max_cells_per_group_for_markers,
          random.seed = 1234,
          densify = FALSE,
          verbose = FALSE
        )
        
        if (nrow(rna_markers) > 0) {
          
          if (!"gene" %in% colnames(rna_markers)) {
            rna_markers$gene <- rownames(
              rna_markers
            )
          }
          
          write_csv_safe(
            rna_markers,
            file.path(
              enrichment_dir,
              "RNA_all_positive_markers.csv"
            )
          )
          
          top_rna_markers <-
            top_markers_per_group(
              rna_markers,
              rna_top_markers_per_celltype
            )
          
          write_csv_safe(
            top_rna_markers,
            file.path(
              enrichment_dir,
              "RNA_top_markers_per_celltype.csv"
            )
          )
          
          top_rna_features <- unique(
            top_rna_markers$gene
          )
          
          top_rna_features <- head(
            top_rna_features,
            rna_top_features_for_dotplot
          )
          
          if (length(top_rna_features) > 0) {
            
            DefaultAssay(rna_marker_object) <-
              "RNA"
            
            p_top_rna <- DotPlot(
              object = rna_marker_object,
              features = top_rna_features,
              group.by =
                "celltype_for_analysis",
              assay = "RNA",
              scale = TRUE
            ) +
              RotatedAxis() +
              labs(
                title = paste0(
                  gse_id,
                  " — top differential RNA markers"
                ),
                x = "RNA marker",
                y = "Annotated cell type"
              ) +
              theme(
                plot.title = element_text(
                  face = "bold"
                )
              )
            
            if (isTRUE(export_web_plot_data)) {
              
              dotplot_data_list[[length(dotplot_data_list) + 1L]] <- standardize_dotplot_data(
                plot_object = p_top_rna,
                dataset_id = gse_id,
                plot_id =
                  "RNA_top_differential_markers",
                plot_file = file.path(
                  "05_RNA_DEG_GO_KEGG",
                  "RNA_top_markers_dotplot.png"
                ),
                modality = "RNA"
              )
            }
            
            save_plot_safe(
              p_top_rna,
              file.path(
                enrichment_dir,
                "RNA_top_markers_dotplot.png"
              ),
              width = max(
                12,
                0.35 *
                  length(top_rna_features) +
                  5
              ),
              height = max(
                7,
                0.35 *
                  length(eligible_groups) +
                  3
              )
            )
          }
          
          if (isTRUE(run_rna_deg_enrichment)) {
            
            enrichment_result <-
              run_celltype_deg_enrichment(
                marker_table =
                  rna_markers,
                object =
                  rna_marker_object,
                dataset_id =
                  gse_id,
                OrgDb =
                  org_db,
                kegg_organism =
                  kegg_organism,
                go_ontologies =
                  go_ontologies,
                p_adjust_threshold =
                  rna_deg_padj_threshold,
                minimum_log2fc =
                  rna_deg_min_log2fc_for_enrichment,
                minimum_pct =
                  rna_deg_min_pct_for_enrichment,
                minimum_genes =
                  min_genes_for_enrichment,
                maximum_genes =
                  max_genes_per_celltype_for_enrichment,
                enrichment_pvalue_cutoff =
                  enrichment_pvalue_cutoff,
                enrichment_qvalue_cutoff =
                  enrichment_qvalue_cutoff
              )
            
            enrichment_table <-
              enrichment_result$enrichment
            
            enrichment_gene_mapping <-
              enrichment_result$gene_mapping
            
            enrichment_input_summary <-
              enrichment_result$input_summary
            
            write_csv_safe(
              enrichment_gene_mapping,
              file.path(
                enrichment_dir,
                "RNA_DEG_gene_to_Entrez_mapping.csv"
              )
            )
            
            write_csv_safe(
              enrichment_input_summary,
              file.path(
                enrichment_dir,
                "RNA_DEG_enrichment_input_summary.csv"
              )
            )
            
            if (nrow(enrichment_table) > 0) {
              
              enrichment_table <- enrichment_table[
                order(
                  enrichment_table$celltype,
                  enrichment_table$database,
                  enrichment_table$p_adjust,
                  na.last = TRUE
                ),
                ,
                drop = FALSE
              ]
              
              write_csv_safe(
                enrichment_table,
                file.path(
                  enrichment_dir,
                  "RNA_DEG_GO_KEGG_enrichment_all.csv"
                )
              )
              
              go_plot_file <- file.path(
                "05_RNA_DEG_GO_KEGG",
                "GO_enrichment_overview.png"
              )
              
              go_plot_data <-
                prepare_enrichment_plot_data(
                  enrichment_table =
                    enrichment_table,
                  database_selector = "GO",
                  number_per_database =
                    enrichment_top_terms_per_database,
                  plot_file =
                    go_plot_file
                )
              
              if (nrow(go_plot_data) > 0) {
                
                p_go_overview <-
                  make_enrichment_overview_plot(
                    plot_data =
                      go_plot_data,
                    dataset_id =
                      gse_id,
                    title_suffix =
                      "GO enrichment of cell-type differential genes",
                    facet_columns = 3L
                  )
                
                save_plot_safe(
                  p_go_overview,
                  file.path(
                    output_dir,
                    go_plot_file
                  ),
                  width = 16,
                  height = max(
                    9,
                    4.2 *
                      ceiling(
                        length(
                          unique(
                            go_plot_data$celltype
                          )
                        ) /
                          3
                      )
                  )
                )
                
                if (
                  isTRUE(
                    export_web_plot_data
                  )
                ) {
                  
                  enrichment_plot_data_list[[length(enrichment_plot_data_list) + 1L]] <- go_plot_data
                }
              }
              
              kegg_plot_file <- file.path(
                "05_RNA_DEG_GO_KEGG",
                "KEGG_enrichment_overview.png"
              )
              
              kegg_plot_data <-
                prepare_enrichment_plot_data(
                  enrichment_table =
                    enrichment_table,
                  database_selector =
                    "KEGG",
                  number_per_database =
                    enrichment_top_terms_per_database,
                  plot_file =
                    kegg_plot_file
                )
              
              if (nrow(kegg_plot_data) > 0) {
                
                p_kegg_overview <-
                  make_enrichment_overview_plot(
                    plot_data =
                      kegg_plot_data,
                    dataset_id =
                      gse_id,
                    title_suffix =
                      "KEGG enrichment of cell-type differential genes",
                    facet_columns = 3L
                  )
                
                save_plot_safe(
                  p_kegg_overview,
                  file.path(
                    output_dir,
                    kegg_plot_file
                  ),
                  width = 16,
                  height = max(
                    9,
                    4.2 *
                      ceiling(
                        length(
                          unique(
                            kegg_plot_data$celltype
                          )
                        ) /
                          3
                      )
                  )
                )
                
                if (
                  isTRUE(
                    export_web_plot_data
                  )
                ) {
                  
                  enrichment_plot_data_list[[length(enrichment_plot_data_list) + 1L]] <- kegg_plot_data
                }
              }
              
              per_celltype_dir <- file.path(
                enrichment_dir,
                "per_celltype"
              )
              
              dir.create(
                per_celltype_dir,
                recursive = TRUE,
                showWarnings = FALSE
              )
              
              enrichment_celltypes <- unique(
                enrichment_table$celltype
              )
              
              for (
                current_celltype in
                enrichment_celltypes
              ) {
                
                celltype_plot_file <-
                  file.path(
                    "05_RNA_DEG_GO_KEGG",
                    "per_celltype",
                    paste0(
                      safe_filename(
                        current_celltype
                      ),
                      "_GO_KEGG_enrichment.png"
                    )
                  )
                
                current_plot_data <-
                  prepare_enrichment_plot_data(
                    enrichment_table =
                      enrichment_table,
                    database_selector =
                      current_celltype,
                    number_per_database =
                      enrichment_top_terms_per_database,
                    plot_file =
                      celltype_plot_file
                  )
                
                if (
                  nrow(current_plot_data) ==
                  0
                ) {
                  next
                }
                
                p_celltype_enrichment <-
                  make_enrichment_overview_plot(
                    plot_data =
                      current_plot_data,
                    dataset_id =
                      gse_id,
                    title_suffix =
                      paste0(
                        current_celltype,
                        " differential-gene GO/KEGG enrichment"
                      ),
                    facet_columns = 1L
                  )
                
                save_plot_safe(
                  p_celltype_enrichment,
                  file.path(
                    output_dir,
                    celltype_plot_file
                  ),
                  width = 11,
                  height = max(
                    7,
                    0.32 *
                      nrow(
                        current_plot_data
                      ) +
                      3
                  )
                )
                
                if (
                  isTRUE(
                    export_web_plot_data
                  )
                ) {
                  
                  enrichment_plot_data_list[[length(enrichment_plot_data_list) + 1L]] <- current_plot_data
                }
              }
              
            } else {
              
              write_text_safe(
                c(
                  "No significant GO or KEGG enrichment term was produced.",
                  "GO uses the installed OrgDb annotation package.",
                  "KEGG may also require an active internet connection."
                ),
                file.path(
                  enrichment_dir,
                  "RNA_DEG_enrichment_not_available.txt"
                )
              )
            }
          }
        }
        
        rm(rna_marker_object)
        gc()
      }
      
      if (
        run_rna_marker_step &&
        nrow(rna_markers) == 0
      ) {
        
        write_text_safe(
          c(
            "No RNA differential-gene table was produced.",
            paste0(
              "Eligible groups with at least ",
              min_cells_for_marker_group,
              " cells: ",
              length(eligible_groups)
            )
          ),
          file.path(
            enrichment_dir,
            "RNA_markers_not_available.txt"
          )
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.7 Matched RNA–protein correlations
      ## ---------------------------------------------------
      
      pair_table <- build_rna_adt_pairs(
        obj,
        species
      )
      
      write_csv_safe(
        pair_table,
        file.path(
          correlation_dir,
          "matched_RNA_ADT_pairs.csv"
        )
      )
      
      correlation_table <- data.frame()
      
      if (nrow(pair_table) > 0) {
        
        correlation_table <-
          compute_rna_adt_correlations(
            object = obj,
            pair_table = pair_table,
            group_column =
              "celltype_for_analysis",
            minimum_cells =
              min_cells_for_correlation,
            maximum_cells_per_group =
              max_cells_per_group_for_correlation,
            include_all_cells =
              include_all_cells_correlation
          )
        
        if (nrow(correlation_table) > 0) {
          
          write_csv_safe(
            correlation_table,
            file.path(
              correlation_dir,
              "RNA_ADT_spearman_correlations.csv"
            )
          )
          
          ranked_pairs <-
            rank_correlation_pairs(
              correlation_table
            )
          
          heatmap_pairs <- head(
            ranked_pairs,
            max_correlation_pairs_in_heatmap
          )
          
          heatmap_correlation_data <-
            correlation_table[
              correlation_table$celltype !=
                "All cells" &
                correlation_table$pair_label %in%
                heatmap_pairs,
              ,
              drop = FALSE
            ]
          
          if (
            nrow(heatmap_correlation_data) >
            0
          ) {
            
            heatmap_correlation_data$pair_label <-
              factor(
                heatmap_correlation_data$pair_label,
                levels = rev(heatmap_pairs)
              )
            
            p_correlation_heatmap <- ggplot(
              heatmap_correlation_data,
              aes(
                x = pair_label,
                y = celltype,
                fill = spearman_rho
              )
            ) +
              geom_tile(
                color = "white",
                linewidth = 0.15
              ) +
              scale_fill_gradient2(
                low = "navy",
                mid = "white",
                high = "firebrick",
                midpoint = 0,
                limits = c(-1, 1),
                na.value = "grey85"
              ) +
              labs(
                title = paste0(
                  gse_id,
                  " — cell-type-specific RNA–protein correlations"
                ),
                x = "Matched RNA–protein pair",
                y = "Annotated cell type",
                fill = "Spearman\nrho"
              ) +
              theme_minimal() +
              theme(
                plot.title = element_text(
                  face = "bold"
                ),
                axis.text.x = element_text(
                  angle = 45,
                  hjust = 1
                )
              )
            
            save_plot_safe(
              p_correlation_heatmap,
              file.path(
                correlation_dir,
                "RNA_ADT_correlation_heatmap.png"
              ),
              width = max(
                12,
                0.46 *
                  length(heatmap_pairs) +
                  5
              ),
              height = max(
                7,
                0.35 *
                  nlevels(
                    obj$celltype_for_analysis
                  ) +
                  3
              )
            )
            
            if (isTRUE(export_web_plot_data)) {
              
              heatmap_export_data <-
                heatmap_correlation_data
              
              heatmap_export_data$pair_label <-
                as.character(
                  heatmap_export_data$pair_label
                )
              
              heatmap_export_data$dataset_id <-
                gse_id
              
              heatmap_export_data$pair_order <-
                match(
                  heatmap_export_data$pair_label,
                  rev(heatmap_pairs)
                )
              
              heatmap_export_data$celltype_order <-
                factor_order_index(
                  heatmap_export_data$celltype
                )
              
              heatmap_export_data$plot_file <-
                normalize_path_for_table(
                  file.path(
                    "03_RNA_ADT_correlation",
                    "RNA_ADT_correlation_heatmap.png"
                  )
                )
              
              write_csv_safe(
                heatmap_export_data,
                file.path(
                  web_data_dir,
                  "RNA_ADT_correlation_heatmap_plot_data.csv"
                )
              )
            }
          }
          
          scatter_pairs <- head(
            ranked_pairs,
            number_of_scatter_pairs
          )
          
          if (length(scatter_pairs) > 0) {
            
            for (
              scatter_index in
              seq_along(scatter_pairs)
            ) {
              
              pair_label_current <-
                scatter_pairs[
                  scatter_index
                ]
              
              pair_row <- pair_table[
                match(
                  pair_label_current,
                  pair_table$pair_label
                ),
                ,
                drop = FALSE
              ]
              
              if (nrow(pair_row) == 0) {
                next
              }
              
              p_scatter <-
                make_pair_scatter_plot(
                  object = obj,
                  pair_row = pair_row,
                  group_column =
                    "celltype_for_analysis",
                  maximum_cells =
                    max_cells_for_scatter
                )
              
              scatter_plot_file <- file.path(
                "03_RNA_ADT_correlation",
                paste0(
                  "scatter_",
                  sprintf(
                    "%02d",
                    scatter_index
                  ),
                  "_",
                  safe_filename(
                    pair_row$marker
                  ),
                  ".png"
                )
              )
              
              if (isTRUE(export_web_plot_data)) {
                
                current_scatter_data <-
                  tryCatch(
                    {
                      as.data.frame(
                        p_scatter$data
                      )
                    },
                    error = function(e) {
                      data.frame()
                    }
                  )
                
                if (
                  nrow(current_scatter_data) >
                  0
                ) {
                  
                  current_scatter_data$dataset_id <-
                    gse_id
                  
                  current_scatter_data$marker <-
                    pair_row$marker
                  
                  current_scatter_data$pair_label <-
                    pair_row$pair_label
                  
                  current_scatter_data$rna_genes <-
                    pair_row$rna_genes
                  
                  current_scatter_data$adt_feature <-
                    pair_row$adt_feature
                  
                  current_scatter_data$plot_file <-
                    normalize_path_for_table(
                      scatter_plot_file
                    )
                  
                  scatter_plot_data_list[[length(scatter_plot_data_list) + 1L]] <- current_scatter_data
                }
              }
              
              save_plot_safe(
                p_scatter,
                file.path(
                  output_dir,
                  scatter_plot_file
                ),
                width = 10,
                height = 7
              )
            }
          }
        }
      }
      
      if (
        nrow(pair_table) == 0 ||
        nrow(correlation_table) == 0
      ) {
        
        write_text_safe(
          c(
            "No RNA–ADT correlation table was produced.",
            paste0(
              "Matched RNA–ADT pairs: ",
              nrow(pair_table)
            ),
            paste0(
              "Minimum cells per cell type: ",
              min_cells_for_correlation
            )
          ),
          file.path(
            correlation_dir,
            "RNA_ADT_correlations_not_available.txt"
          )
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.8 WNN modality weights
      ## ---------------------------------------------------
      
      weight_detection <-
        detect_modality_weight_columns(
          obj
        )
      
      obj <- weight_detection$object
      
      rna_weight_column <-
        weight_detection$rna_column
      
      adt_weight_column <-
        weight_detection$adt_column
      
      weight_long <- data.frame()
      
      weight_summary <-
        summarize_modality_weights(
          object = obj,
          group_column =
            "celltype_for_analysis",
          rna_weight_column =
            rna_weight_column,
          adt_weight_column =
            adt_weight_column
        )
      
      if (nrow(weight_summary) > 0) {
        
        write_csv_safe(
          weight_summary,
          file.path(
            weight_dir,
            "WNN_modality_weight_summary_by_celltype.csv"
          )
        )
        
        weight_long <-
          make_modality_weight_long_data(
            object = obj,
            group_column =
              "celltype_for_analysis",
            rna_weight_column =
              rna_weight_column,
            adt_weight_column =
              adt_weight_column,
            maximum_per_group =
              max_cells_per_group_for_weight_plot
          )
        
        if (nrow(weight_long) > 0) {
          
          if (isTRUE(export_web_plot_data)) {
            
            weight_long$dataset_id <- gse_id
            weight_long$plot_file <-
              normalize_path_for_table(
                file.path(
                  "04_WNN_modality_weights",
                  "WNN_modality_weight_violin.png"
                )
              )
            
            write_csv_safe(
              weight_long,
              file.path(
                web_data_dir,
                "WNN_modality_weight_violin_plot_data.csv"
              )
            )
          }
          
          p_weight_violin <- ggplot(
            weight_long,
            aes(
              x = CellType,
              y = Weight
            )
          ) +
            geom_violin(
              scale = "width",
              trim = TRUE
            ) +
            geom_boxplot(
              width = 0.13,
              outlier.shape = NA
            ) +
            facet_wrap(
              ~ Modality,
              ncol = 1
            ) +
            labs(
              title = paste0(
                gse_id,
                " — WNN modality weights"
              ),
              x = "Annotated cell type",
              y = "Modality weight"
            ) +
            theme_classic() +
            theme(
              plot.title = element_text(
                face = "bold"
              ),
              axis.text.x = element_text(
                angle = 45,
                hjust = 1
              )
            )
          
          save_plot_safe(
            p_weight_violin,
            file.path(
              weight_dir,
              "WNN_modality_weight_violin.png"
            ),
            width = max(
              11,
              0.45 *
                nlevels(
                  obj$celltype_for_analysis
                ) +
                5
            ),
            height = 9
          )
        }
        
        feature_weight_columns <- c(
          rna_weight_column,
          adt_weight_column
        )
        
        feature_weight_columns <-
          feature_weight_columns[
            !is.na(
              feature_weight_columns
            )
          ]
        
        if (
          length(feature_weight_columns) >
          0
        ) {
          
          p_weight_umap <- FeaturePlot(
            object = obj,
            features =
              feature_weight_columns,
            reduction = reduction_name,
            ncol =
              length(feature_weight_columns),
            raster = NULL
          ) &
            theme(
              plot.title = element_text(
                face = "bold"
              )
            )
          
          save_plot_safe(
            p_weight_umap,
            file.path(
              weight_dir,
              "WNN_modality_weights_on_UMAP.png"
            ),
            width = 6 *
              length(feature_weight_columns),
            height = 6
          )
        }
      } else {
        
        write_text_safe(
          c(
            "No WNN modality-weight columns were found.",
            paste0(
              "Available metadata columns: ",
              paste(
                colnames(obj@meta.data),
                collapse = ", "
              )
            )
          ),
          file.path(
            weight_dir,
            "WNN_weights_not_available.txt"
          )
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.9 Web-ready plotting data
      ## ---------------------------------------------------
      
      if (isTRUE(export_web_plot_data)) {
        
        cell_embedding_data <-
          build_cell_embedding_plot_data(
            object = obj,
            dataset_id = gse_id,
            reduction_name =
              reduction_name,
            cluster_column =
              cluster_column,
            annotation_column =
              annotation_column,
            rna_weight_column =
              rna_weight_column,
            adt_weight_column =
              adt_weight_column
          )
        
        write_csv_safe(
          cell_embedding_data,
          file.path(
            web_data_dir,
            "cell_embedding_metadata_plot_data.csv"
          )
        )
        
        combined_dotplot_data <-
          rbind_fill_data_frames(
            dotplot_data_list
          )
        
        if (nrow(combined_dotplot_data) > 0) {
          
          write_csv_safe(
            combined_dotplot_data,
            file.path(
              web_data_dir,
              "marker_dotplot_all_plot_data.csv"
            )
          )
        }
        
        combined_scatter_data <-
          rbind_fill_data_frames(
            scatter_plot_data_list
          )
        
        if (nrow(combined_scatter_data) > 0) {
          
          write_csv_safe(
            combined_scatter_data,
            file.path(
              web_data_dir,
              "RNA_ADT_scatter_all_plot_data.csv"
            )
          )
        }
        
        combined_enrichment_plot_data <-
          rbind_fill_data_frames(
            enrichment_plot_data_list
          )
        
        if (
          nrow(
            combined_enrichment_plot_data
          ) > 0
        ) {
          
          write_csv_safe(
            combined_enrichment_plot_data,
            file.path(
              web_data_dir,
              "RNA_DEG_GO_KEGG_enrichment_plot_data.csv"
            )
          )
        }
        
        plot_data_manifest <-
          build_plot_data_manifest(
            output_dir =
              output_dir,
            web_data_dir =
              web_data_dir
          )
        
        manifest_file <- file.path(
          web_data_dir,
          "plot_data_manifest.csv"
        )
        
        write_csv_safe(
          plot_data_manifest,
          manifest_file
        )
        
        write_plot_data_readme(
          file = file.path(
            web_data_dir,
            "README_plot_data.txt"
          ),
          manifest_file =
            manifest_file
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.10 Optional analysis-ready object
      ## ---------------------------------------------------
      
      if (isTRUE(save_analysis_ready_rds)) {
        
        analysis_ready_file <- file.path(
          output_dir,
          paste0(
            gse_id,
            "_CITE_analysis_ready.rds"
          )
        )
        
        saveRDS(
          obj,
          analysis_ready_file,
          compress = FALSE
        )
      }
      
      
      ## ---------------------------------------------------
      ## 13.11 Save run information and completion flag
      ## ---------------------------------------------------
      
      run_information <- c(
        paste0("GSE: ", gse_id),
        paste0(
          "Input RDS: ",
          input_rds_file
        ),
        paste0(
          "Completed: ",
          timestamp_now()
        ),
        paste0(
          "Cells: ",
          ncol(obj)
        ),
        paste0(
          "RNA features: ",
          nrow(obj[["RNA"]])
        ),
        paste0(
          "ADT features: ",
          nrow(obj[["ADT"]])
        ),
        paste0(
          "Annotation source: ",
          annotation_column
        ),
        paste0(
          "Cluster column: ",
          cluster_column
        ),
        paste0(
          "Reduction: ",
          reduction_name
        ),
        paste0(
          "Matched RNA–ADT pairs: ",
          nrow(pair_table)
        ),
        paste0(
          "RNA weight column: ",
          ifelse(
            is.na(rna_weight_column),
            "not found",
            rna_weight_column
          )
        ),
        paste0(
          "ADT weight column: ",
          ifelse(
            is.na(adt_weight_column),
            "not found",
            adt_weight_column
          )
        ),
        "",
        "Session information:",
        capture.output(sessionInfo())
      )
      
      write_text_safe(
        run_information,
        file.path(
          output_dir,
          "run_information.txt"
        )
      )
      
      write_text_safe(
        paste0(
          "Completed successfully at ",
          timestamp_now()
        ),
        completion_flag
      )
      
      write(
        paste0(
          "[",
          timestamp_now(),
          "] Completed successfully."
        ),
        file = log_file,
        append = TRUE
      )
      
      output_row <- data.frame(
        GSE = gse_id,
        status = "COMPLETED",
        input_rds = input_rds_file,
        output_dir = output_dir,
        cells = ncol(obj),
        RNA_features =
          nrow(obj[["RNA"]]),
        ADT_features =
          nrow(obj[["ADT"]]),
        annotation_column =
          annotation_column,
        cluster_column =
          cluster_column,
        reduction =
          reduction_name,
        matched_RNA_ADT_pairs =
          nrow(pair_table),
        message = "CITE-specific analysis completed.",
        stringsAsFactors = FALSE
      )
      
      rm(obj)
      gc()
      
      message(
        "[DONE] ",
        gse_id,
        "\nOutput: ",
        output_dir
      )
      
      output_row
    },
    
    error = function(e) {
      
      error_message <- conditionMessage(
        e
      )
      
      error_text <- c(
        paste0("GSE: ", gse_id),
        paste0(
          "Input RDS: ",
          input_rds_file
        ),
        paste0(
          "Failed: ",
          timestamp_now()
        ),
        paste0(
          "Error: ",
          error_message
        ),
        "",
        "Traceback:",
        capture.output(
          traceback()
        ),
        "",
        "Session information:",
        capture.output(
          sessionInfo()
        )
      )
      
      write_text_safe(
        error_text,
        file.path(
          output_dir,
          "analysis_error.txt"
        )
      )
      
      write(
        paste0(
          "[",
          timestamp_now(),
          "] Failed: ",
          error_message
        ),
        file = log_file,
        append = TRUE
      )
      
      message(
        "[FAILED] ",
        gse_id,
        ": ",
        error_message
      )
      
      if (exists("obj")) {
        rm(obj)
      }
      
      gc()
      
      data.frame(
        GSE = gse_id,
        status = "FAILED",
        input_rds = input_rds_file,
        output_dir = output_dir,
        cells = NA_integer_,
        RNA_features = NA_integer_,
        ADT_features = NA_integer_,
        annotation_column =
          NA_character_,
        cluster_column =
          NA_character_,
        reduction =
          NA_character_,
        matched_RNA_ADT_pairs =
          NA_integer_,
        message = error_message,
        stringsAsFactors = FALSE
      )
    }
  )
  
  result
}


## =========================================================
## 15. Find all annotated RDS files
## =========================================================

input_rds_files <- list.files(
  path = root_dir,
  pattern = input_rds_pattern,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

## Prevent accidental re-analysis of RDS files saved inside
## an output folder.
input_rds_files <- input_rds_files[
  !grepl(
    paste0(
      "[/\\\\]",
      output_folder_name,
      "[/\\\\]"
    ),
    input_rds_files,
    ignore.case = TRUE
  )
]

input_rds_files <- sort(
  unique(
    normalizePath(
      input_rds_files,
      winslash = "/",
      mustWork = TRUE
    )
  )
)

if (length(input_rds_files) == 0) {
  stop(
    "No annotated RDS files matched this pattern:\n",
    input_rds_pattern,
    "\nUnder:\n",
    root_dir
  )
}

if (
  is.finite(test_first_n) &&
  test_first_n > 0
) {
  
  input_rds_files <- head(
    input_rds_files,
    as.integer(test_first_n)
  )
}

cat(
  "\nAnnotated RDS files selected: ",
  length(input_rds_files),
  "\n",
  sep = ""
)

print(
  input_rds_files
)


## =========================================================
## 16. Batch loop
## Each RDS is processed independently. Errors are recorded
## and the next RDS continues.
## =========================================================

batch_results <- vector(
  mode = "list",
  length = length(input_rds_files)
)

for (
  file_index in
  seq_along(input_rds_files)
) {
  
  cat(
    "\n\n[",
    file_index,
    "/",
    length(input_rds_files),
    "] ",
    basename(
      input_rds_files[file_index]
    ),
    "\n",
    sep = ""
  )
  
  batch_results[[file_index]] <-
    analyze_one_annotated_rds(
      input_rds_files[
        file_index
      ]
    )
  
  gc()
}


## =========================================================
## 17. Save batch summary
## =========================================================

batch_summary <- do.call(
  rbind,
  batch_results
)

rownames(batch_summary) <- NULL

batch_summary_file <- file.path(
  root_dir,
  "CITE_specific_analysis_batch_summary.csv"
)

write_csv_safe(
  batch_summary,
  batch_summary_file
)

cat(
  "\n\n",
  paste(rep("=", 72), collapse = ""),
  "\n",
  "Batch CITE-specific analysis finished.\n",
  "Summary file:\n",
  batch_summary_file,
  "\n\n",
  "Status counts:\n",
  sep = ""
)

print(
  table(
    batch_summary$status,
    useNA = "ifany"
  )
)

cat(
  "\nDetailed summary:\n"
)

print(
  batch_summary,
  row.names = FALSE
)

gc()
