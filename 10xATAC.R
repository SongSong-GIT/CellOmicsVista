# ============================================================
# GSM8289638 RNA + ATAC matrix and GSM8289639 fragments
# Signac + Seurat v5
#
# OFFLINE-SAFE FIXED VERSION
#
# Input files expected in raw_dir:
#   GSM8289638_barcodes.tsv.gz
#   GSM8289638_features.tsv.gz
#   GSM8289638_matrix.mtx.gz
#   GSM8289639_atac_fragments.tsv.gz
#
# Key fixes:
#   1) Reads the prefixed 10x MTX/features/barcodes files directly.
#   2) Splits Gene Expression and Peaks using features.tsv.gz column 3.
#   3) Creates the fragments .tbi index automatically when possible.
#   4) Does NOT call seqlevelsStyle() or GetGRangesFromEnsDb().
#   5) Builds genome coverage plots and tables for the top 5 RNA genes.
# ============================================================

gc()
set.seed(1234)

# ============================================================
# 0. Libraries
# ============================================================

library(
  Seurat,
  lib.loc = "D:/R/R-4.3.1/library/seuratV5"
)

# required_pkgs <- c(
#   "Signac",
#   "GenomicRanges",
#   "GenomeInfoDb",
#   "IRanges",
#   "Matrix",
#   "ggplot2",
#   "patchwork",
#   "dplyr",
#   "tidyr",
#   "tibble",
#   "stringr",
#   "data.table",
#   "Rsamtools",
#   "EnsDb.Hsapiens.v86",
#   "ensembldb",
#   "clusterProfiler",
#   "org.Hs.eg.db",
#   "AnnotationDbi",
#   "ggrepel"
# )
# 
# missing_pkgs <- required_pkgs[
#   !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
# ]
# 
# if (length(missing_pkgs) > 0) {
#   stop(
#     "Missing required packages:\n",
#     paste(missing_pkgs, collapse = ", "),
#     "\n\nInstall CRAN packages with install.packages(), and Bioconductor packages with BiocManager::install()."
#   )
# }

suppressPackageStartupMessages({
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(Rsamtools)
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggrepel)
})

message("Seurat version: ", as.character(packageVersion("Seurat")))
message("Signac version: ", as.character(packageVersion("Signac")))

# ============================================================
# 1. Paths, parameters, and output folders
# ============================================================

# IMPORTANT: change only this line to the folder containing the four files.
raw_dir <- "F:\\老师\\数据库\\新一轮\\scATAC-seq 7月6日\\GSE268321\\GSE268321_RAW"

out_dir <- file.path(
  dirname(raw_dir),
  "GSM8289638_GSM8289639_signac_result_complete"
)
plot_dir <- file.path(out_dir, "01_example_figures")
web_dir <- file.path(out_dir, "02_web_plot_data")
analysis_dir <- file.path(out_dir, "03_analysis_tables")
object_dir <- file.path(out_dir, "04_seurat_object")
checkpoint_dir <- file.path(out_dir, "05_checkpoints")

for (x in c(out_dir, plot_dir, web_dir, analysis_dir, object_dir, checkpoint_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

dataset_id <- "GSM8289638_GSM8289639"

qc_min_nCount_peaks <- 9000
qc_max_nCount_peaks <- 100000
qc_max_blacklist_ratio <- 0.01
qc_max_nucleosome_signal <- 4
qc_min_TSS_enrichment <- 4

n_rna_variable_features <- 2000
n_rna_pcs <- 30
n_lsi_dims <- 30
cluster_resolution <- 0.8

n_top_peaks_web <- 30
n_da_peaks_web <- 10
fragment_hist_max_length <- 1000
fragment_hist_region <- GRanges("chr1", IRanges(1, 2000000))
coverage_bin_width <- 100
coverage_extend_bp <- 1000

# Genome-track genes:
# - Leave coverage_genes_manual empty to automatically use the first five
#   annotated Seurat RNA variable features.
# - Or set, for example, c("MAP2", "TUBB3", "SYP", "GFAP", "MBP").
coverage_top_n_genes <- 5
coverage_genes_manual <- character(0)

da_ident_1_requested <- "0"
da_ident_2_requested <- "1"

go_ontology <- "BP"
enrichment_pvalue_cutoff <- 0.05
enrichment_qvalue_cutoff <- 0.05

# ============================================================
# 2. Sample information and file validation
# ============================================================

sample_info <- data.frame(
  sample_id = "GSM8289638_multiome",
  timepoint = "sample1",
  barcode_file = "GSM8289638_barcodes.tsv.gz",
  feature_file = "GSM8289638_features.tsv.gz",
  matrix_file = "GSM8289638_matrix.mtx.gz",
  fragment_file = "GSM8289639_atac_fragments.tsv.gz",
  stringsAsFactors = FALSE
)

sample_info$barcode_path <- file.path(raw_dir, sample_info$barcode_file)
sample_info$feature_path <- file.path(raw_dir, sample_info$feature_file)
sample_info$matrix_path <- file.path(raw_dir, sample_info$matrix_file)
sample_info$fragment_path <- file.path(raw_dir, sample_info$fragment_file)
sample_info$fragment_index_path <- paste0(sample_info$fragment_path, ".tbi")

print(sample_info)

required_input_paths <- c(
  sample_info$barcode_path,
  sample_info$feature_path,
  sample_info$matrix_path,
  sample_info$fragment_path
)

missing_input_paths <- required_input_paths[!file.exists(required_input_paths)]

if (length(missing_input_paths) > 0) {
  stop(
    "The following input files are missing:\n",
    paste(missing_input_paths, collapse = "\n"),
    "\n\nPlease change raw_dir so that it points to the folder shown in your screenshot."
  )
}

ensure_tabix_index <- function(fragment_path) {
  index_path <- paste0(fragment_path, ".tbi")
  
  if (file.exists(index_path)) {
    message("Fragment index found: ", index_path)
    return(index_path)
  }
  
  message(
    "Fragment index was not found. Attempting to create: ",
    index_path
  )
  
  index_error <- NULL
  
  tryCatch(
    {
      Rsamtools::indexTabix(
        file = fragment_path,
        format = "bed"
      )
    },
    error = function(e) {
      index_error <<- conditionMessage(e)
    }
  )
  
  if (!file.exists(index_path)) {
    stop(
      "Could not create the Tabix index for:\n",
      fragment_path,
      "\n\nThe fragment file must be coordinate-sorted and BGZF-compressed.",
      "\nDownload the matching .tbi file or recompress/index the fragments with bgzip + tabix.",
      if (!is.null(index_error)) paste0("\nRsamtools error: ", index_error) else ""
    )
  }
  
  message("Fragment index created: ", index_path)
  index_path
}

sample_info$fragment_index_path <- vapply(
  sample_info$fragment_path,
  ensure_tabix_index,
  character(1)
)

extract_gsm <- function(x) {
  stringr::str_extract(basename(x), "GSM[0-9]+")
}

matrix_gsm <- extract_gsm(sample_info$matrix_file)
fragment_gsm <- extract_gsm(sample_info$fragment_file)

message(
  "Configured matrix files from ", matrix_gsm,
  " and ATAC fragments from ", fragment_gsm, "."
)

# ============================================================
# 3. Helper functions
# ============================================================

normalize_barcode <- function(x) {
  x <- as.character(x)
  x <- sub("^.*_", "", x)
  x <- sub("-1$", "", x)
  x
}

normalize_seqname_offline <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x == "MT"] <- "M"
  paste0("chr", x)
}

standardize_granges_offline <- function(
    gr,
    keep_standard = TRUE,
    clear_genome = TRUE
) {
  stopifnot(inherits(gr, "GRanges"))
  
  standard_chr <- c(
    paste0("chr", 1:22),
    "chrX",
    "chrY",
    "chrM"
  )
  
  new_seqnames <- normalize_seqname_offline(
    as.character(seqnames(gr))
  )
  
  keep <- rep(TRUE, length(gr))
  if (keep_standard) {
    keep <- new_seqnames %in% standard_chr
  }
  
  gr <- gr[keep]
  new_seqnames <- new_seqnames[keep]
  
  gr_out <- GRanges(
    seqnames = new_seqnames,
    ranges = ranges(gr),
    strand = strand(gr)
  )
  
  mcols(gr_out) <- mcols(gr)
  names(gr_out) <- names(gr)
  
  desired_levels <- standard_chr[
    standard_chr %in% unique(new_seqnames)
  ]
  
  if (length(desired_levels) > 0) {
    gr_out <- keepSeqlevels(
      gr_out,
      desired_levels,
      pruning.mode = "coarse"
    )
  }
  
  if (clear_genome) {
    suppressWarnings(genome(gr_out) <- NA_character_)
  }
  
  gr_out
}

load_local_hg38_annotation <- function() {
  message("Loading hg38 annotation locally from EnsDb.Hsapiens.v86...")
  
  ann <- ensembldb::genes(
    EnsDb.Hsapiens.v86,
    return.type = "GRanges"
  )
  
  ann <- standardize_granges_offline(
    ann,
    keep_standard = TRUE,
    clear_genome = TRUE
  )
  
  if (!"gene_id" %in% colnames(mcols(ann))) {
    stop("The local EnsDb annotation does not contain gene_id.")
  }
  
  if (!"gene_name" %in% colnames(mcols(ann))) {
    ann$gene_name <- ann$gene_id
  }
  
  ann$gene_name[
    is.na(ann$gene_name) |
      !nzchar(as.character(ann$gene_name))
  ] <- ann$gene_id[
    is.na(ann$gene_name) |
      !nzchar(as.character(ann$gene_name))
  ]
  
  if (!"gene_biotype" %in% colnames(mcols(ann))) {
    ann$gene_biotype <- NA_character_
  }
  
  ann$type <- "gene"
  
  ann <- sort(ann)
  
  message(
    "Local annotation loaded: ",
    length(ann),
    " genes on ",
    paste(seqlevels(ann), collapse = ", "),
    "."
  )
  
  ann
}

save_plot <- function(plot_object, filename, width, height, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

theme_web <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

safe_cor_test <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  if (
    length(x) < 3 ||
    !is.finite(stats::sd(x)) ||
    !is.finite(stats::sd(y)) ||
    stats::sd(x) == 0 ||
    stats::sd(y) == 0
  ) {
    return(
      data.frame(
        method = tools::toTitleCase(method),
        correlation = NA_real_,
        p_value = NA_real_,
        n = length(x),
        stringsAsFactors = FALSE
      )
    )
  }
  
  result <- suppressWarnings(
    stats::cor.test(
      x = x,
      y = y,
      method = method,
      exact = FALSE
    )
  )
  
  data.frame(
    method = tools::toTitleCase(method),
    correlation = unname(result$estimate),
    p_value = result$p.value,
    n = length(x),
    stringsAsFactors = FALSE
  )
}

ratio_to_numeric <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, "/", fixed = TRUE)
  vapply(
    parts,
    function(z) {
      if (length(z) != 2) return(NA_real_)
      numerator <- suppressWarnings(as.numeric(z[1]))
      denominator <- suppressWarnings(as.numeric(z[2]))
      if (!is.finite(denominator) || denominator == 0) return(NA_real_)
      numerator / denominator
    },
    numeric(1)
  )
}

scan_fragment_region <- function(fragment_path, region_gr) {
  if (!file.exists(fragment_path)) {
    stop("Fragment file does not exist: ", fragment_path)
  }
  if (!file.exists(paste0(fragment_path, ".tbi"))) {
    stop("Tabix index does not exist: ", paste0(fragment_path, ".tbi"))
  }
  
  tabix_file <- Rsamtools::TabixFile(fragment_path)
  lines <- Rsamtools::scanTabix(
    tabix_file,
    param = region_gr
  )[[1]]
  
  if (length(lines) == 0) {
    return(
      data.table(
        chromosome = character(),
        start = integer(),
        end = integer(),
        barcode = character(),
        count = numeric(),
        fragment_length = integer(),
        barcode_clean = character()
      )
    )
  }
  
  fields <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
  
  if (length(fields) < 4) {
    stop("Fragment file has fewer than four tab-separated columns.")
  }
  
  count_values <- if (length(fields) >= 5) {
    suppressWarnings(as.numeric(fields[[5]]))
  } else {
    rep(1, length(fields[[1]]))
  }
  
  count_values[!is.finite(count_values) | count_values <= 0] <- 1
  
  fragment_dt <- data.table(
    chromosome = fields[[1]],
    start = as.integer(fields[[2]]) + 1L,
    end = as.integer(fields[[3]]),
    barcode = fields[[4]],
    count = count_values
  )
  
  fragment_dt[, fragment_length := end - start + 1L]
  fragment_dt[, barcode_clean := normalize_barcode(barcode)]
  
  fragment_dt
}

make_umap_plot <- function(
    data,
    x_col,
    y_col,
    color_col,
    label_col,
    title,
    x_label,
    y_label
) {
  plot_data <- data %>%
    filter(
      is.finite(.data[[x_col]]),
      is.finite(.data[[y_col]]),
      !is.na(.data[[color_col]])
    )
  
  n_points <- nrow(plot_data)
  point_size <- max(0.05, min(0.6, 2500 / max(n_points, 1)))
  
  label_data <- plot_data %>%
    group_by(.data[[label_col]]) %>%
    summarise(
      label_x = median(.data[[x_col]], na.rm = TRUE),
      label_y = median(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(
    plot_data,
    aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[color_col]]
    )
  ) +
    geom_point(size = point_size, alpha = 0.75) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(
        x = label_x,
        y = label_y,
        label = .data[[label_col]]
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3,
      fontface = "bold",
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      color = color_col
    ) +
    theme_web() +
    theme(legend.position = "right")
}

make_enrichment_barplot <- function(enrichment_table, type_value, title) {
  plot_data <- enrichment_table %>%
    filter(
      .data$enrichment_type == type_value,
      is.finite(p.adjust)
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(unique(Description))
      ),
      neg_log10_padj = -log10(
        pmax(p.adjust, .Machine$double.xmin)
      )
    )
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  ggplot(
    plot_data,
    aes(
      x = neg_log10_padj,
      y = Description,
      size = Count
    )
  ) +
    geom_point(alpha = 0.85) +
    labs(
      title = title,
      x = "-log10(adjusted P value)",
      y = NULL,
      size = "Gene count"
    ) +
    theme_web()
}

# ============================================================
# 4. Read one paired multiome sample from MTX/features/barcodes
# ============================================================

read_feature_metadata <- function(feature_path) {
  feature_meta <- utils::read.delim(
    gzfile(feature_path),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  source_ncol <- ncol(feature_meta)
  
  if (source_ncol < 1) {
    stop("The features file is empty: ", feature_path)
  }
  
  feature_meta$feature_id <- as.character(feature_meta[[1]])
  
  if (source_ncol >= 2) {
    feature_meta$feature_name <- as.character(feature_meta[[2]])
  } else {
    feature_meta$feature_name <- feature_meta$feature_id
  }
  
  if (source_ncol >= 3) {
    feature_meta$feature_type <- as.character(feature_meta[[3]])
  } else {
    feature_meta$feature_type <- NA_character_
  }
  
  attr(feature_meta, "source_ncol") <- source_ncol
  feature_meta
}

read_one_atac_sample <- function(
    sample_id,
    timepoint,
    matrix_path,
    barcode_path,
    feature_path,
    fragment_path
) {
  message("========================================")
  message("Reading multiome sample: ", sample_id)
  
  feature_meta <- read_feature_metadata(feature_path)
  
  feature_column <- if (
    attr(feature_meta, "source_ncol") >= 2 &&
    mean(
      !is.na(feature_meta$feature_name) &
      nzchar(trimws(feature_meta$feature_name))
    ) >= 0.95
  ) {
    2L
  } else {
    1L
  }
  
  message(
    "Reading 10x MTX matrix with feature column ",
    feature_column,
    "..."
  )
  
  counts_raw <- ReadMtx(
    mtx = matrix_path,
    cells = barcode_path,
    features = feature_path,
    cell.column = 1,
    feature.column = feature_column,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  
  if (nrow(counts_raw) != nrow(feature_meta)) {
    stop(
      "The matrix row count (", nrow(counts_raw),
      ") does not match the features row count (",
      nrow(feature_meta), ")."
    )
  }
  
  selected_feature_names <- if (feature_column == 2L) {
    feature_meta$feature_name
  } else {
    feature_meta$feature_id
  }
  
  feature_type_lower <- tolower(
    trimws(as.character(feature_meta$feature_type))
  )
  
  coordinate_like <- grepl(
    "^(chr)?([0-9]+|x|y|m|mt):[0-9]+-[0-9]+$",
    selected_feature_names,
    ignore.case = TRUE
  )
  
  rna_index <- which(
    grepl(
      "gene expression|gene_expression|^rna$|^gex$",
      feature_type_lower
    )
  )
  
  atac_index <- which(
    grepl(
      "peak|atac|chromatin accessibility",
      feature_type_lower
    )
  )
  
  if (length(atac_index) == 0) {
    atac_index <- which(coordinate_like)
  }
  
  if (length(rna_index) == 0) {
    rna_index <- which(!coordinate_like)
  }
  
  if (length(rna_index) == 0 || length(atac_index) == 0) {
    feature_type_preview <- paste(
      head(unique(feature_meta$feature_type), 20),
      collapse = ", "
    )
    
    stop(
      "Could not split Gene Expression and Peaks from features.tsv.gz.",
      "\nDetected feature types: ", feature_type_preview,
      "\nThe third column should contain values such as 'Gene Expression' and 'Peaks'."
    )
  }
  
  if (length(intersect(rna_index, atac_index)) > 0) {
    stop("Some features were classified as both RNA and ATAC.")
  }
  
  rna_counts <- counts_raw[rna_index, , drop = FALSE]
  atac_counts <- counts_raw[atac_index, , drop = FALSE]
  
  message(
    "RNA matrix: ",
    nrow(rna_counts), " genes x ",
    ncol(rna_counts), " cells"
  )
  
  message(
    "Raw ATAC matrix: ",
    nrow(atac_counts), " peaks x ",
    ncol(atac_counts), " cells"
  )
  
  message("Standardizing peak chromosome names offline...")
  
  peak_names <- rownames(atac_counts)
  peak_chr <- sub(":.*$", "", peak_names)
  peak_coord <- sub("^[^:]+:", "", peak_names)
  peak_start <- suppressWarnings(
    as.integer(sub("-.*$", "", peak_coord))
  )
  peak_end <- suppressWarnings(
    as.integer(sub("^.*-", "", peak_coord))
  )
  
  peak_chr <- normalize_seqname_offline(peak_chr)
  
  standard_chr <- c(
    paste0("chr", 1:22),
    "chrX",
    "chrY",
    "chrM"
  )
  
  keep_peaks <- (
    peak_chr %in% standard_chr &
      is.finite(peak_start) &
      is.finite(peak_end) &
      peak_end >= peak_start
  )
  
  if (sum(keep_peaks) == 0) {
    stop(
      "No valid peak coordinates were found after parsing row names.",
      "\nExpected peak names such as chr1:10000-10500."
    )
  }
  
  atac_counts <- atac_counts[keep_peaks, , drop = FALSE]
  peak_chr <- peak_chr[keep_peaks]
  peak_start <- peak_start[keep_peaks]
  peak_end <- peak_end[keep_peaks]
  
  standardized_peak_names <- paste0(
    peak_chr, ":", peak_start, "-", peak_end
  )
  
  keep_unique <- !duplicated(standardized_peak_names)
  
  if (sum(!keep_unique) > 0) {
    warning(
      sum(!keep_unique),
      " duplicated peak coordinates were removed."
    )
  }
  
  atac_counts <- atac_counts[keep_unique, , drop = FALSE]
  rownames(atac_counts) <- standardized_peak_names[keep_unique]
  
  message(
    "ATAC matrix after standard chromosome filtering: ",
    nrow(atac_counts), " peaks x ",
    ncol(atac_counts), " cells"
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fragment_path,
    min.cells = 10,
    min.features = 200,
    validate.fragments = TRUE
  )
  
  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    project = sample_id
  )
  
  common_barcodes <- intersect(
    colnames(obj),
    colnames(rna_counts)
  )
  
  if (length(common_barcodes) == 0) {
    stop(
      "No common barcodes were found between the ATAC and RNA matrices."
    )
  }
  
  obj <- subset(obj, cells = common_barcodes)
  rna_counts <- rna_counts[, common_barcodes, drop = FALSE]
  
  obj[["RNA"]] <- SeuratObject::CreateAssay5Object(
    counts = rna_counts
  )
  
  obj$dataset_id <- dataset_id
  obj$sample_id <- sample_id
  obj$timepoint <- timepoint
  obj$orig.ident <- sample_id
  
  obj
}

# ============================================================
# 5. Load object
# ============================================================

gc()

pbmc <- read_one_atac_sample(
  sample_id = sample_info$sample_id[1],
  timepoint = sample_info$timepoint[1],
  matrix_path = sample_info$matrix_path[1],
  barcode_path = sample_info$barcode_path[1],
  feature_path = sample_info$feature_path[1],
  fragment_path = sample_info$fragment_path[1]
)

print(pbmc)

saveRDS(
  pbmc,
  file.path(
    checkpoint_dir,
    "01_after_loading.rds"
  )
)

# ============================================================
# 6. Load local hg38 annotation without any online request
# ============================================================

DefaultAssay(pbmc) <- "peaks"

annotations <- load_local_hg38_annotation()

Annotation(pbmc) <- annotations

message("Annotation successfully assigned without seqlevelsStyle().")

saveRDS(
  pbmc,
  file.path(
    checkpoint_dir,
    "02_after_annotation.rds"
  )
)

# ============================================================
# 7. QC metrics and pre-filter metadata
# ============================================================

message("Computing nucleosome signal...")
pbmc <- NucleosomeSignal(object = pbmc)

message("Computing TSS enrichment...")
pbmc <- TSSEnrichment(object = pbmc)

peak_counts_all <- GetAssayData(
  pbmc,
  assay = "peaks",
  layer = "counts"
)

pbmc$peak_region_fragments <- Matrix::colSums(
  peak_counts_all
)

message("Counting fragments...")
frag_counts <- CountFragments(
  fragments = sample_info$fragment_path[1]
)

frag_counts$barcode_clean <- normalize_barcode(
  frag_counts$CB
)

object_barcode_clean <- normalize_barcode(
  colnames(pbmc)
)

frag_match <- match(
  object_barcode_clean,
  frag_counts$barcode_clean
)

pbmc$passed_filters <- frag_counts$frequency_count[frag_match]

if ("mononucleosome" %in% names(frag_counts)) {
  pbmc$mononucleosome_fragments <- frag_counts$mononucleosome[frag_match]
}
if ("nucleosome_free" %in% names(frag_counts)) {
  pbmc$nucleosome_free_fragments <- frag_counts$nucleosome_free[frag_match]
}
if ("reads_count" %in% names(frag_counts)) {
  pbmc$fragment_reads_count <- frag_counts$reads_count[frag_match]
}

missing_fragment_metric <- is.na(pbmc$passed_filters)

if (any(missing_fragment_metric)) {
  warning(
    sum(missing_fragment_metric),
    " cells could not be matched to CountFragments output. ",
    "passed_filters falls back to peak_region_fragments."
  )
  
  pbmc$passed_filters[missing_fragment_metric] <-
    pbmc$peak_region_fragments[missing_fragment_metric]
}

pbmc$passed_filters[
  !is.finite(pbmc$passed_filters) |
    pbmc$passed_filters <= 0
] <- 1

pbmc$pct_reads_in_peaks <-
  pbmc$peak_region_fragments /
  pbmc$passed_filters *
  100

# -----------------------------
# Blacklist ratio, standardized offline
# -----------------------------

blacklist_env <- new.env(parent = emptyenv())

suppressWarnings(
  try(
    data(
      "blacklist_hg38_unified",
      package = "Signac",
      envir = blacklist_env
    ),
    silent = TRUE
  )
)

if (exists("blacklist_hg38_unified", envir = blacklist_env)) {
  blacklist_regions <- get(
    "blacklist_hg38_unified",
    envir = blacklist_env
  )
  
  blacklist_regions <- standardize_granges_offline(
    blacklist_regions,
    keep_standard = TRUE,
    clear_genome = TRUE
  )
  
  pbmc$blacklist_ratio <- FractionCountsInRegion(
    object = pbmc,
    assay = "peaks",
    regions = blacklist_regions
  )
} else {
  warning(
    "blacklist_hg38_unified could not be loaded. ",
    "blacklist_ratio is set to 0."
  )
  pbmc$blacklist_ratio <- 0
}

pbmc$pass_qc <- with(
  pbmc@meta.data,
  nCount_peaks > qc_min_nCount_peaks &
    nCount_peaks < qc_max_nCount_peaks &
    blacklist_ratio < qc_max_blacklist_ratio &
    nucleosome_signal < qc_max_nucleosome_signal &
    TSS.enrichment > qc_min_TSS_enrichment
)

qc_meta_all <- pbmc@meta.data %>%
  tibble::rownames_to_column("barcode")

n_cells_before_qc <- ncol(pbmc)
n_cells_after_qc <- sum(pbmc$pass_qc, na.rm = TRUE)

message(
  "QC retained ",
  n_cells_after_qc,
  " / ",
  n_cells_before_qc,
  " cells."
)

if (n_cells_after_qc < 10) {
  stop(
    "Fewer than 10 cells passed QC. ",
    "Please inspect QC metrics and thresholds."
  )
}

pbmc <- subset(
  pbmc,
  cells = colnames(pbmc)[pbmc$pass_qc]
)

saveRDS(
  pbmc,
  file.path(
    checkpoint_dir,
    "03_after_qc.rds"
  )
)

gc()

# ============================================================
# 8. RNA and ATAC dimension reduction + WNN
# ============================================================

DefaultAssay(pbmc) <- "RNA"

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = n_rna_variable_features,
  verbose = FALSE
)

pbmc <- ScaleData(
  pbmc,
  features = VariableFeatures(pbmc),
  verbose = FALSE
)

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(pbmc),
  npcs = n_rna_pcs,
  verbose = FALSE
)

rna_npcs_available <- ncol(
  Embeddings(pbmc, reduction = "pca")
)

rna_dims <- seq_len(
  min(n_rna_pcs, rna_npcs_available)
)

pbmc <- RunUMAP(
  pbmc,
  reduction = "pca",
  dims = rna_dims,
  reduction.name = "umap.rna",
  reduction.key = "rnaUMAP_",
  verbose = FALSE
)

DefaultAssay(pbmc) <- "peaks"

pbmc <- RunTFIDF(pbmc, verbose = FALSE)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(
  pbmc,
  n = n_lsi_dims,
  verbose = FALSE
)

lsi_dims_available <- ncol(
  Embeddings(pbmc, reduction = "lsi")
)

if (lsi_dims_available < 2) {
  stop("LSI returned fewer than two dimensions.")
}

atac_dims <- 2:min(
  n_lsi_dims,
  lsi_dims_available
)

pbmc <- RunUMAP(
  pbmc,
  reduction = "lsi",
  dims = atac_dims,
  reduction.name = "umap.atac",
  reduction.key = "atacUMAP_",
  verbose = FALSE
)

pbmc <- FindMultiModalNeighbors(
  pbmc,
  reduction.list = list("pca", "lsi"),
  dims.list = list(rna_dims, atac_dims),
  modality.weight.name = "RNA.weight",
  verbose = FALSE
)

pbmc <- RunUMAP(
  pbmc,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = FALSE
)

pbmc <- FindClusters(
  pbmc,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = cluster_resolution,
  verbose = FALSE
)

pbmc$seurat_clusters <- as.character(Idents(pbmc))

saveRDS(
  pbmc,
  file.path(
    checkpoint_dir,
    "04_after_wnn.rds"
  )
)

# ============================================================
# 9. ATAC gene activity
# ============================================================

DefaultAssay(pbmc) <- "peaks"

message("Computing gene activity matrix...")

gene_activities <- GeneActivity(
  object = pbmc,
  assay = "peaks",
  extend.upstream = 2000,
  extend.downstream = 0,
  gene.id = FALSE,
  verbose = TRUE
)

pbmc[["ACTIVITY"]] <- SeuratObject::CreateAssay5Object(
  counts = gene_activities
)

activity_library_size <- Matrix::colSums(
  GetAssayData(
    pbmc,
    assay = "ACTIVITY",
    layer = "counts"
  )
)

pbmc$nCount_ACTIVITY <- activity_library_size

activity_scale_factor <- stats::median(
  activity_library_size[
    is.finite(activity_library_size) &
      activity_library_size > 0
  ]
)

if (
  length(activity_scale_factor) == 0 ||
  !is.finite(activity_scale_factor) ||
  activity_scale_factor <= 0
) {
  activity_scale_factor <- 10000
}

pbmc <- NormalizeData(
  pbmc,
  assay = "ACTIVITY",
  normalization.method = "LogNormalize",
  scale.factor = activity_scale_factor,
  verbose = FALSE
)

# ============================================================
# 10. RNA and ATAC state scores
# ============================================================

state_genes_list <- list(
  Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "CDK1"),
  Stemness = c("SOX2", "NANOG", "POU5F1", "PROM1", "ALDH1A1"),
  Immune_Activity = c("CD2", "CD3D", "CD3E", "CD4", "CD8A", "GZMB", "PRF1"),
  Inflammation = c("IL6", "TNF", "IL1B", "CXCL8", "CCL2"),
  Hypoxia = c("HIF1A", "VEGFA", "LDHA", "ENO1", "PDK1")
)

calculate_mean_gene_score <- function(
    object,
    assay,
    genes,
    score_name
) {
  assay_data <- GetAssayData(
    object,
    assay = assay,
    layer = "data"
  )
  
  genes_found <- intersect(
    genes,
    rownames(assay_data)
  )
  
  if (length(genes_found) == 0) {
    warning(
      "No genes were found for ",
      score_name,
      " in assay ",
      assay,
      ". The score is set to 0."
    )
    return(rep(0, ncol(object)))
  }
  
  as.numeric(
    Matrix::colMeans(
      assay_data[
        genes_found,
        ,
        drop = FALSE
      ]
    )
  )
}

for (state_name in names(state_genes_list)) {
  rna_col <- paste0("RNA_Score_", state_name)
  atac_col <- paste0("ATAC_Score_", state_name)
  
  pbmc[[rna_col]] <- calculate_mean_gene_score(
    object = pbmc,
    assay = "RNA",
    genes = state_genes_list[[state_name]],
    score_name = state_name
  )
  
  pbmc[[atac_col]] <- calculate_mean_gene_score(
    object = pbmc,
    assay = "ACTIVITY",
    genes = state_genes_list[[state_name]],
    score_name = state_name
  )
}

# ============================================================
# 11. Candidate marker scoring and automatic annotation
# ============================================================

candidate_markers <- list(
  Keratinocyte = c("KRT5", "KRT14", "KRT1", "KRT10"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
  Endothelial = c("PECAM1", "VWF", "CDH5"),
  T_Cell = c("CD3D", "CD3E", "CD3G", "TRAC"),
  B_Cell = c("CD19", "MS4A1", "CD79A"),
  Myeloid = c("LYZ", "CD14", "AIF1"),
  Melanocyte = c("PMEL", "TYR", "MLANA"),
  Mast_Cell = c("TPSAB1", "TPSB2", "CPA3"),
  Pericyte = c("ACTA2", "PDGFRB"),
  NK_Cell = c("NCAM1", "NCR1", "KLRD1"),
  Plasma_Cell = c("MZB1", "IGKC", "IGHG1"),
  Dendritic_Cell = c("CD1C", "CLEC9A", "LILRA4"),
  Neuron = c("MAP2", "TUBB3", "SYP"),
  Astrocyte = c("GFAP", "S100B"),
  Oligodendrocyte = c("MBP", "PLP1", "SOX10"),
  Smooth_Muscle = c("MYH11", "CNN1"),
  Epithelial = c("EPCAM", "KRT8", "KRT18")
)

DefaultAssay(pbmc) <- "RNA"

rna_data <- GetAssayData(
  pbmc,
  assay = "RNA",
  layer = "data"
)

marker_score_matrix <- matrix(
  0,
  nrow = ncol(pbmc),
  ncol = length(candidate_markers),
  dimnames = list(
    colnames(pbmc),
    names(candidate_markers)
  )
)

for (celltype_name in names(candidate_markers)) {
  genes_found <- intersect(
    candidate_markers[[celltype_name]],
    rownames(rna_data)
  )
  
  if (length(genes_found) > 0) {
    marker_score_matrix[, celltype_name] <-
      as.numeric(
        Matrix::colMeans(
          rna_data[
            genes_found,
            ,
            drop = FALSE
          ]
        )
      )
  }
}

marker_score_meta_names <- paste0(
  "MarkerScore_",
  colnames(marker_score_matrix)
)

stopifnot(length(marker_score_meta_names) == ncol(marker_score_matrix))
stopifnot(nrow(marker_score_matrix) == ncol(pbmc))
stopifnot(all(colnames(pbmc) %in% rownames(marker_score_matrix)))

marker_score_matrix <- marker_score_matrix[
  colnames(pbmc),
  ,
  drop = FALSE
]

for (j in seq_along(marker_score_meta_names)) {
  pbmc[[marker_score_meta_names[j]]] <-
    as.numeric(marker_score_matrix[, j])
}

marker_score_df <- as.data.frame(
  marker_score_matrix,
  check.names = FALSE
)

marker_score_df$seurat_clusters <-
  as.character(pbmc$seurat_clusters)

cluster_marker_means <- marker_score_df %>%
  group_by(seurat_clusters) %>%
  summarise(
    across(
      all_of(names(candidate_markers)),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

annotation_rows <- lapply(
  seq_len(nrow(cluster_marker_means)),
  function(i) {
    score_vector <- unlist(
      cluster_marker_means[
        i,
        names(candidate_markers),
        drop = TRUE
      ],
      use.names = TRUE
    )
    
    score_vector <- as.numeric(score_vector)
    names(score_vector) <- names(candidate_markers)
    
    ordered_types <- names(
      sort(score_vector, decreasing = TRUE)
    )
    
    top1 <- ordered_types[1]
    top2 <- if (length(ordered_types) >= 2) ordered_types[2] else NA_character_
    top3 <- if (length(ordered_types) >= 3) ordered_types[3] else NA_character_
    
    top1_score <- score_vector[top1]
    top2_score <- if (!is.na(top2)) score_vector[top2] else NA_real_
    
    assigned_type <- if (
      !is.finite(top1_score) ||
      top1_score <= 0
    ) {
      "Unassigned"
    } else {
      top1
    }
    
    data.frame(
      seurat_clusters = as.character(
        cluster_marker_means$seurat_clusters[i]
      ),
      cell_type = assigned_type,
      top2_cell_type = top2,
      top3_cell_type = top3,
      top1_score = unname(top1_score),
      top2_score = unname(top2_score),
      score_margin = unname(top1_score - top2_score),
      stringsAsFactors = FALSE
    )
  }
)

auto_annotation <- dplyr::bind_rows(annotation_rows)

cluster_annotation_scores <- cluster_marker_means

names(cluster_annotation_scores)[
  names(cluster_annotation_scores) %in% names(candidate_markers)
] <- paste0(
  "mean_marker_score_",
  names(cluster_annotation_scores)[
    names(cluster_annotation_scores) %in% names(candidate_markers)
  ]
)

cluster_annotation_scores <- auto_annotation %>%
  left_join(
    cluster_annotation_scores,
    by = "seurat_clusters"
  )

cluster_to_celltype <- setNames(
  auto_annotation$cell_type,
  auto_annotation$seurat_clusters
)

pbmc$cell_type <- unname(
  cluster_to_celltype[
    as.character(pbmc$seurat_clusters)
  ]
)

pbmc$cell_type[is.na(pbmc$cell_type)] <- "Unassigned"
Idents(pbmc) <- "cell_type"

# ============================================================
# 12. Unified per-cell web table
# ============================================================

rna_embedding <- Embeddings(pbmc, reduction = "umap.rna")
atac_embedding <- Embeddings(pbmc, reduction = "umap.atac")
wnn_embedding <- Embeddings(pbmc, reduction = "wnn.umap")

final_meta <- pbmc@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(
    seurat_clusters = as.character(seurat_clusters),
    cell_type = as.character(cell_type)
  )

final_meta$rna_umap_1 <- rna_embedding[final_meta$barcode, 1]
final_meta$rna_umap_2 <- rna_embedding[final_meta$barcode, 2]
final_meta$atac_umap_1 <- atac_embedding[final_meta$barcode, 1]
final_meta$atac_umap_2 <- atac_embedding[final_meta$barcode, 2]
final_meta$wnn_umap_1 <- wnn_embedding[final_meta$barcode, 1]
final_meta$wnn_umap_2 <- wnn_embedding[final_meta$barcode, 2]

final_additional_columns <- setdiff(
  names(final_meta),
  names(qc_meta_all)
)

cell_metadata_web <- qc_meta_all %>%
  left_join(
    final_meta %>%
      select(
        barcode,
        all_of(final_additional_columns)
      ),
    by = "barcode"
  )

preferred_cell_columns <- c(
  "barcode",
  "dataset_id",
  "sample_id",
  "timepoint",
  "orig.ident",
  "pass_qc",
  "nCount_RNA",
  "nFeature_RNA",
  "nCount_peaks",
  "nFeature_peaks",
  "peak_region_fragments",
  "passed_filters",
  "pct_reads_in_peaks",
  "TSS.enrichment",
  "nucleosome_signal",
  "blacklist_ratio",
  "mononucleosome_fragments",
  "nucleosome_free_fragments",
  "fragment_reads_count",
  "nCount_ACTIVITY",
  "seurat_clusters",
  "cell_type",
  "rna_umap_1",
  "rna_umap_2",
  "atac_umap_1",
  "atac_umap_2",
  "wnn_umap_1",
  "wnn_umap_2"
)

preferred_cell_columns <- intersect(
  preferred_cell_columns,
  names(cell_metadata_web)
)

other_cell_columns <- setdiff(
  names(cell_metadata_web),
  preferred_cell_columns
)

cell_metadata_web <- cell_metadata_web %>%
  select(
    all_of(preferred_cell_columns),
    all_of(other_cell_columns)
  )

data.table::fwrite(
  cell_metadata_web,
  file.path(web_dir, "01_cell_metadata_web.csv")
)

data.table::fwrite(
  cluster_annotation_scores,
  file.path(web_dir, "02_cluster_annotation_scores_web.csv")
)

# ============================================================
# 13. Unified RNA + ATAC gene table
# ============================================================

select_genes <- c(
  "MS4A1", "CD3D", "CD3E", "CD4", "CD8A", "LYZ", "EPCAM", "MKI67",
  "COL1A1", "PECAM1", "MAP2", "TUBB3", "PMEL", "TPSAB1", "ACTA2",
  "GFAP", "MBP", "IL6", "VEGFA", "SOX2", "TOP2A", "PCNA", "NANOG",
  "PROM1", "GZMB", "PRF1", "TNF", "IL1B", "HIF1A", "LDHA"
)

web_genes_requested <- unique(
  c(
    select_genes,
    unlist(state_genes_list, use.names = FALSE),
    unlist(candidate_markers, use.names = FALSE)
  )
)

rna_data <- GetAssayData(
  pbmc,
  assay = "RNA",
  layer = "data"
)

activity_data <- GetAssayData(
  pbmc,
  assay = "ACTIVITY",
  layer = "data"
)

web_genes_found <- web_genes_requested[
  web_genes_requested %in%
    union(
      rownames(rna_data),
      rownames(activity_data)
    )
]

if (length(web_genes_found) == 0) {
  stop(
    "None of the requested web genes were found. ",
    "Check whether RNA row names are gene symbols or Ensembl IDs."
  )
}

gene_table_list <- lapply(
  web_genes_found,
  function(gene_name) {
    rna_value <- if (gene_name %in% rownames(rna_data)) {
      as.numeric(rna_data[gene_name, ])
    } else {
      rep(NA_real_, ncol(pbmc))
    }
    
    atac_value <- if (gene_name %in% rownames(activity_data)) {
      as.numeric(activity_data[gene_name, ])
    } else {
      rep(NA_real_, ncol(pbmc))
    }
    
    data.frame(
      barcode = colnames(pbmc),
      gene = gene_name,
      rna_value = rna_value,
      atac_value = atac_value,
      seurat_clusters = as.character(pbmc$seurat_clusters),
      cell_type = as.character(pbmc$cell_type),
      stringsAsFactors = FALSE
    )
  }
)

gene_multiome_long_web <- data.table::rbindlist(
  gene_table_list,
  use.names = TRUE,
  fill = TRUE
)

data.table::fwrite(
  gene_multiome_long_web,
  file.path(web_dir, "03_gene_multiome_long_web.csv")
)

# ============================================================
# 14. Unified correlation summary
# ============================================================

correlation_rows <- list()
correlation_index <- 1L

for (gene_name in web_genes_found) {
  gene_df <- gene_multiome_long_web[gene == gene_name]
  
  for (method_name in c("pearson", "spearman")) {
    cor_result <- safe_cor_test(
      gene_df$rna_value,
      gene_df$atac_value,
      method = method_name
    )
    
    correlation_rows[[correlation_index]] <- data.frame(
      relationship = "RNA_vs_ATAC",
      state = NA_character_,
      gene = gene_name,
      cor_result,
      stringsAsFactors = FALSE
    )
    
    correlation_index <- correlation_index + 1L
  }
}

for (state_name in names(state_genes_list)) {
  score_col <- paste0("RNA_Score_", state_name)
  
  genes_to_test <- intersect(
    state_genes_list[[state_name]],
    rownames(rna_data)
  )
  
  for (gene_name in genes_to_test) {
    gene_expression <- as.numeric(
      rna_data[gene_name, ]
    )
    
    state_score <- pbmc@meta.data[[score_col]]
    
    for (method_name in c("pearson", "spearman")) {
      cor_result <- safe_cor_test(
        gene_expression,
        state_score,
        method = method_name
      )
      
      correlation_rows[[correlation_index]] <- data.frame(
        relationship = "State_score_vs_RNA",
        state = state_name,
        gene = gene_name,
        cor_result,
        stringsAsFactors = FALSE
      )
      
      correlation_index <- correlation_index + 1L
    }
  }
}

correlation_summary_web <- dplyr::bind_rows(
  correlation_rows
)

data.table::fwrite(
  correlation_summary_web,
  file.path(web_dir, "04_correlation_summary_web.csv")
)

# ============================================================
# 15. Top peak accessibility table
# ============================================================

DefaultAssay(pbmc) <- "peaks"

peak_counts <- GetAssayData(
  pbmc,
  assay = "peaks",
  layer = "counts"
)

peak_means <- Matrix::rowMeans(peak_counts)

n_top_peaks_actual <- min(
  n_top_peaks_web,
  length(peak_means)
)

top_peaks <- names(
  sort(peak_means, decreasing = TRUE)
)[seq_len(n_top_peaks_actual)]

top_peak_table_list <- lapply(
  top_peaks,
  function(peak_name) {
    peak_value <- as.numeric(
      peak_counts[peak_name, ]
    )
    
    data.frame(
      barcode = colnames(pbmc),
      peak = peak_name,
      peak_count = peak_value,
      log1p_peak_count = log1p(peak_value),
      seurat_clusters = as.character(pbmc$seurat_clusters),
      cell_type = as.character(pbmc$cell_type),
      stringsAsFactors = FALSE
    )
  }
)

top_peak_accessibility_long_web <- data.table::rbindlist(
  top_peak_table_list,
  use.names = TRUE,
  fill = TRUE
)

data.table::fwrite(
  top_peak_accessibility_long_web,
  file.path(web_dir, "05_top_peak_accessibility_long_web.csv")
)

# ============================================================
# 16. Fragment-length histogram table
# ============================================================

message("Building fragment-length histogram data...")

fragment_hist_raw <- scan_fragment_region(
  fragment_path = sample_info$fragment_path[1],
  region_gr = fragment_hist_region
)

cell_group_lookup <- data.table(
  barcode_clean = normalize_barcode(colnames(pbmc)),
  cell_type = as.character(pbmc$cell_type)
)

fragment_hist_data <- merge(
  fragment_hist_raw,
  cell_group_lookup,
  by = "barcode_clean",
  all = FALSE
)

fragment_histogram_web <- fragment_hist_data[
  fragment_length >= 1 &
    fragment_length <= fragment_hist_max_length,
  .(
    fragment_count = sum(count, na.rm = TRUE)
  ),
  by = .(
    cell_type,
    fragment_length
  )
]

if (nrow(fragment_histogram_web) > 0) {
  fragment_histogram_web[
    ,
    fraction_within_cell_type :=
      fragment_count / sum(fragment_count),
    by = cell_type
  ]
  
  data.table::setorder(
    fragment_histogram_web,
    cell_type,
    fragment_length
  )
}

data.table::fwrite(
  fragment_histogram_web,
  file.path(web_dir, "06_fragment_length_histogram_web.csv")
)

# ============================================================
# 17. Differential accessibility
# ============================================================

message("Running differential accessibility analysis...")

DefaultAssay(pbmc) <- "peaks"
Idents(pbmc) <- "seurat_clusters"

available_clusters <- levels(Idents(pbmc))

if (
  da_ident_1_requested %in% available_clusters &&
  da_ident_2_requested %in% available_clusters
) {
  da_ident_1 <- da_ident_1_requested
  da_ident_2 <- da_ident_2_requested
} else if (length(available_clusters) >= 2) {
  da_ident_1 <- available_clusters[1]
  da_ident_2 <- available_clusters[2]
  
  warning(
    "Requested DA clusters were not both present. Using ",
    da_ident_1, " versus ", da_ident_2, "."
  )
} else {
  stop("Differential accessibility requires at least two clusters.")
}

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = da_ident_1,
  ident.2 = da_ident_2,
  test.use = "LR",
  latent.vars = "nCount_peaks",
  min.pct = 0.1
)

da_peaks_output <- da_peaks %>%
  tibble::rownames_to_column("peak")

fc_col <- intersect(
  c("avg_log2FC", "avg_logFC"),
  names(da_peaks_output)
)[1]

if (is.na(fc_col)) {
  stop("Cannot find avg_log2FC or avg_logFC in FindMarkers output.")
}

da_peaks_output <- da_peaks_output %>%
  mutate(
    comparison = paste0(da_ident_1, "_vs_", da_ident_2),
    significance = case_when(
      p_val_adj < 0.05 &
        .data[[fc_col]] > 0 ~
        paste0("Higher_in_", da_ident_1),
      p_val_adj < 0.05 &
        .data[[fc_col]] < 0 ~
        paste0("Higher_in_", da_ident_2),
      TRUE ~ "Not_significant"
    )
  )

data.table::fwrite(
  da_peaks_output,
  file.path(analysis_dir, "07_da_peak_statistics.csv")
)

da_peak_cell_values_web <- data.frame(
  barcode = character(),
  peak = character(),
  peak_count = numeric(),
  log1p_peak_count = numeric(),
  seurat_clusters = character(),
  cell_type = character(),
  atac_umap_1 = numeric(),
  atac_umap_2 = numeric(),
  stringsAsFactors = FALSE
)

if (nrow(da_peaks_output) > 0) {
  da_ranked <- da_peaks_output %>%
    mutate(
      rank_p = ifelse(is.finite(p_val_adj), p_val_adj, 1),
      rank_fc = abs(.data[[fc_col]])
    ) %>%
    arrange(rank_p, desc(rank_fc))
  
  top_da_peaks <- head(
    da_ranked$peak,
    n_da_peaks_web
  )
  
  da_value_list <- lapply(
    top_da_peaks,
    function(peak_name) {
      values <- as.numeric(
        peak_counts[peak_name, ]
      )
      
      data.frame(
        barcode = colnames(pbmc),
        peak = peak_name,
        peak_count = values,
        log1p_peak_count = log1p(values),
        seurat_clusters = as.character(pbmc$seurat_clusters),
        cell_type = as.character(pbmc$cell_type),
        atac_umap_1 = atac_embedding[colnames(pbmc), 1],
        atac_umap_2 = atac_embedding[colnames(pbmc), 2],
        stringsAsFactors = FALSE
      )
    }
  )
  
  da_peak_cell_values_web <- data.table::rbindlist(
    da_value_list,
    use.names = TRUE,
    fill = TRUE
  )
}

data.table::fwrite(
  da_peak_cell_values_web,
  file.path(web_dir, "08_da_peak_cell_values_web.csv")
)

# ============================================================
# 18. Closest genes + GO and KEGG enrichment
# ============================================================

message("Running GO and KEGG enrichment analysis...")

closest_genes_to_da <- data.frame()

functional_enrichment <- data.frame(
  enrichment_type = character(),
  ID = character(),
  Description = character(),
  GeneRatio = character(),
  BgRatio = character(),
  pvalue = numeric(),
  p.adjust = numeric(),
  qvalue = numeric(),
  geneID = character(),
  Count = integer(),
  gene_ratio_numeric = numeric(),
  background_ratio_numeric = numeric(),
  neg_log10_p_adjust = numeric(),
  stringsAsFactors = FALSE
)

if (nrow(da_peaks_output) > 0) {
  significant_open_peaks <- da_peaks_output %>%
    filter(
      p_val_adj < 0.05,
      .data[[fc_col]] > 0
    ) %>%
    pull(peak)
  
  if (length(significant_open_peaks) < 5) {
    warning(
      "Fewer than five significant positive DA peaks were found. ",
      "Using the top positive DA peaks for demonstration."
    )
    
    significant_open_peaks <- da_peaks_output %>%
      filter(.data[[fc_col]] > 0) %>%
      arrange(
        p_val_adj,
        desc(.data[[fc_col]])
      ) %>%
      slice_head(n = 200) %>%
      pull(peak)
  }
  
  if (length(significant_open_peaks) >= 5) {
    closest_genes_to_da <- ClosestFeature(
      pbmc,
      regions = significant_open_peaks
    )
    
    closest_genes_to_da$peak <- rownames(
      closest_genes_to_da
    )
    
    data.table::fwrite(
      closest_genes_to_da,
      file.path(
        analysis_dir,
        "11_closest_genes_to_DA_peaks.csv"
      )
    )
    
    all_tested_closest <- ClosestFeature(
      pbmc,
      regions = da_peaks_output$peak
    )
    
    foreground_ensembl <- unique(
      sub(
        "\\..*$",
        "",
        as.character(closest_genes_to_da$gene_id)
      )
    )
    
    universe_ensembl <- unique(
      sub(
        "\\..*$",
        "",
        as.character(all_tested_closest$gene_id)
      )
    )
    
    foreground_ensembl <- foreground_ensembl[
      !is.na(foreground_ensembl) &
        nzchar(foreground_ensembl)
    ]
    
    universe_ensembl <- universe_ensembl[
      !is.na(universe_ensembl) &
        nzchar(universe_ensembl)
    ]
    
    foreground_map <- suppressMessages(
      clusterProfiler::bitr(
        foreground_ensembl,
        fromType = "ENSEMBL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db
      )
    )
    
    universe_map <- suppressMessages(
      clusterProfiler::bitr(
        universe_ensembl,
        fromType = "ENSEMBL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db
      )
    )
    
    foreground_entrez <- unique(
      foreground_map$ENTREZID
    )
    universe_entrez <- unique(
      universe_map$ENTREZID
    )
    
    foreground_entrez <- foreground_entrez[
      !is.na(foreground_entrez)
    ]
    universe_entrez <- universe_entrez[
      !is.na(universe_entrez)
    ]
    
    go_result <- NULL
    kegg_result <- NULL
    
    if (length(foreground_entrez) >= 5) {
      go_result <- tryCatch(
        enrichGO(
          gene = foreground_entrez,
          universe = universe_entrez,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = go_ontology,
          pAdjustMethod = "BH",
          pvalueCutoff = enrichment_pvalue_cutoff,
          qvalueCutoff = enrichment_qvalue_cutoff,
          readable = TRUE
        ),
        error = function(e) {
          warning(
            "GO enrichment failed: ",
            conditionMessage(e)
          )
          NULL
        }
      )
      
      kegg_result <- tryCatch(
        enrichKEGG(
          gene = foreground_entrez,
          universe = universe_entrez,
          organism = "hsa",
          pAdjustMethod = "BH",
          pvalueCutoff = enrichment_pvalue_cutoff,
          qvalueCutoff = enrichment_qvalue_cutoff
        ),
        error = function(e) {
          warning(
            "KEGG enrichment failed, probably because KEGG was unreachable: ",
            conditionMessage(e)
          )
          NULL
        }
      )
      
      if (
        !is.null(kegg_result) &&
        nrow(as.data.frame(kegg_result)) > 0
      ) {
        kegg_result <- tryCatch(
          setReadable(
            kegg_result,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID"
          ),
          error = function(e) {
            warning(
              "KEGG setReadable failed; Entrez IDs are retained: ",
              conditionMessage(e)
            )
            kegg_result
          }
        )
      }
    }
    
    go_df <- if (
      !is.null(go_result) &&
      nrow(as.data.frame(go_result)) > 0
    ) {
      as.data.frame(go_result) %>%
        mutate(
          enrichment_type = paste0("GO_", go_ontology)
        )
    } else {
      data.frame()
    }
    
    kegg_df <- if (
      !is.null(kegg_result) &&
      nrow(as.data.frame(kegg_result)) > 0
    ) {
      as.data.frame(kegg_result) %>%
        mutate(enrichment_type = "KEGG")
    } else {
      data.frame()
    }
    
    functional_enrichment <- bind_rows(
      go_df,
      kegg_df
    )
    
    if (nrow(functional_enrichment) > 0) {
      functional_enrichment <- functional_enrichment %>%
        mutate(
          gene_ratio_numeric = ratio_to_numeric(GeneRatio),
          background_ratio_numeric = ratio_to_numeric(BgRatio),
          neg_log10_p_adjust = -log10(
            pmax(p.adjust, .Machine$double.xmin)
          )
        )
    }
  }
}

data.table::fwrite(
  functional_enrichment,
  file.path(
    web_dir,
    "09_functional_enrichment_GO_KEGG.csv"
  )
)

# ============================================================
# 19. Top-5 RNA gene genome-track tables
# ============================================================

message("Building top-5 gene web coverage tracks...")

sanitize_filename_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "gene")
}

empty_genome_track_table <- function() {
  data.frame(
    gene_rank = integer(),
    gene = character(),
    selection_method = character(),
    selection_score = numeric(),
    track_type = character(),
    chromosome = character(),
    start = integer(),
    end = integer(),
    position = numeric(),
    group = character(),
    value = numeric(),
    value_name = character(),
    label = character(),
    strand = character(),
    stringsAsFactors = FALSE
  )
}

ann <- Annotation(pbmc)

gene_col_candidates <- c(
  "gene_name",
  "symbol",
  "external_gene_name",
  "gene_id"
)

gene_col <- gene_col_candidates[
  gene_col_candidates %in%
    colnames(mcols(ann))
][1]

if (is.na(gene_col)) {
  stop(
    "No gene-name column was found in the annotation. ",
    "Expected one of: ",
    paste(gene_col_candidates, collapse = ", ")
  )
}

annotation_gene_names <- as.character(
  mcols(ann)[[gene_col]]
)

annotation_gene_names <- annotation_gene_names[
  !is.na(annotation_gene_names) &
    nzchar(annotation_gene_names)
]

annotation_gene_names_upper <- unique(
  toupper(annotation_gene_names)
)

rna_gene_names <- rownames(rna_data)
rna_gene_names_upper <- toupper(rna_gene_names)

rna_variable_genes <- VariableFeatures(pbmc[["RNA"]])

eligible_variable_genes <- rna_variable_genes[
  toupper(rna_variable_genes) %in%
    annotation_gene_names_upper
]

manual_rna_match <- match(
  toupper(coverage_genes_manual),
  rna_gene_names_upper
)

manual_genes_found <- unique(
  rna_gene_names[
    manual_rna_match[
      !is.na(manual_rna_match)
    ]
  ]
)

manual_genes_found <- manual_genes_found[
  toupper(manual_genes_found) %in%
    annotation_gene_names_upper
]

if (length(coverage_genes_manual) > 0) {
  missing_manual_genes <- setdiff(
    toupper(coverage_genes_manual),
    toupper(manual_genes_found)
  )
  
  if (length(missing_manual_genes) > 0) {
    warning(
      "The following manually requested coverage genes were not found ",
      "in both RNA and annotation: ",
      paste(missing_manual_genes, collapse = ", ")
    )
  }
}

if (length(manual_genes_found) > 0) {
  selected_gene_candidates <- unique(
    c(
      manual_genes_found,
      eligible_variable_genes
    )
  )
  selection_method <- "manual_then_Seurat_RNA_variable_feature_rank"
} else {
  selected_gene_candidates <- unique(eligible_variable_genes)
  selection_method <- "Seurat_RNA_variable_feature_rank"
}

if (length(selected_gene_candidates) < coverage_top_n_genes) {
  mean_rna_all <- Matrix::rowMeans(rna_data)
  
  fallback_order <- names(
    sort(mean_rna_all, decreasing = TRUE)
  )
  
  fallback_genes <- fallback_order[
    toupper(fallback_order) %in%
      annotation_gene_names_upper
  ]
  
  selected_gene_candidates <- unique(
    c(selected_gene_candidates, fallback_genes)
  )
  
  if (selection_method != "manual") {
    selection_method <-
      "Seurat_RNA_variable_feature_rank_then_mean_expression"
  }
}

top_coverage_genes <- head(
  selected_gene_candidates,
  coverage_top_n_genes
)

if (length(top_coverage_genes) == 0) {
  stop(
    "No RNA genes could be matched to the local annotation ",
    "for genome-track construction."
  )
}

top_gene_mean_expression <- Matrix::rowMeans(
  rna_data[top_coverage_genes, , drop = FALSE]
)

top_gene_selection <- data.frame(
  gene_rank = seq_along(top_coverage_genes),
  gene = top_coverage_genes,
  selection_method = selection_method,
  mean_rna_expression = as.numeric(
    top_gene_mean_expression[top_coverage_genes]
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  top_gene_selection,
  file.path(
    web_dir,
    "10_top5_gene_selection.csv"
  )
)

build_one_gene_genome_track <- function(
    gene_name,
    gene_rank,
    selection_method,
    selection_score
) {
  gene_ann <- ann[
    toupper(
      as.character(
        mcols(ann)[[gene_col]]
      )
    ) == toupper(gene_name)
  ]
  
  if (length(gene_ann) == 0) {
    warning("Annotation was not found for coverage gene: ", gene_name)
    return(empty_genome_track_table())
  }
  
  gene_widths <- width(gene_ann)
  representative_index <- which.max(gene_widths)
  gene_chromosome <- as.character(
    seqnames(gene_ann)
  )[representative_index]
  
  gene_ann <- gene_ann[
    as.character(seqnames(gene_ann)) == gene_chromosome
  ]
  
  gene_start <- min(start(gene_ann))
  gene_end <- max(end(gene_ann))
  
  region_start <- max(
    1,
    gene_start - coverage_extend_bp
  )
  region_end <- gene_end + coverage_extend_bp
  
  gene_region <- GRanges(
    seqnames = gene_chromosome,
    ranges = IRanges(
      start = region_start,
      end = region_end
    )
  )
  
  gene_fragments <- scan_fragment_region(
    sample_info$fragment_path[1],
    gene_region
  )
  
  gene_fragments <- merge(
    gene_fragments,
    cell_group_lookup,
    by = "barcode_clean",
    all = FALSE
  )
  
  n_bins <- ceiling(
    (region_end - region_start + 1) /
      coverage_bin_width
  )
  
  cell_type_sizes <- data.table(
    cell_type = as.character(pbmc$cell_type)
  )[
    ,
    .(
      n_cells = .N
    ),
    by = cell_type
  ]
  
  coverage_grid <- data.table::CJ(
    cell_type = sort(
      unique(as.character(pbmc$cell_type))
    ),
    bin_id = seq_len(n_bins),
    unique = TRUE
  )
  
  if (nrow(gene_fragments) > 0) {
    insertion_start <- gene_fragments[
      ,
      .(
        position = start,
        weight = count,
        cell_type
      )
    ]
    
    insertion_end <- gene_fragments[
      ,
      .(
        position = end,
        weight = count,
        cell_type
      )
    ]
    
    insertion_events <- rbind(
      insertion_start,
      insertion_end,
      use.names = TRUE
    )
    
    insertion_events <- insertion_events[
      position >= region_start &
        position <= region_end
    ]
    
    insertion_events[
      ,
      bin_id := pmin(
        n_bins,
        floor(
          (position - region_start) /
            coverage_bin_width
        ) +
          1L
      )
    ]
    
    coverage_counts <- insertion_events[
      ,
      .(
        raw_insertions = sum(
          weight,
          na.rm = TRUE
        )
      ),
      by = .(
        cell_type,
        bin_id
      )
    ]
  } else {
    coverage_counts <- data.table(
      cell_type = character(),
      bin_id = integer(),
      raw_insertions = numeric()
    )
  }
  
  coverage_grid <- merge(
    coverage_grid,
    coverage_counts,
    by = c("cell_type", "bin_id"),
    all.x = TRUE
  )
  
  coverage_grid <- merge(
    coverage_grid,
    cell_type_sizes,
    by = "cell_type",
    all.x = TRUE
  )
  
  coverage_grid[
    is.na(raw_insertions),
    raw_insertions := 0
  ]
  
  coverage_grid[
    ,
    bin_start :=
      region_start +
      (bin_id - 1L) * coverage_bin_width
  ]
  
  coverage_grid[
    ,
    bin_end := pmin(
      region_end,
      bin_start + coverage_bin_width - 1L
    )
  ]
  
  coverage_grid[
    ,
    position := (bin_start + bin_end) / 2
  ]
  
  coverage_grid[
    ,
    insertions_per_cell :=
      raw_insertions / pmax(n_cells, 1)
  ]
  
  coverage_track <- coverage_grid %>%
    transmute(
      gene_rank = as.integer(gene_rank),
      gene = gene_name,
      selection_method = selection_method,
      selection_score = selection_score,
      track_type = "coverage",
      chromosome = gene_chromosome,
      start = bin_start,
      end = bin_end,
      position = position,
      group = cell_type,
      value = insertions_per_cell,
      value_name = "insertions_per_cell",
      label = NA_character_,
      strand = NA_character_
    )
  
  peak_granges <- granges(pbmc[["peaks"]])
  peak_overlap <- subsetByOverlaps(
    peak_granges,
    gene_region
  )
  
  peak_track <- if (length(peak_overlap) > 0) {
    data.frame(
      gene_rank = as.integer(gene_rank),
      gene = gene_name,
      selection_method = selection_method,
      selection_score = selection_score,
      track_type = "peak",
      chromosome = as.character(seqnames(peak_overlap)),
      start = start(peak_overlap),
      end = end(peak_overlap),
      position = (
        start(peak_overlap) + end(peak_overlap)
      ) / 2,
      group = NA_character_,
      value = NA_real_,
      value_name = NA_character_,
      label = paste0(
        as.character(seqnames(peak_overlap)),
        ":",
        start(peak_overlap),
        "-",
        end(peak_overlap)
      ),
      strand = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    empty_genome_track_table()
  }
  
  gene_strand <- as.character(
    unique(strand(gene_ann))
  )[1]
  
  gene_track <- data.frame(
    gene_rank = as.integer(gene_rank),
    gene = gene_name,
    selection_method = selection_method,
    selection_score = selection_score,
    track_type = "gene",
    chromosome = gene_chromosome,
    start = gene_start,
    end = gene_end,
    position = (
      gene_start + gene_end
    ) / 2,
    group = NA_character_,
    value = NA_real_,
    value_name = NA_character_,
    label = gene_name,
    strand = gene_strand,
    stringsAsFactors = FALSE
  )
  
  bind_rows(
    coverage_track,
    peak_track,
    gene_track
  )
}

genome_track_list <- lapply(
  seq_len(nrow(top_gene_selection)),
  function(i) {
    build_one_gene_genome_track(
      gene_name = top_gene_selection$gene[i],
      gene_rank = top_gene_selection$gene_rank[i],
      selection_method = top_gene_selection$selection_method[i],
      selection_score = top_gene_selection$mean_rna_expression[i]
    )
  }
)

genome_tracks_web <- bind_rows(genome_track_list)

data.table::fwrite(
  genome_tracks_web,
  file.path(
    web_dir,
    "10_top5_genome_tracks_web.csv"
  )
)

for (i in seq_len(nrow(top_gene_selection))) {
  gene_name <- top_gene_selection$gene[i]
  gene_rank <- top_gene_selection$gene_rank[i]
  gene_file_label <- sanitize_filename_component(gene_name)
  
  one_gene_table <- genome_tracks_web %>%
    filter(
      .data$gene_rank == .env$gene_rank,
      .data$gene == .env$gene_name
    )
  
  data.table::fwrite(
    one_gene_table,
    file.path(
      web_dir,
      sprintf(
        "10_top5_genome_tracks_rank%02d_%s.csv",
        gene_rank,
        gene_file_label
      )
    )
  )
}

# ============================================================
# 20. QC and parameter summary
# ============================================================

qc_summary <- data.frame(
  metric = c(
    "cells_before_qc",
    "cells_after_qc",
    "retention_fraction",
    "qc_min_nCount_peaks",
    "qc_max_nCount_peaks",
    "qc_max_blacklist_ratio",
    "qc_max_nucleosome_signal",
    "qc_min_TSS_enrichment",
    "cluster_resolution",
    "DA_ident_1",
    "DA_ident_2"
  ),
  value = c(
    n_cells_before_qc,
    n_cells_after_qc,
    n_cells_after_qc / n_cells_before_qc,
    qc_min_nCount_peaks,
    qc_max_nCount_peaks,
    qc_max_blacklist_ratio,
    qc_max_nucleosome_signal,
    qc_min_TSS_enrichment,
    cluster_resolution,
    da_ident_1,
    da_ident_2
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  qc_summary,
  file.path(
    analysis_dir,
    "12_pipeline_parameters_and_qc_summary.csv"
  )
)

# ============================================================
# 21. Draw example figures
# ============================================================

p_qc_density <- ggplot(
  cell_metadata_web,
  aes(
    x = nCount_peaks,
    y = TSS.enrichment,
    color = pass_qc
  )
) +
  geom_point(alpha = 0.35, size = 0.45) +
  scale_x_log10() +
  labs(
    title = "ATAC QC: Peak Fragments versus TSS Enrichment",
    x = "Peak-region fragments (log10 scale)",
    y = "TSS enrichment",
    color = "Pass QC"
  ) +
  theme_web()

save_plot(
  p_qc_density,
  "01_qc_density_nCount_TSS.png",
  7,
  5.5
)

qc_metrics <- c(
  "nCount_peaks",
  "TSS.enrichment",
  "blacklist_ratio",
  "nucleosome_signal",
  "pct_reads_in_peaks"
)

qc_violin_long <- cell_metadata_web %>%
  select(
    timepoint,
    pass_qc,
    all_of(qc_metrics)
  ) %>%
  pivot_longer(
    cols = all_of(qc_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(is.finite(value))

p_qc_violin <- ggplot(
  qc_violin_long,
  aes(
    x = timepoint,
    y = value,
    fill = pass_qc
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.75
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    title = "ATAC QC Metrics",
    x = "Timepoint",
    y = "Value",
    fill = "Pass QC"
  ) +
  theme_web()

save_plot(
  p_qc_violin,
  "02_qc_violin_all_metrics.png",
  12,
  7
)

if (nrow(fragment_histogram_web) > 0) {
  p_fragment <- ggplot(
    fragment_histogram_web,
    aes(
      x = fragment_length,
      y = fraction_within_cell_type,
      color = cell_type
    )
  ) +
    geom_line(linewidth = 0.45, alpha = 0.85) +
    labs(
      title = "ATAC Fragment-length Distribution",
      x = "Fragment length (bp)",
      y = "Fraction within cell type",
      color = "Cell type"
    ) +
    theme_web()
  
  save_plot(
    p_fragment,
    "03_fragment_length_histogram.png",
    9,
    5.5
  )
}

plot_cells <- cell_metadata_web %>%
  filter(pass_qc %in% TRUE)

p_rna_umap <- make_umap_plot(
  data = plot_cells,
  x_col = "rna_umap_1",
  y_col = "rna_umap_2",
  color_col = "seurat_clusters",
  label_col = "seurat_clusters",
  title = "RNA UMAP Clusters",
  x_label = "RNA UMAP 1",
  y_label = "RNA UMAP 2"
)

save_plot(
  p_rna_umap,
  "04_rna_umap_clusters.png",
  7.5,
  6
)

p_atac_umap <- make_umap_plot(
  data = plot_cells,
  x_col = "atac_umap_1",
  y_col = "atac_umap_2",
  color_col = "cell_type",
  label_col = "cell_type",
  title = "ATAC UMAP: Auto-annotated Cell Types",
  x_label = "ATAC UMAP 1",
  y_label = "ATAC UMAP 2"
)

save_plot(
  p_atac_umap,
  "05_atac_umap_celltypes.png",
  8.5,
  6
)

p_wnn_umap <- make_umap_plot(
  data = plot_cells,
  x_col = "wnn_umap_1",
  y_col = "wnn_umap_2",
  color_col = "cell_type",
  label_col = "cell_type",
  title = "WNN UMAP: Integrated RNA and ATAC",
  x_label = "WNN UMAP 1",
  y_label = "WNN UMAP 2"
)

save_plot(
  p_wnn_umap,
  "06_wnn_umap_celltypes.png",
  8.5,
  6
)

marker_mean_cols <- grep(
  "^mean_marker_score_",
  names(cluster_annotation_scores),
  value = TRUE
)

cluster_marker_long <- cluster_annotation_scores %>%
  select(
    seurat_clusters,
    all_of(marker_mean_cols)
  ) %>%
  pivot_longer(
    cols = all_of(marker_mean_cols),
    names_to = "candidate_cell_type",
    values_to = "mean_marker_score"
  ) %>%
  mutate(
    candidate_cell_type = sub(
      "^mean_marker_score_",
      "",
      candidate_cell_type
    )
  )

p_marker_heatmap <- ggplot(
  cluster_marker_long,
  aes(
    x = candidate_cell_type,
    y = factor(
      seurat_clusters,
      levels = sort(unique(seurat_clusters))
    ),
    fill = mean_marker_score
  )
) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "Cluster-level Candidate Cell-type Marker Scores",
    x = "Candidate cell type",
    y = "Seurat cluster",
    fill = "Mean marker score"
  ) +
  theme_web() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  p_marker_heatmap,
  "07_cluster_marker_score_heatmap.png",
  11,
  6
)

violin_genes <- intersect(
  c("MAP2", "CD3D", "LYZ"),
  unique(gene_multiome_long_web$gene)
)

if (length(violin_genes) > 0) {
  gene_violin_data <- gene_multiome_long_web %>%
    filter(gene %in% violin_genes) %>%
    pivot_longer(
      cols = c(rna_value, atac_value),
      names_to = "modality",
      values_to = "value"
    ) %>%
    mutate(
      modality = recode(
        modality,
        rna_value = "RNA expression",
        atac_value = "ATAC gene activity"
      )
    ) %>%
    filter(is.finite(value))
  
  p_gene_violin <- ggplot(
    gene_violin_data,
    aes(
      x = cell_type,
      y = value,
      fill = cell_type
    )
  ) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      alpha = 0.8
    ) +
    geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.55
    ) +
    facet_grid(
      modality ~ gene,
      scales = "free_y"
    ) +
    labs(
      title = "RNA Expression and ATAC Gene Activity",
      x = "Cell type",
      y = "Normalized value"
    ) +
    theme_web() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  save_plot(
    p_gene_violin,
    "08_gene_RNA_ATAC_violin.png",
    14,
    8
  )
}

feature_gene_candidates <- c(
  "MAP2",
  "CD3D",
  "LYZ",
  "MKI67"
)

feature_gene <- feature_gene_candidates[
  feature_gene_candidates %in%
    unique(gene_multiome_long_web$gene)
][1]

if (!is.na(feature_gene)) {
  feature_data <- gene_multiome_long_web %>%
    filter(gene == feature_gene) %>%
    left_join(
      plot_cells %>%
        select(
          barcode,
          rna_umap_1,
          rna_umap_2,
          atac_umap_1,
          atac_umap_2
        ),
      by = "barcode"
    )
  
  p_feature_rna <- ggplot(
    feature_data,
    aes(
      x = rna_umap_1,
      y = rna_umap_2,
      color = rna_value
    )
  ) +
    geom_point(size = 0.35, alpha = 0.8) +
    scale_color_viridis_c(na.value = "grey90") +
    labs(
      title = paste0(feature_gene, " RNA expression"),
      x = "RNA UMAP 1",
      y = "RNA UMAP 2",
      color = "RNA"
    ) +
    theme_web()
  
  p_feature_atac <- ggplot(
    feature_data,
    aes(
      x = atac_umap_1,
      y = atac_umap_2,
      color = atac_value
    )
  ) +
    geom_point(size = 0.35, alpha = 0.8) +
    scale_color_viridis_c(na.value = "grey90") +
    labs(
      title = paste0(feature_gene, " ATAC gene activity"),
      x = "ATAC UMAP 1",
      y = "ATAC UMAP 2",
      color = "ATAC"
    ) +
    theme_web()
  
  save_plot(
    p_feature_rna + p_feature_atac,
    "09_gene_RNA_ATAC_featureplots.png",
    12,
    5.5
  )
}

state_score_cols <- c(
  paste0("RNA_Score_", names(state_genes_list)),
  paste0("ATAC_Score_", names(state_genes_list))
)

state_score_cols <- intersect(
  state_score_cols,
  names(plot_cells)
)

state_score_long <- plot_cells %>%
  select(
    barcode,
    cell_type,
    all_of(state_score_cols)
  ) %>%
  pivot_longer(
    cols = all_of(state_score_cols),
    names_to = "score_name",
    values_to = "score"
  ) %>%
  mutate(
    modality = ifelse(
      grepl("^RNA_", score_name),
      "RNA",
      "ATAC"
    ),
    state = sub(
      "^(RNA|ATAC)_Score_",
      "",
      score_name
    )
  )

p_state_box <- ggplot(
  state_score_long,
  aes(
    x = cell_type,
    y = score,
    fill = cell_type
  )
) +
  geom_boxplot(
    outlier.size = 0.25,
    linewidth = 0.35
  ) +
  facet_grid(
    modality ~ state,
    scales = "free_y"
  ) +
  labs(
    title = "RNA and ATAC Cell-state Scores",
    x = "Cell type",
    y = "State score"
  ) +
  theme_web() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

save_plot(
  p_state_box,
  "10_state_score_boxplot_all.png",
  18,
  8
)

state_score_heatmap_data <- state_score_long %>%
  group_by(
    cell_type,
    modality,
    state
  ) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    state_modality = paste(
      modality,
      state,
      sep = ": "
    )
  )

p_state_heatmap <- ggplot(
  state_score_heatmap_data,
  aes(
    x = state_modality,
    y = cell_type,
    fill = mean_score
  )
) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "Mean Cell-state Scores by Cell Type",
    x = "Modality and state",
    y = "Cell type",
    fill = "Mean score"
  ) +
  theme_web() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  p_state_heatmap,
  "11_state_score_heatmap.png",
  12,
  7
)

state_gene_examples <- data.frame(
  state = c(
    "Proliferation",
    "Hypoxia",
    "Immune_Activity"
  ),
  gene = c(
    "MKI67",
    "HIF1A",
    "GZMB"
  ),
  stringsAsFactors = FALSE
)

state_gene_plot_list <- list()

for (i in seq_len(nrow(state_gene_examples))) {
  state_name <- state_gene_examples$state[i]
  gene_name <- state_gene_examples$gene[i]
  score_col <- paste0("RNA_Score_", state_name)
  
  if (
    gene_name %in% unique(gene_multiome_long_web$gene) &&
    score_col %in% names(plot_cells)
  ) {
    scatter_data <- gene_multiome_long_web %>%
      filter(gene == gene_name) %>%
      select(
        barcode,
        gene,
        rna_value,
        cell_type
      ) %>%
      left_join(
        plot_cells %>%
          select(
            barcode,
            all_of(score_col)
          ),
        by = "barcode"
      )
    
    names(scatter_data)[
      names(scatter_data) == score_col
    ] <- "state_score"
    
    p_one <- ggplot(
      scatter_data,
      aes(
        x = rna_value,
        y = state_score,
        color = cell_type
      )
    ) +
      geom_point(size = 0.45, alpha = 0.45) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        color = "black",
        linewidth = 0.6
      ) +
      labs(
        title = paste0(state_name, ": ", gene_name),
        x = paste0(gene_name, " RNA expression"),
        y = paste0(state_name, " RNA score"),
        color = "Cell type"
      ) +
      theme_web() +
      theme(legend.position = "none")
    
    state_gene_plot_list[[
      length(state_gene_plot_list) + 1
    ]] <- p_one
  }
}

if (length(state_gene_plot_list) > 0) {
  save_plot(
    patchwork::wrap_plots(
      state_gene_plot_list,
      ncol = length(state_gene_plot_list)
    ),
    "12_state_gene_scatter_examples.png",
    5.5 * length(state_gene_plot_list),
    5
  )
}

rna_atac_example_genes <- intersect(
  c("MAP2", "CD3D", "LYZ", "MKI67"),
  unique(gene_multiome_long_web$gene)
)

rna_atac_plot_list <- lapply(
  rna_atac_example_genes,
  function(gene_name) {
    scatter_data <- gene_multiome_long_web %>%
      filter(
        gene == gene_name,
        is.finite(rna_value),
        is.finite(atac_value)
      )
    
    ggplot(
      scatter_data,
      aes(
        x = atac_value,
        y = rna_value,
        color = cell_type
      )
    ) +
      geom_point(size = 0.45, alpha = 0.45) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        color = "black",
        linewidth = 0.6
      ) +
      labs(
        title = gene_name,
        x = "ATAC gene activity",
        y = "RNA expression",
        color = "Cell type"
      ) +
      theme_web() +
      theme(legend.position = "none")
  }
)

if (length(rna_atac_plot_list) > 0) {
  save_plot(
    patchwork::wrap_plots(
      rna_atac_plot_list,
      ncol = min(2, length(rna_atac_plot_list))
    ),
    "13_RNA_ATAC_gene_scatter_examples.png",
    11,
    5 * ceiling(length(rna_atac_plot_list) / 2)
  )
}

correlation_plot_data <- correlation_summary_web %>%
  filter(
    method == "Pearson",
    is.finite(correlation)
  ) %>%
  mutate(
    abs_correlation = abs(correlation),
    feature_label = ifelse(
      is.na(state),
      gene,
      paste0(state, ": ", gene)
    )
  ) %>%
  group_by(relationship) %>%
  slice_max(
    order_by = abs_correlation,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup()

if (nrow(correlation_plot_data) > 0) {
  p_correlation <- ggplot(
    correlation_plot_data,
    aes(
      x = reorder(feature_label, correlation),
      y = correlation,
      fill = relationship
    )
  ) +
    geom_col() +
    coord_flip() +
    facet_wrap(
      ~ relationship,
      scales = "free_y"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    labs(
      title = "Top Pearson Correlations",
      x = NULL,
      y = "Pearson correlation"
    ) +
    theme_web() +
    theme(legend.position = "none")
  
  save_plot(
    p_correlation,
    "14_correlation_summary.png",
    11,
    8
  )
}

top_peak_heatmap_data <-
  top_peak_accessibility_long_web %>%
  group_by(
    peak,
    cell_type
  ) %>%
  summarise(
    mean_log1p_peak_count = mean(
      log1p_peak_count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p_top_peak_heatmap <- ggplot(
  top_peak_heatmap_data,
  aes(
    x = cell_type,
    y = peak,
    fill = mean_log1p_peak_count
  )
) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "Top Peak Accessibility by Cell Type",
    x = "Cell type",
    y = "Peak",
    fill = "Mean log1p count"
  ) +
  theme_web() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6)
  )

save_plot(
  p_top_peak_heatmap,
  "15_top_peak_celltype_heatmap.png",
  12,
  9
)

if (nrow(da_peaks_output) > 0) {
  da_volcano_data <- da_peaks_output %>%
    mutate(
      neg_log10_padj = -log10(
        pmax(p_val_adj, .Machine$double.xmin)
      )
    )
  
  p_da_volcano <- ggplot(
    da_volcano_data,
    aes(
      x = .data[[fc_col]],
      y = neg_log10_padj,
      color = significance
    )
  ) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(
      title = paste0(
        "Differential Accessibility: Cluster ",
        da_ident_1,
        " versus ",
        da_ident_2
      ),
      x = fc_col,
      y = "-log10(adjusted P value)",
      color = "Significance"
    ) +
    theme_web()
  
  save_plot(
    p_da_volcano,
    "16_DA_peaks_volcano.png",
    7,
    5.5
  )
  
  if (nrow(da_peak_cell_values_web) > 0) {
    best_peak <- unique(
      da_peak_cell_values_web$peak
    )[1]
    
    best_peak_data <- da_peak_cell_values_web %>%
      filter(peak == best_peak)
    
    p_da_violin <- ggplot(
      best_peak_data,
      aes(
        x = seurat_clusters,
        y = log1p_peak_count,
        fill = seurat_clusters
      )
    ) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        alpha = 0.55
      ) +
      labs(
        title = paste0("Best DA Peak: ", best_peak),
        x = "Seurat cluster",
        y = "log1p peak count"
      ) +
      theme_web() +
      theme(legend.position = "none")
    
    save_plot(
      p_da_violin,
      "17_DA_best_peak_violin.png",
      8,
      5.5
    )
    
    p_da_feature <- ggplot(
      best_peak_data,
      aes(
        x = atac_umap_1,
        y = atac_umap_2,
        color = log1p_peak_count
      )
    ) +
      geom_point(size = 0.4, alpha = 0.8) +
      scale_color_viridis_c() +
      labs(
        title = paste0("ATAC UMAP: ", best_peak),
        x = "ATAC UMAP 1",
        y = "ATAC UMAP 2",
        color = "log1p count"
      ) +
      theme_web()
    
    save_plot(
      p_da_feature,
      "18_DA_best_peak_featureplot.png",
      7,
      5.5
    )
  }
}

if (nrow(functional_enrichment) > 0) {
  go_type <- paste0("GO_", go_ontology)
  
  p_go <- make_enrichment_barplot(
    functional_enrichment,
    type_value = go_type,
    title = paste0("GO ", go_ontology, " Enrichment")
  )
  
  if (!is.null(p_go)) {
    save_plot(
      p_go,
      "19_GO_enrichment_barplot.png",
      9,
      6
    )
  }
  
  p_kegg <- make_enrichment_barplot(
    functional_enrichment,
    type_value = "KEGG",
    title = "KEGG Pathway Enrichment"
  )
  
  if (!is.null(p_kegg)) {
    save_plot(
      p_kegg,
      "20_KEGG_enrichment_barplot.png",
      9,
      6
    )
  }
  
  combined_enrichment_plot_data <-
    functional_enrichment %>%
    filter(is.finite(p.adjust)) %>%
    group_by(enrichment_type) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(unique(Description))
      )
    )
  
  if (nrow(combined_enrichment_plot_data) > 0) {
    p_enrichment_combined <- ggplot(
      combined_enrichment_plot_data,
      aes(
        x = gene_ratio_numeric,
        y = Description,
        size = Count,
        color = neg_log10_p_adjust
      )
    ) +
      geom_point(alpha = 0.9) +
      facet_wrap(
        ~ enrichment_type,
        scales = "free_y"
      ) +
      scale_color_viridis_c() +
      labs(
        title = "GO and KEGG Enrichment",
        x = "Gene ratio",
        y = NULL,
        size = "Gene count",
        color = "-log10 adjusted P"
      ) +
      theme_web()
    
    save_plot(
      p_enrichment_combined,
      "21_GO_KEGG_enrichment_dotplot.png",
      13,
      8
    )
  }
}

coverage_figure_mapping <- tibble::tibble(
  figure_file = character(),
  table_file = character(),
  key_fields = character(),
  description = character()
)

if (nrow(genome_tracks_web) > 0) {
  for (i in seq_len(nrow(top_gene_selection))) {
    gene_name <- top_gene_selection$gene[i]
    gene_rank <- top_gene_selection$gene_rank[i]
    gene_file_label <- sanitize_filename_component(gene_name)
    
    one_gene_tracks <- genome_tracks_web %>%
      filter(
        .data$gene_rank == .env$gene_rank,
        .data$gene == .env$gene_name
      )
    
    coverage_plot_data <- one_gene_tracks %>%
      filter(track_type == "coverage")
    
    peak_plot_data <- one_gene_tracks %>%
      filter(track_type == "peak")
    
    gene_plot_data <- one_gene_tracks %>%
      filter(track_type == "gene")
    
    p_coverage <- ggplot(
      coverage_plot_data,
      aes(
        x = position,
        y = value,
        color = group
      )
    ) +
      geom_line(linewidth = 0.55) +
      labs(
        title = paste0(
          "Top ", gene_rank, " gene ",
          gene_name,
          ": Aggregated ATAC Insertions"
        ),
        x = NULL,
        y = "Insertions per cell",
        color = "Cell type"
      ) +
      theme_web()
    
    if (nrow(peak_plot_data) > 0) {
      p_peak_track <- ggplot(
        peak_plot_data
      ) +
        geom_segment(
          aes(
            x = start,
            xend = end,
            y = 1,
            yend = 1
          ),
          linewidth = 2
        ) +
        labs(
          x = NULL,
          y = "Peaks"
        ) +
        theme_web() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
    } else {
      p_peak_track <- ggplot() +
        labs(
          x = NULL,
          y = "Peaks"
        ) +
        theme_web()
    }
    
    p_gene_track <- ggplot(
      gene_plot_data
    ) +
      geom_segment(
        aes(
          x = start,
          xend = end,
          y = 1,
          yend = 1
        ),
        linewidth = 2
      ) +
      geom_text(
        aes(
          x = position,
          y = 1.15,
          label = label
        ),
        fontface = "bold"
      ) +
      labs(
        x = "Genomic position",
        y = "Gene"
      ) +
      theme_web() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    
    p_gene_tracks <-
      p_coverage /
      p_peak_track /
      p_gene_track +
      plot_layout(
        heights = c(4, 1, 1)
      )
    
    figure_filename <- sprintf(
      "22_top5_gene_rank%02d_%s_coverage_tracks.png",
      gene_rank,
      gene_file_label
    )
    
    table_filename <- sprintf(
      "10_top5_genome_tracks_rank%02d_%s.csv",
      gene_rank,
      gene_file_label
    )
    
    save_plot(
      p_gene_tracks,
      figure_filename,
      13,
      8
    )
    
    coverage_figure_mapping <- bind_rows(
      coverage_figure_mapping,
      tibble::tibble(
        figure_file = figure_filename,
        table_file = table_filename,
        key_fields = paste(
          "gene_rank; gene; track_type; start; end;",
          "position; group; value"
        ),
        description = paste0(
          "Top-", gene_rank, " RNA gene ",
          gene_name,
          " coverage, peak, and gene tracks"
        )
      )
    )
  }
}

# ============================================================
# 22. Figure-to-table mapping
# ============================================================

figure_table_mapping <- tibble::tribble(
  ~figure_file, ~table_file, ~key_fields, ~description,
  "01_qc_density_nCount_TSS.png",
  "01_cell_metadata_web.csv",
  "nCount_peaks; TSS.enrichment; pass_qc",
  "QC density scatter",
  "02_qc_violin_all_metrics.png",
  "01_cell_metadata_web.csv",
  "timepoint; pass_qc; QC metric columns",
  "QC violin plots",
  "03_fragment_length_histogram.png",
  "06_fragment_length_histogram_web.csv",
  "cell_type; fragment_length; fraction_within_cell_type",
  "Fragment-length distribution",
  "04_rna_umap_clusters.png",
  "01_cell_metadata_web.csv",
  "rna_umap_1; rna_umap_2; seurat_clusters",
  "RNA UMAP",
  "05_atac_umap_celltypes.png",
  "01_cell_metadata_web.csv",
  "atac_umap_1; atac_umap_2; cell_type",
  "ATAC UMAP",
  "06_wnn_umap_celltypes.png",
  "01_cell_metadata_web.csv",
  "wnn_umap_1; wnn_umap_2; cell_type",
  "WNN UMAP",
  "07_cluster_marker_score_heatmap.png",
  "02_cluster_annotation_scores_web.csv",
  "seurat_clusters; mean_marker_score_*",
  "Automatic annotation marker-score heatmap",
  "08_gene_RNA_ATAC_violin.png",
  "03_gene_multiome_long_web.csv",
  "gene; cell_type; rna_value; atac_value",
  "RNA and ATAC violin plots",
  "09_featureplot_RNA_ATAC.png",
  "01_cell_metadata_web.csv + 03_gene_multiome_long_web.csv",
  "barcode; RNA/ATAC UMAP coordinates; rna_value; atac_value",
  "RNA and ATAC feature plots",
  "10_marker_score_distribution.png",
  "01_cell_metadata_web.csv",
  "MarkerScore_*; seurat_clusters; cell_type",
  "Cell-level marker-score distributions",
  "11_state_score_heatmap.png",
  "01_cell_metadata_web.csv",
  "cell_type; RNA_Score_*; ATAC_Score_*",
  "State-score heatmap",
  "12_state_gene_scatter_examples.png",
  "01_cell_metadata_web.csv + 03_gene_multiome_long_web.csv",
  "barcode; RNA_Score_*; gene; rna_value",
  "State score versus RNA expression",
  "13_RNA_ATAC_gene_scatter_examples.png",
  "03_gene_multiome_long_web.csv",
  "gene; rna_value; atac_value; cell_type",
  "RNA versus ATAC scatter",
  "14_correlation_summary.png",
  "04_correlation_summary_web.csv",
  "relationship; state; gene; method; correlation",
  "Correlation summary",
  "15_top_peak_celltype_heatmap.png",
  "05_top_peak_accessibility_long_web.csv",
  "peak; cell_type; log1p_peak_count",
  "Top peak heatmap",
  "16_DA_peaks_volcano.png",
  "07_da_peak_statistics.csv",
  paste0("peak; ", fc_col, "; p_val_adj; significance"),
  "DA volcano plot",
  "17_DA_best_peak_violin.png",
  "08_da_peak_cell_values_web.csv",
  "peak; seurat_clusters; log1p_peak_count",
  "Best DA peak violin",
  "18_DA_best_peak_featureplot.png",
  "08_da_peak_cell_values_web.csv",
  "atac_umap_1; atac_umap_2; log1p_peak_count",
  "Best DA peak feature plot",
  "19_GO_enrichment_barplot.png",
  "09_functional_enrichment_GO_KEGG.csv",
  "enrichment_type; Description; p.adjust; Count",
  "GO enrichment",
  "20_KEGG_enrichment_barplot.png",
  "09_functional_enrichment_GO_KEGG.csv",
  "enrichment_type; Description; p.adjust; Count",
  "KEGG enrichment",
  "21_GO_KEGG_enrichment_dotplot.png",
  "09_functional_enrichment_GO_KEGG.csv",
  "enrichment_type; Description; gene_ratio_numeric; Count; neg_log10_p_adjust",
  "Combined GO and KEGG dot plot"
)

figure_table_mapping <- bind_rows(
  figure_table_mapping,
  coverage_figure_mapping
)

data.table::fwrite(
  figure_table_mapping,
  file.path(
    web_dir,
    "00_figure_table_mapping.csv"
  )
)

print(figure_table_mapping, n = Inf)

# ============================================================
# 23. Save object and session information
# ============================================================

saveRDS(
  pbmc,
  file = file.path(
    object_dir,
    "GSM8289638_GSM8289639_multiome_processed_final.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo.txt"
  )
)

message("============================================================")
message("Finished.")
message("Example figures: ", plot_dir)
message("Web plotting data: ", web_dir)
message("Analysis tables: ", analysis_dir)
message("Seurat object: ", object_dir)
message(
  "Figure-to-table mapping: ",
  file.path(web_dir, "00_figure_table_mapping.csv")
)
message("============================================================")
