`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

row_dataset_key <- function(row) {
  key <- row$dataset_key %||% row$accession
  as.character(key[[1L]])
}

as_dgc <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  as(as(x, "CsparseMatrix"), "dgCMatrix")
}

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
  flush.console()
}

sample_stratified <- function(cells, strata, n) {
  cells <- as.character(cells)
  strata <- as.character(strata)
  strata[is.na(strata) | !nzchar(strata)] <- "Unknown"
  if (length(cells) <= n) return(cells)
  groups <- split(cells, strata)
  quota <- max(1L, floor(n / length(groups)))
  selected <- unlist(lapply(groups, function(z) sample(z, min(length(z), quota))),
                     use.names = FALSE)
  if (length(selected) < n) {
    remaining <- setdiff(cells, selected)
    selected <- c(selected, sample(remaining, min(n - length(selected), length(remaining))))
  }
  sample(selected, n)
}

collapse_rows <- function(x, new_names) {
  x <- as_dgc(x)
  keep <- !is.na(new_names) & nzchar(new_names)
  x <- x[keep, , drop = FALSE]
  new_names <- as.character(new_names[keep])
  target <- unique(new_names)
  map <- Matrix::sparseMatrix(i = match(new_names, target),
                              j = seq_along(new_names), x = 1,
                              dims = c(length(target), length(new_names)))
  answer <- as_dgc(map %*% x)
  dimnames(answer) <- list(target, colnames(x))
  answer
}

expand_rows <- function(x, features) {
  x <- as_dgc(x)
  map <- Matrix::sparseMatrix(i = match(rownames(x), features),
                              j = seq_len(nrow(x)), x = 1,
                              dims = c(length(features), nrow(x)))
  answer <- as_dgc(map %*% x)
  dimnames(answer) <- list(features, colnames(x))
  answer
}

bind_feature_union <- function(matrices) {
  matrices <- Filter(function(x) !is.null(x) && ncol(x) > 0L, matrices)
  if (!length(matrices)) return(NULL)
  features <- unique(unlist(lapply(matrices, rownames), use.names = FALSE))
  answer <- do.call(cbind, lapply(matrices, expand_rows, features = features))
  as_dgc(answer)
}

log_normalize <- function(x, scale_factor = 1e4) {
  x <- as_dgc(x)
  original_dimnames <- dimnames(x)
  lib <- Matrix::colSums(x)
  answer <- as_dgc(x %*% Matrix::Diagonal(x = scale_factor / pmax(lib, 1)))
  answer@x <- log1p(answer@x)
  dimnames(answer) <- original_dimnames
  answer
}

clr_normalize <- function(x) {
  x <- as_dgc(x)
  if (!nrow(x)) return(x)
  dense <- as.matrix(x)
  gm <- exp(colMeans(log1p(dense)))
  answer <- Matrix::Matrix(log1p(sweep(dense, 2, gm, "/")), sparse = TRUE)
  dimnames(answer) <- dimnames(x)
  answer
}

map_annotation_ids <- function(ids, database, keytype) {
  keys <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!length(keys)) return(rep(NA_character_, length(ids)))
  mapped <- suppressMessages(AnnotationDbi::mapIds(
    database, keys = keys, column = "SYMBOL", keytype = keytype,
    multiVals = "first"
  ))
  unname(mapped[ids])
}

canonicalize_human <- function(x) {
  identifiers <- sub("-RNA$", "", rownames(x), ignore.case = TRUE)
  if (is.null(identifiers) || length(identifiers) != nrow(x))
    stop("Gene matrix lacks feature names")
  stripped <- sub("\\.[0-9]+$", "", identifiers)
  target <- toupper(identifiers)
  ensembl <- grepl("^ENSG[0-9]+$", stripped, ignore.case = TRUE)
  entrez <- grepl("^[0-9]+$", stripped)
  if (any(ensembl))
    target[ensembl] <- map_annotation_ids(toupper(stripped[ensembl]),
                                           org.Hs.eg.db::org.Hs.eg.db, "ENSEMBL")
  if (any(entrez))
    target[entrez] <- map_annotation_ids(stripped[entrez],
                                         org.Hs.eg.db::org.Hs.eg.db, "ENTREZID")
  collapse_rows(x, target)
}

map_mouse_to_human <- function(x) {
  identifiers <- sub("-RNA$", "", rownames(x), ignore.case = TRUE)
  if (is.null(identifiers) || length(identifiers) != nrow(x))
    stop("Gene matrix lacks feature names")
  stripped <- sub("\\.[0-9]+$", "", identifiers)
  mouse_symbol <- identifiers
  ensembl <- grepl("^ENSMUSG[0-9]+$", stripped, ignore.case = TRUE)
  entrez <- grepl("^[0-9]+$", stripped)
  if (any(ensembl))
    mouse_symbol[ensembl] <- map_annotation_ids(toupper(stripped[ensembl]),
                                                 org.Mm.eg.db::org.Mm.eg.db,
                                                 "ENSEMBL")
  if (any(entrez))
    mouse_symbol[entrez] <- map_annotation_ids(stripped[entrez],
                                               org.Mm.eg.db::org.Mm.eg.db,
                                               "ENTREZID")
  query <- unique(mouse_symbol[!is.na(mouse_symbol) & nzchar(mouse_symbol)])
  mapped <- tryCatch(
    suppressWarnings(babelgene::orthologs(query, species = "mouse",
                                          human = FALSE, min_support = 2, top = TRUE)),
    error = function(e) NULL
  )
  if (is.null(mapped) || !nrow(mapped)) {
    warning("Mouse ortholog lookup failed; falling back to upper-case symbols")
    return(list(matrix = collapse_rows(x, toupper(mouse_symbol)),
                method = "mouse_symbol_uppercase_fallback"))
  }
  mapped <- mapped[order(-mapped$support_n, mapped$symbol, mapped$human_symbol), ]
  lookup <- setNames(mapped$human_symbol, mapped$symbol)
  lookup <- lookup[!duplicated(names(lookup))]
  target <- unname(lookup[mouse_symbol])
  used <- unique(data.frame(mouse_feature = identifiers, human_symbol = target,
                            stringsAsFactors = FALSE))
  used <- used[!is.na(used$human_symbol) & nzchar(used$human_symbol), ]
  list(matrix = collapse_rows(x, toupper(target)), method = "babelgene",
       map = used)
}

infer_species <- function(source_root, hint = "auto", genes = character()) {
  hint <- tolower(trimws(hint))
  if (hint %in% c("human", "mouse")) return(hint)
  stripped_genes <- toupper(sub("\\.[0-9]+$", "", genes))
  if (any(grepl("^ENSMUSG[0-9]+$", stripped_genes))) return("mouse")
  if (any(grepl("^ENSG[0-9]+$", stripped_genes))) return("human")
  path <- tolower(source_root)
  if (grepl("mouse|小鼠", path)) return("mouse")
  if (grepl("human|重新跑人", path)) return("human")
  series <- list.files(source_root, pattern = "series_matrix.*\\.txt\\.gz$",
                       recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  for (file in head(series, 3L)) {
    organism <- tryCatch({
      con <- gzfile(file, "rt")
      on.exit(close(con), add = TRUE)
      lines <- readLines(con, n = 300L, warn = FALSE)
      paste(lines[grepl("^!Sample_organism_ch1", lines)], collapse = " ")
    }, error = function(e) "")
    if (grepl("mus musculus|mouse", organism, ignore.case = TRUE)) return("mouse")
    if (grepl("homo sapiens|human", organism, ignore.case = TRUE)) return("human")
  }
  genes <- genes[grepl("^[A-Za-z][A-Za-z0-9.-]*$", genes)]
  if (length(genes)) {
    upper_fraction <- mean(genes == toupper(genes))
    if (upper_fraction >= 0.8) return("human")
    if (upper_fraction <= 0.55) return("mouse")
  }
  "human"
}

standardize_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  answer <- rep("Unassigned", length(x))
  answer[grepl("T[_ ]?cell|CD4|CD8|Treg|NK|NKT", x, ignore.case = TRUE)] <- "T_NK"
  answer[grepl("B[_ ]?cell|Plasma", x, ignore.case = TRUE)] <- "B_cell"
  answer[grepl("Myeloid|Mono|Macro|Dendritic|Neutro", x, ignore.case = TRUE)] <- "Myeloid"
  answer[grepl("Endothelial", x, ignore.case = TRUE)] <- "Endothelial"
  answer[grepl("Fibro|Stromal|Pericyte|Smooth", x, ignore.case = TRUE)] <- "Fibroblast"
  answer[grepl("Epithelial|Keratin|Melanocyte|Tumou?r|Malignant", x,
               ignore.case = TRUE)] <- "Epithelial_or_tumor"
  answer
}

marker_label <- function(expr) {
  markers <- list(
    T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY"),
    B_cell = c("MS4A1", "CD79A", "CD74", "CD37"),
    Myeloid = c("LYZ", "LST1", "FCER1G", "CTSS", "TYROBP"),
    Endothelial = c("PECAM1", "VWF", "KDR", "EMCN"),
    Fibroblast = c("COL1A1", "COL1A2", "DCN", "COL3A1"),
    Epithelial_or_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19")
  )
  scores <- vapply(markers, function(g) {
    g <- intersect(g, rownames(expr))
    if (!length(g)) rep(0, ncol(expr)) else Matrix::colMeans(expr[g, , drop = FALSE])
  }, numeric(ncol(expr)))
  if (is.null(dim(scores))) scores <- matrix(scores, ncol = 1L)
  best <- max.col(scores, ties.method = "first")
  answer <- colnames(scores)[best]
  answer[apply(scores, 1L, max) < 0.05] <- "Unassigned"
  answer
}

pick_column <- function(metadata, patterns, fallback = NULL) {
  for (pattern in patterns) {
    hit <- grep(pattern, names(metadata), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[[1L]])
  }
  fallback
}

choose_seurat_cells <- function(object, gene_assay, n) {
  metadata <- object[[]]
  cells <- rownames(metadata)
  feature_column <- paste0("nFeature_", gene_assay)
  count_column <- paste0("nCount_", gene_assay)
  good <- rep(TRUE, nrow(metadata))
  if (feature_column %in% names(metadata))
    good <- good & is.finite(metadata[[feature_column]]) & metadata[[feature_column]] >= 100
  if (count_column %in% names(metadata))
    good <- good & is.finite(metadata[[count_column]]) & metadata[[count_column]] >= 200
  mt <- grep("percent[._]?mt", names(metadata), ignore.case = TRUE, value = TRUE)
  if (length(mt)) good <- good & (is.na(metadata[[mt[[1L]]]]) | metadata[[mt[[1L]]]] < 30)
  candidates <- cells[good]
  if (!length(candidates)) candidates <- cells
  sample_column <- pick_column(metadata, c("^sample_id$", "^sample_name$", "^sample$",
                                            "^orig.ident$"))
  label_column <- pick_column(metadata, c("cell.*type", "annotation", "singler", "cluster"))
  strata <- rep("all", length(candidates))
  if (!is.null(sample_column)) strata <- as.character(metadata[candidates, sample_column])
  if (!is.null(label_column))
    strata <- paste(strata, metadata[candidates, label_column], sep = "__")
  sample_stratified(candidates, strata, n)
}

assay_matrix <- function(object, assay, cells = NULL) {
  aa <- object[[assay]]
  layers <- SeuratObject::Layers(aa)
  count_layers <- grep("^counts($|\\.)", layers, value = TRUE)
  chosen <- if (length(count_layers)) count_layers else grep("^data($|\\.)", layers, value = TRUE)
  if (!length(chosen)) chosen <- layers[[1L]]
  parts <- lapply(chosen, function(layer) {
    value <- tryCatch(SeuratObject::LayerData(aa, layer = layer), error = function(e) NULL)
    if (is.null(value)) return(NULL)
    if (!is.null(cells)) value <- value[, intersect(cells, colnames(value)), drop = FALSE]
    if (!ncol(value)) return(NULL)
    as_dgc(value)
  })
  matrix <- bind_feature_union(parts)
  if (is.null(matrix)) stop("No readable layer in assay ", assay)
  if (!is.null(cells)) matrix <- matrix[, intersect(cells, colnames(matrix)), drop = FALSE]
  list(matrix = matrix, is_counts = any(chosen %in% count_layers), layers = chosen)
}

extract_seurat_coords <- function(object, cells) {
  images <- tryCatch(names(object@images), error = function(e) character())
  for (image in images) {
    value <- tryCatch(Seurat::GetTissueCoordinates(object[[image]]),
                      error = function(e) NULL)
    if (is.null(value) || !nrow(value)) next
    common <- intersect(cells, rownames(value))
    if (!length(common)) next
    x_name <- grep("imagecol|pxl_col|^x$|col", names(value), ignore.case = TRUE, value = TRUE)[1L]
    y_name <- grep("imagerow|pxl_row|^y$|row", names(value), ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(x_name) || is.na(y_name)) next
    answer <- as.matrix(value[common, c(x_name, y_name), drop = FALSE])
    colnames(answer) <- c("x", "y")
    return(answer)
  }
  NULL
}

tfidf_peaks <- function(counts, max_peaks = 1000L) {
  counts <- as_dgc(counts)
  prevalence <- Matrix::rowMeans(counts > 0)
  score <- prevalence * (1 - prevalence)
  selected <- head(order(score, decreasing = TRUE, na.last = NA), max_peaks)
  counts <- counts[selected, , drop = FALSE]
  original <- rownames(counts)
  tf <- counts %*% Matrix::Diagonal(x = 1 / pmax(Matrix::colSums(counts), 1))
  idf <- log1p(ncol(counts) / pmax(Matrix::rowSums(counts > 0), 1))
  answer <- as_dgc(Matrix::Diagonal(x = idf) %*% tf)
  answer@x <- answer@x * 1e4
  colnames(answer) <- colnames(counts)
  rownames(answer) <- paste0("ATACPEAK__", original)
  list(matrix = answer, original_features = original)
}

make_protein_links <- function(protein, genes) {
  cleaned <- sub("^PROTEIN__", "", rownames(protein))
  cleaned <- gsub("(?i)anti[-_ ]?(human|mouse)?[-_ ]?", "", cleaned, perl = TRUE)
  cleaned <- sub("[-_ ](ADT|TotalSeq.*)$", "", cleaned, ignore.case = TRUE)
  cleaned <- toupper(gsub("[^A-Za-z0-9.-].*$", "", cleaned))
  keep <- cleaned %in% genes
  data.frame(protein = rownames(protein)[keep], gene = cleaned[keep],
             weight = rep(1, sum(keep)),
             stringsAsFactors = FALSE)
}

make_peak_links <- function(object, assay, selected_original) {
  ranges <- tryCatch(GenomicRanges::granges(object[[assay]]), error = function(e) NULL)
  annotation <- tryCatch(Signac::Annotation(object), error = function(e) NULL)
  if (is.null(ranges) || is.null(annotation) || !length(annotation)) return(NULL)
  ranges <- ranges[match(selected_original, rownames(object[[assay]]))]
  nearest <- GenomicRanges::distanceToNearest(ranges, annotation, ignore.strand = TRUE)
  genes <- as.character(S4Vectors::mcols(annotation)$gene_name[S4Vectors::subjectHits(nearest)])
  distance <- S4Vectors::mcols(nearest)$distance
  answer <- data.frame(
    peak = paste0("ATACPEAK__", selected_original[S4Vectors::queryHits(nearest)]),
    gene = toupper(genes), distance = distance,
    weight = exp(-distance / 50000), stringsAsFactors = FALSE
  )
  answer[!is.na(answer$gene) & nzchar(answer$gene) & answer$distance <= 100000, ]
}

finalize_gene_matrix <- function(x, is_counts, species) {
  x <- as_dgc(x)
  if (is_counts) x <- log_normalize(x)
  if (species == "mouse") {
    symbols <- rownames(x)
    symbols <- symbols[grepl("^[A-Za-z][A-Za-z0-9.-]*$", symbols)]
    upper_fraction <- if (length(symbols)) mean(symbols == toupper(symbols)) else 0
    contains_mouse_ensembl <- any(grepl("^ENSMUSG", rownames(x), ignore.case = TRUE))
    if (upper_fraction >= 0.8 && !contains_mouse_ensembl) {
      return(list(matrix = canonicalize_human(x),
                  method = "precanonicalized_mouse_symbols"))
    }
    return(map_mouse_to_human(x))
  }
  list(matrix = canonicalize_human(x), method = "human_symbols")
}

seurat_to_standard <- function(object, row, max_cells = 1500L) {
  assays <- SeuratObject::Assays(object)
  gene_assay <- if (row$omics == "Visium" && "Spatial" %in% assays) "Spatial" else
    if ("RNA" %in% assays) "RNA" else if ("ACTIVITY" %in% assays) "ACTIVITY" else
      if ("SCT" %in% assays) "SCT" else NA_character_
  if (is.na(gene_assay)) stop("No RNA, Spatial, ACTIVITY or SCT assay")
  cells <- choose_seurat_cells(object, gene_assay, max_cells)
  metadata <- object[[]][cells, , drop = FALSE]
  coords <- if (row$omics == "Visium") extract_seurat_coords(object, cells) else NULL
  subset_object <- subset(object, cells = cells)
  gene_raw <- assay_matrix(subset_object, gene_assay, cells)
  species <- infer_species(row$source_root, row$species_hint, rownames(gene_raw$matrix))
  gene_result <- finalize_gene_matrix(gene_raw$matrix, gene_raw$is_counts, species)
  gene <- gene_result$matrix
  common_cells <- colnames(gene)

  protein <- protein_links <- NULL
  protein_assay <- grep("ADT|protein|CITE", assays, ignore.case = TRUE, value = TRUE)[1L]
  if (!is.na(protein_assay)) {
    protein_raw <- assay_matrix(subset_object, protein_assay, common_cells)$matrix
    common_cells <- intersect(common_cells, colnames(protein_raw))
    protein <- clr_normalize(protein_raw[, common_cells, drop = FALSE])
    rownames(protein) <- paste0("PROTEIN__", sub("^PROTEIN__", "", rownames(protein)))
  }

  peak <- peak_links <- NULL
  peak_assay <- grep("^peaks$|ATAC", assays, ignore.case = TRUE, value = TRUE)[1L]
  if (!is.na(peak_assay)) {
    peak_raw <- assay_matrix(subset_object, peak_assay, common_cells)$matrix
    common_cells <- intersect(common_cells, colnames(peak_raw))
    peak_result <- tfidf_peaks(peak_raw[, common_cells, drop = FALSE])
    peak <- peak_result$matrix
    peak_links <- make_peak_links(subset_object, peak_assay, peak_result$original_features)
  }

  gene <- gene[, common_cells, drop = FALSE]
  if (!is.null(protein)) protein <- protein[, common_cells, drop = FALSE]
  if (!is.null(peak)) peak <- peak[, common_cells, drop = FALSE]
  if (!is.null(coords)) {
    common_cells <- intersect(common_cells, rownames(coords))
    gene <- gene[, common_cells, drop = FALSE]
    if (!is.null(protein)) protein <- protein[, common_cells, drop = FALSE]
    if (!is.null(peak)) peak <- peak[, common_cells, drop = FALSE]
    coords <- coords[common_cells, , drop = FALSE]
  }
  if (!length(common_cells)) stop("No aligned cells remain")

  sample_column <- pick_column(metadata, c("^sample_id$", "^sample_name$", "^sample$",
                                            "^orig.ident$"))
  label_column <- pick_column(metadata, c("cell.*type", "annotation", "singler", "cluster"))
  sample_value <- if (is.null(sample_column)) rep(row$accession, length(common_cells)) else
    as.character(metadata[common_cells, sample_column])
  labels <- if (is.null(label_column)) rep("Unassigned", length(common_cells)) else
    standardize_label(metadata[common_cells, label_column])
  fallback <- marker_label(gene)
  labels[labels == "Unassigned"] <- fallback[labels == "Unassigned"]

  dataset_key <- row_dataset_key(row)
  cell_id <- paste0(dataset_key, "__", make.unique(common_cells))
  colnames(gene) <- cell_id
  if (!is.null(protein)) {
    colnames(protein) <- cell_id
    protein_links <- make_protein_links(protein, rownames(gene))
  }
  if (!is.null(peak)) colnames(peak) <- cell_id
  if (!is.null(coords)) rownames(coords) <- cell_id
  output_metadata <- data.frame(
    cell_id = cell_id, original_cell_id = common_cells, dataset = dataset_key,
    accession = row$accession,
    disease = row$disease, omics = row$omics, species = species,
    sample = sample_value, broad_celltype = labels, stringsAsFactors = FALSE
  )
  rownames(output_metadata) <- cell_id
  compact(list(gene = gene, protein = protein, protein_links = protein_links,
               peak = peak, peak_links = peak_links, coords = coords,
               metadata = output_metadata,
               provenance = list(reader = "Seurat", gene_assay = gene_assay,
                                 ortholog_method = gene_result$method)))
}

compact <- function(x) x[!vapply(x, is.null, logical(1))]

read_table_auto <- function(path, header = FALSE) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else path
  is_csv <- grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)
  separator <- if (is_csv) "," else "\t"
  quote_character <- if (is_csv) "\"" else ""
  value <- read.table(connection, header = header, sep = separator,
                      stringsAsFactors = FALSE, check.names = FALSE,
                      quote = quote_character, comment.char = "")
  if (is_csv && !header && nrow(value) > 1L) {
    first <- tolower(trimws(as.character(unlist(value[1L, , drop = TRUE]))))
    csv_header <- any(is.na(value[1L, ])) || any(!nzchar(first)) ||
      any(first %in% c("x", "cell_id", "cell", "barcode", "feature", "gene"))
    if (csv_header) value <- value[-1L, , drop = FALSE]
  }
  value
}

find_pair_file <- function(matrix_file, type = c("features", "barcodes")) {
  type <- match.arg(type)
  name <- basename(matrix_file)
  directory <- dirname(matrix_file)
  replacements <- if (type == "features") c("features.tsv.gz", "genes.tsv.gz",
                                              "features.tsv", "genes.tsv",
                                              "features.csv.gz", "genes.csv.gz",
                                              "features.csv", "genes.csv") else
    c("barcodes.tsv.gz", "barcodes.tsv", "barcodes.csv.gz", "barcodes.csv")
  stem <- sub("matrix\\.mtx(\\.gz)?$", "", name, ignore.case = TRUE)
  candidates <- file.path(directory, paste0(stem, replacements))
  candidates <- c(candidates, file.path(directory, replacements))
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) return(existing[[1L]])
  pattern <- if (type == "features") "(features|genes).*\\.(tsv|csv)(\\.gz)?$" else
    "barcodes.*\\.(tsv|csv)(\\.gz)?$"
  found <- list.files(directory, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!length(found)) return(NA_character_)
  found[[1L]]
}

read_mtx_10x <- function(matrix_file, fallback_root = NULL) {
  feature_file <- find_pair_file(matrix_file, "features")
  barcode_file <- find_pair_file(matrix_file, "barcodes")
  if (is.na(feature_file) && !is.null(fallback_root)) {
    shared <- list.files(fallback_root, pattern = "(features|genes).*\\.(tsv|csv)(\\.gz)?$",
                         recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(shared)) feature_file <- shared[[1L]]
  }
  if (is.na(feature_file) || is.na(barcode_file)) stop("Incomplete 10X triplet")
  connection <- if (grepl("\\.gz$", matrix_file, ignore.case = TRUE))
    gzfile(matrix_file, "rt") else matrix_file
  matrix <- as_dgc(Matrix::readMM(connection))
  features <- read_table_auto(feature_file)
  barcode_table <- read_table_auto(barcode_file)
  barcodes <- barcode_table[[if (ncol(barcode_table) >= 2L) ncol(barcode_table) else 1L]]
  if (nrow(matrix) == nrow(features) && ncol(matrix) == length(barcodes)) {
    # already features by cells
  } else if (ncol(matrix) == nrow(features) && nrow(matrix) == length(barcodes)) {
    matrix <- as_dgc(Matrix::t(matrix))
  } else stop("10X matrix dimensions do not match features/barcodes")
  feature_names <- if (ncol(features) >= 2L) features[[2L]] else features[[1L]]
  rownames(matrix) <- make.unique(as.character(feature_names))
  colnames(matrix) <- make.unique(as.character(barcodes))
  feature_type <- if (ncol(features) >= 3L) as.character(features[[3L]]) else
    rep("Gene Expression", nrow(features))
  protein_rows <- grepl("antibody|protein|ADT", feature_type, ignore.case = TRUE)
  gene_rows <- grepl("gene expression|RNA", feature_type, ignore.case = TRUE)
  if (!any(gene_rows)) gene_rows <- !protein_rows
  compact(list(gene = matrix[gene_rows, , drop = FALSE],
               protein = if (any(protein_rows)) matrix[protein_rows, , drop = FALSE] else NULL,
               original_barcodes = colnames(matrix)))
}

read_h5_10x <- function(path) {
  value <- Seurat::Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
  if (inherits(value, "Matrix") || is.matrix(value))
    return(list(gene = as_dgc(value), original_barcodes = colnames(value)))
  names_lower <- tolower(names(value))
  gene_index <- grep("gene expression|rna", names_lower)[1L]
  if (is.na(gene_index)) gene_index <- 1L
  protein_index <- grep("antibody|protein|adt", names_lower)[1L]
  compact(list(gene = as_dgc(value[[gene_index]]),
               protein = if (!is.na(protein_index)) as_dgc(value[[protein_index]]) else NULL,
               original_barcodes = colnames(value[[gene_index]])))
}

read_rds_any <- function(path) {
  if (!grepl("\\.rds\\.gz$", path, ignore.case = TRUE)) return(readRDS(path))
  probe <- gzfile(path, "rb")
  signature <- readBin(probe, "raw", n = 2L)
  close(probe)
  nested_gzip <- length(signature) == 2L &&
    identical(as.integer(signature), c(0x1fL, 0x8bL))
  if (!nested_gzip) return(readRDS(gzfile(path, "rb")))
  outer <- gzfile(path, "rb")
  inner <- gzcon(outer, allowNonCompressed = FALSE, text = FALSE)
  on.exit(close(inner), add = TRUE)
  readRDS(inner)
}

matrix_object_to_raw <- function(object) {
  if (inherits(object, "Matrix") || is.matrix(object)) {
    matrix <- as_dgc(object)
    features <- rownames(matrix)
    protein_rows <- grepl("(-ADT$|^ADT[_-]|antibody|protein)", features, ignore.case = TRUE)
    gene_rows <- !protein_rows
    rownames(matrix) <- sub("-RNA$", "", rownames(matrix), ignore.case = TRUE)
    return(compact(list(gene = matrix[gene_rows, , drop = FALSE],
                        protein = if (any(protein_rows)) matrix[protein_rows, , drop = FALSE] else NULL,
                        original_barcodes = colnames(matrix))))
  }
  if (is.list(object)) {
    collect_matrices <- function(value, prefix = "") {
      if (inherits(value, "Matrix") || is.matrix(value))
        return(setNames(list(value), prefix))
      if (!is.list(value)) return(list())
      answer <- list()
      for (name in names(value)) {
        child <- collect_matrices(value[[name]], paste0(prefix, if (nzchar(prefix)) "$" else "", name))
        answer <- c(answer, child)
      }
      answer
    }
    possible <- collect_matrices(object)
    if (!length(possible)) stop("RDS list contains no matrix")
    names_lower <- tolower(names(possible))
    score <- 20 * grepl("gene expression|rna", names_lower) +
      15 * grepl("umicount", names_lower) + 10 * grepl("exon", names_lower) +
      5 * grepl("\\$all$", names_lower) - 5 * grepl("downsampl", names_lower)
    gene_index <- which.max(score)
    protein_index <- grep("antibody|protein|adt", names_lower)[1L]
    gene <- as_dgc(possible[[gene_index]])
    protein <- if (!is.na(protein_index)) as_dgc(possible[[protein_index]]) else NULL
    return(compact(list(gene = gene, protein = protein,
                        original_barcodes = colnames(gene))))
  }
  stop("Unsupported RDS object class: ", paste(class(object), collapse = ","))
}

read_coordinate_file <- function(path) {
  value <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) NULL)
  if (is.null(value)) return(NULL)
  if (!any(grepl("barcode", names(value), ignore.case = TRUE))) {
    value <- tryCatch(read.csv(path, header = FALSE, stringsAsFactors = FALSE,
                               check.names = FALSE), error = function(e) NULL)
    if (is.null(value) || ncol(value) < 6L) return(NULL)
    names(value)[1:6] <- c("barcode", "in_tissue", "array_row", "array_col",
                           "pxl_row_in_fullres", "pxl_col_in_fullres")
  }
  barcode_name <- grep("barcode", names(value), ignore.case = TRUE, value = TRUE)[1L]
  x_name <- grep("pxl_col|imagecol|^x$", names(value), ignore.case = TRUE, value = TRUE)[1L]
  y_name <- grep("pxl_row|imagerow|^y$", names(value), ignore.case = TRUE, value = TRUE)[1L]
  if (any(is.na(c(barcode_name, x_name, y_name)))) return(NULL)
  answer <- as.matrix(value[, c(x_name, y_name), drop = FALSE])
  storage.mode(answer) <- "double"
  rownames(answer) <- as.character(value[[barcode_name]])
  colnames(answer) <- c("x", "y")
  answer
}

find_coordinates <- function(search_root, barcodes, data_path = search_root) {
  files <- list.files(search_root,
                      pattern = "tissue_positions.*\\.csv(\\.gz)?$|spatial.*coordinates.*\\.csv(\\.gz)?$",
                      recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) return(NULL)
  best <- NULL
  best_score <- -Inf
  data_directory <- normalizePath(dirname(data_path), winslash = "/", mustWork = FALSE)
  for (file in files) {
    value <- read_coordinate_file(file)
    if (is.null(value)) next
    overlap <- sum(barcodes %in% rownames(value))
    if (!overlap) next
    coordinate_directory <- normalizePath(dirname(file), winslash = "/", mustWork = FALSE)
    left <- strsplit(data_directory, "", fixed = TRUE)[[1L]]
    right <- strsplit(coordinate_directory, "", fixed = TRUE)[[1L]]
    common_length <- min(length(left), length(right))
    shared_prefix <- if (common_length) {
      mismatch <- which(left[seq_len(common_length)] != right[seq_len(common_length)])[1L]
      if (is.na(mismatch)) common_length else mismatch - 1L
    } else 0L
    score <- overlap * 1e5 + shared_prefix
    if (score > best_score) {
      best <- value[barcodes[barcodes %in% rownames(value)], , drop = FALSE]
      best_score <- score
    }
  }
  best
}

qc_sample_raw <- function(raw, sample_id, target, spatial = FALSE, search_root = NULL,
                          data_path = NULL) {
  gene <- as_dgc(raw$gene)
  nfeature <- diff(gene@p)
  counts <- Matrix::colSums(gene)
  good <- which(nfeature >= 100 & counts >= 200)
  if (!length(good)) good <- which(nfeature > 0 & counts > 0)
  if (!length(good)) stop("No nonempty cells")
  selected <- sort(sample(good, min(target, length(good))))
  original <- colnames(gene)[selected]
  prefix <- paste0(sample_id, "__")
  cells <- paste0(prefix, make.unique(original))
  gene <- gene[, selected, drop = FALSE]
  colnames(gene) <- cells
  protein <- raw$protein
  if (!is.null(protein)) {
    protein <- as_dgc(protein)[, original, drop = FALSE]
    colnames(protein) <- cells
  }
  coords <- NULL
  if (spatial && !is.null(search_root)) {
    coords <- find_coordinates(search_root, original, data_path)
    if (!is.null(coords)) {
      keep <- original %in% rownames(coords)
      gene <- gene[, keep, drop = FALSE]
      if (!is.null(protein)) protein <- protein[, keep, drop = FALSE]
      original <- original[keep]
      cells <- cells[keep]
      coords <- coords[original, , drop = FALSE]
      rownames(coords) <- cells
    }
  }
  compact(list(gene = gene, protein = protein, coords = coords,
               sample = rep(sample_id, ncol(gene)), original = original,
               nCount = counts[selected][seq_len(ncol(gene))],
               nFeature = nfeature[selected][seq_len(ncol(gene))]))
}

raw_parts_to_standard <- function(parts, row, max_cells = 1500L) {
  parts <- Filter(function(x) !is.null(x) && ncol(x$gene) > 0L, parts)
  if (!length(parts)) stop("No readable matrix parts")
  gene_counts <- bind_feature_union(lapply(parts, `[[`, "gene"))
  if (ncol(gene_counts) > max_cells) {
    all_samples <- unlist(lapply(parts, `[[`, "sample"), use.names = FALSE)
    selected <- sample_stratified(colnames(gene_counts), all_samples, max_cells)
    gene_counts <- gene_counts[, selected, drop = FALSE]
  }
  cells <- colnames(gene_counts)
  sample_value <- unlist(lapply(parts, `[[`, "sample"), use.names = FALSE)
  names(sample_value) <- unlist(lapply(parts, function(x) colnames(x$gene)), use.names = FALSE)
  original <- unlist(lapply(parts, `[[`, "original"), use.names = FALSE)
  names(original) <- names(sample_value)

  species <- infer_species(row$source_root, row$species_hint, rownames(gene_counts))
  gene_result <- finalize_gene_matrix(gene_counts, TRUE, species)
  gene <- gene_result$matrix
  protein_parts <- lapply(parts, `[[`, "protein")
  protein <- if (any(!vapply(protein_parts, is.null, logical(1)))) {
    features <- unique(unlist(lapply(Filter(Negate(is.null), protein_parts), rownames),
                              use.names = FALSE))
    complete <- lapply(seq_along(parts), function(i) {
      value <- protein_parts[[i]]
      if (is.null(value)) {
        Matrix::sparseMatrix(i = integer(), j = integer(),
                             dims = c(length(features), ncol(parts[[i]]$gene)),
                             dimnames = list(features, colnames(parts[[i]]$gene)))
      } else expand_rows(value, features)
    })
    value <- as_dgc(do.call(cbind, complete))[, cells, drop = FALSE]
    value <- clr_normalize(value)
    rownames(value) <- paste0("PROTEIN__", sub("^PROTEIN__", "", rownames(value)))
    value
  } else NULL

  coords <- NULL
  if (row$omics == "Visium") {
    coord_parts <- lapply(parts, `[[`, "coords")
    if (any(vapply(coord_parts, is.null, logical(1))))
      stop("One or more spatial samples lack tissue coordinates")
    coords <- do.call(rbind, coord_parts)
    if (!all(cells %in% rownames(coords))) stop("Spatial coordinates do not cover all spots")
    coords <- coords[cells, , drop = FALSE]
  }
  dataset_key <- row_dataset_key(row)
  cell_id <- paste0(dataset_key, "__", cells)
  colnames(gene) <- cell_id
  if (!is.null(protein)) colnames(protein) <- cell_id
  if (!is.null(coords)) rownames(coords) <- cell_id
  metadata <- data.frame(
    cell_id = cell_id, original_cell_id = unname(original[cells]), dataset = dataset_key,
    accession = row$accession,
    disease = row$disease, omics = row$omics, species = species,
    sample = unname(sample_value[cells]), broad_celltype = marker_label(gene),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- cell_id
  compact(list(gene = gene, protein = protein,
               protein_links = if (!is.null(protein)) make_protein_links(protein, rownames(gene)) else NULL,
               coords = coords, metadata = metadata,
               provenance = list(reader = "matrix", ortholog_method = gene_result$method)))
}

combine_standard_sources <- function(values, row, max_cells = 1500L) {
  values <- Filter(Negate(is.null), values)
  if (!length(values)) stop("No readable RDS sample objects")
  gene <- bind_feature_union(lapply(values, `[[`, "gene"))
  metadata <- do.call(rbind, lapply(values, `[[`, "metadata"))
  rownames(metadata) <- metadata$cell_id
  metadata <- metadata[colnames(gene), , drop = FALSE]
  metadata$dataset <- row_dataset_key(row)
  metadata$accession <- row$accession
  metadata$disease <- row$disease
  metadata$omics <- row$omics

  bind_optional <- function(name) {
    present <- lapply(values, `[[`, name)
    if (!any(!vapply(present, is.null, logical(1)))) return(NULL)
    features <- unique(unlist(lapply(Filter(Negate(is.null), present), rownames),
                              use.names = FALSE))
    complete <- lapply(seq_along(values), function(i) {
      value <- present[[i]]
      cells <- colnames(values[[i]]$gene)
      if (is.null(value)) {
        Matrix::sparseMatrix(i = integer(), j = integer(),
                             dims = c(length(features), length(cells)),
                             dimnames = list(features, cells))
      } else expand_rows(value, features)
    })
    as_dgc(do.call(cbind, complete))[, colnames(gene), drop = FALSE]
  }
  protein <- bind_optional("protein")
  peak <- bind_optional("peak")
  coords <- NULL
  if (row$omics == "Visium") {
    coordinate_parts <- lapply(values, `[[`, "coords")
    if (any(vapply(coordinate_parts, is.null, logical(1))))
      stop("At least one RDS spatial sample lacks coordinates")
    coords <- do.call(rbind, coordinate_parts)[colnames(gene), , drop = FALSE]
  }
  if (ncol(gene) > max_cells) {
    selected <- sample_stratified(colnames(gene),
                                  paste(metadata$sample, metadata$broad_celltype), max_cells)
    gene <- gene[, selected, drop = FALSE]
    metadata <- metadata[selected, , drop = FALSE]
    if (!is.null(protein)) protein <- protein[, selected, drop = FALSE]
    if (!is.null(peak)) peak <- peak[, selected, drop = FALSE]
    if (!is.null(coords)) coords <- coords[selected, , drop = FALSE]
  }
  links <- function(name) {
    answer <- unique(do.call(rbind, Filter(Negate(is.null), lapply(values, `[[`, name))))
    if (is.null(answer) || !nrow(answer)) NULL else answer
  }
  compact(list(gene = gene, protein = protein, protein_links = links("protein_links"),
               peak = peak, peak_links = links("peak_links"), coords = coords,
               metadata = metadata,
               provenance = list(reader = "multiple_RDS_samples",
                                 sample_objects = length(values))))
}

convert_rds_sample_collection <- function(files, row, max_cells) {
  allocation <- max(50L, ceiling(max_cells / min(length(files), 30L)) + 20L)
  values <- list()
  for (file in head(files, 30L)) {
    sample_id <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(file)))
    temporary_row <- row
    temporary_row$dataset_key <- paste0(row_dataset_key(row), "_", sample_id)
    log_message("Reading RDS sample", basename(file))
    value <- tryCatch({
      object <- read_rds_any(file)
      if (inherits(object, "Seurat")) {
        answer <- seurat_to_standard(object, temporary_row, allocation)
      } else {
        raw <- matrix_object_to_raw(object)
        part <- qc_sample_raw(raw, sample_id, allocation,
                              spatial = row$omics == "Visium",
                              search_root = row$source_root, data_path = file)
        answer <- raw_parts_to_standard(list(part), temporary_row, allocation)
      }
      validate_standard_input(answer, row$omics == "Visium")
    }, error = function(e) {
      log_message("RDS sample skipped:", conditionMessage(e)); NULL
    })
    if (!is.null(value)) values[[length(values) + 1L]] <- value
    if (exists("object", inherits = FALSE)) rm(object)
    rm(value); invisible(gc())
  }
  combine_standard_sources(values, row, max_cells)
}

discover_matrix_inputs <- function(root, fallback_root = NULL) {
  mtx <- list.files(root, pattern = "matrix.*\\.mtx(\\.gz)?$|.*_matrix\\.mtx(\\.gz)?$",
                    recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  shared_features <- if (is.null(fallback_root)) character() else
    list.files(fallback_root, pattern = "(features|genes).*\\.tsv(\\.gz)?$",
               recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  mtx <- mtx[vapply(mtx, function(x)
    (!is.na(find_pair_file(x, "features")) || length(shared_features)) &&
      !is.na(find_pair_file(x, "barcodes")), logical(1))]
  h5 <- list.files(root, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE,
                   ignore.case = TRUE)
  h5 <- h5[!grepl("molecule_info|analysis", basename(h5), ignore.case = TRUE)]
  score <- function(x) {
    10 * grepl("filtered", x, ignore.case = TRUE) -
      5 * grepl("raw_feature", x, ignore.case = TRUE) -
      nchar(x) / 10000
  }
  files <- c(mtx, h5)
  files[order(vapply(files, score, numeric(1)), decreasing = TRUE)]
}

infer_sample_id <- function(path, root = NULL) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  components <- strsplit(normalized, "/", fixed = TRUE)[[1L]]
  gsm <- grep("^GSM[0-9]+", components, ignore.case = TRUE, value = TRUE)
  if (length(gsm)) return(sub("^(GSM[0-9]+).*$", "\\1", gsm[[length(gsm)]],
                              ignore.case = TRUE))
  stem <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(path)))
  stem <- sub("[_-]?matrix.*$|[_-]?filtered.*$|[_-]?raw_feature.*$", "", stem,
              ignore.case = TRUE)
  if (!nzchar(stem) || tolower(stem) %in% c("outs", "spatial")) stem <- basename(dirname(path))
  gsub("[^A-Za-z0-9._-]+", "_", stem)
}

convert_matrix_inputs <- function(root, row, max_cells = 1500L) {
  files <- discover_matrix_inputs(root, row$source_root)
  if (!length(files)) stop("No 10X matrix or H5 input found")
  # Keep one preferred representation per biological sample.
  sample_ids <- vapply(files, infer_sample_id, character(1), root = root)
  keys <- if (anyDuplicated(sample_ids)) sample_ids else
    paste(dirname(files), sub("(matrix.*|filtered.*|raw_feature.*)", "", basename(files),
                             ignore.case = TRUE))
  files <- files[!duplicated(keys)]
  sample_ids <- sample_ids[!duplicated(keys)]
  if (length(files) > 30L) {
    files <- files[seq_len(30L)]
    sample_ids <- sample_ids[seq_len(30L)]
  }
  allocation <- max(25L, ceiling(max_cells / min(length(files), 30L)) + 20L)
  parts <- list()
  for (i in seq_along(files)) {
    file <- files[[i]]
    sample_id <- sample_ids[[i]]
    log_message("Reading", basename(file), "sample", sample_id)
    raw <- tryCatch(if (grepl("\\.h5$", file, ignore.case = TRUE)) read_h5_10x(file)
                    else read_mtx_10x(file, row$source_root), error = function(e) {
                      log_message("Skipped", file, "-", conditionMessage(e)); NULL
                    })
    if (is.null(raw)) next
    part <- tryCatch(qc_sample_raw(raw, sample_id, allocation,
                                   spatial = row$omics == "Visium",
                                   search_root = root, data_path = file),
                     error = function(e) {
                       log_message("Sample rejected:", conditionMessage(e)); NULL
                     })
    if (!is.null(part) && (row$omics != "Visium" || !is.null(part$coords)))
      parts[[length(parts) + 1L]] <- part
    rm(raw, part); invisible(gc())
  }
  raw_parts_to_standard(parts, row, max_cells)
}

discover_dense_inputs <- function(root) {
  files <- list.files(root,
                      pattern = "\\.(txt|tsv|csv)(\\.gz)?$",
                      recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  exclude <- grepl("series_matrix|meta[_-]?data|(^|[_-])meta([_-]|\\.)|cell2sample|annotation|barcode|feature|gene[_-]?name|tcr|rename[_-]?cluster",
                   basename(files), ignore.case = TRUE) |
    grepl("tissue_positions|spatial_enrichment", basename(files), ignore.case = TRUE)
  files <- files[!exclude]
  score <- 20 * grepl("dge|expr|count|umi", basename(files), ignore.case = TRUE) +
    10 * grepl("raw", basename(files), ignore.case = TRUE) -
    file.info(files)$size / 1024^4
  files[order(score, decreasing = TRUE)]
}

read_dense_expression <- function(path, target = 200L) {
  header <- data.table::fread(path, nrows = 0L, showProgress = FALSE,
                              data.table = FALSE, check.names = FALSE)
  column_names <- names(header)
  if (length(column_names) < 3L) stop("Dense table has fewer than three columns")
  probe <- data.table::fread(path, select = 1L, nrows = 200L,
                             showProgress = FALSE, data.table = FALSE,
                             check.names = FALSE)[[1L]]
  candidate_features <- column_names[-1L]
  barcode_like <- function(x)
    grepl("(^|[_:-])[ACGTN]{12,}(-[0-9]+)?$", x, ignore.case = TRUE) |
      nchar(x) > 45L
  metadata_name <- grepl(
    "sample|legend|patient|cell.*id|barcode|phenograph|class|cluster|annotation|cell.*type|condition|group|batch",
    candidate_features, ignore.case = TRUE
  )
  feature_like <- grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$|^[A-Za-z][A-Za-z0-9.-]{0,29}$",
                        candidate_features) & !barcode_like(candidate_features) &
    !metadata_name
  probe_cell_like <- grepl("^[0-9]+(\\.[0-9]+)?$", as.character(probe)) |
    barcode_like(as.character(probe)) |
    grepl("[_:-].*[_:-]", as.character(probe))
  cells_by_genes <- mean(feature_like, na.rm = TRUE) >= 0.70 &&
    mean(probe_cell_like, na.rm = TRUE) >= 0.50

  if (cells_by_genes) {
    value <- data.table::fread(path, nrows = target + 50L, showProgress = FALSE,
                               data.table = FALSE, check.names = FALSE)
    cell_column <- grep("cell.*id|barcode", names(value), ignore.case = TRUE)[1L]
    if (is.na(cell_column)) cell_column <- 1L
    cells <- make.unique(as.character(value[[cell_column]]))
    gene_columns <- which(feature_like) + 1L
    numeric_columns <- vapply(value[gene_columns], is.numeric, logical(1))
    gene_columns <- gene_columns[numeric_columns]
    if (length(gene_columns) < 100L)
      stop("Cell-by-gene table has fewer than 100 numeric gene columns")
    matrix <- as.matrix(value[, gene_columns, drop = FALSE])
    storage.mode(matrix) <- "double"
    matrix <- Matrix::Matrix(t(matrix), sparse = TRUE)
    rownames(matrix) <- make.unique(names(value)[gene_columns])
    colnames(matrix) <- cells
  } else if (length(column_names) > 50L) {
    chosen <- sort(sample(seq.int(2L, length(column_names)),
                          min(length(column_names) - 1L, target * 2L)))
    value <- data.table::fread(path, select = c(1L, chosen), showProgress = FALSE,
                               data.table = FALSE, check.names = FALSE)
    genes <- as.character(value[[1L]])
    matrix <- as.matrix(value[, -1L, drop = FALSE])
    storage.mode(matrix) <- "double"
    matrix <- Matrix::Matrix(matrix, sparse = TRUE)
    rownames(matrix) <- make.unique(genes)
    colnames(matrix) <- make.unique(names(value)[-1L])
  } else {
    value <- data.table::fread(path, showProgress = FALSE, data.table = FALSE,
                               check.names = FALSE)
    identifiers <- as.character(value[[1L]])
    numeric <- as.matrix(value[, -1L, drop = FALSE])
    storage.mode(numeric) <- "double"
    matrix <- Matrix::Matrix(numeric, sparse = TRUE)
    rownames(matrix) <- make.unique(identifiers)
    colnames(matrix) <- make.unique(names(value)[-1L])
  }
  list(gene = as_dgc(matrix), original_barcodes = colnames(matrix))
}

convert_dense_inputs <- function(root, row, max_cells = 1500L) {
  files <- discover_dense_inputs(root)
  files <- files[!grepl("(^|[_-])(ab|adt|hto|antibody)([_-]|$)", basename(files),
                        ignore.case = TRUE)]
  if (!length(files)) stop("No dense expression table found")
  sample_ids <- vapply(files, infer_sample_id, character(1), root = root)
  keep <- !duplicated(sample_ids)
  files <- files[keep]
  sample_ids <- sample_ids[keep]
  allocation <- max(25L, ceiling(max_cells / min(length(files), 30L)) + 20L)
  parts <- list()
  for (file in head(files, 30L)) {
    sample_id <- sample_ids[[match(file, files)]]
    log_message("Reading dense expression", basename(file), "sample", sample_id)
    part <- tryCatch({
      raw <- read_dense_expression(file, allocation)
      qc_sample_raw(raw, sample_id, allocation, spatial = row$omics == "Visium",
                    search_root = root, data_path = file)
    }, error = function(e) {
      log_message("Dense table skipped:", conditionMessage(e)); NULL
    })
    if (!is.null(part) && (row$omics != "Visium" || !is.null(part$coords)))
      parts[[length(parts) + 1L]] <- part
    rm(part); invisible(gc())
  }
  raw_parts_to_standard(parts, row, max_cells)
}

validate_standard_input <- function(value, require_spatial = FALSE) {
  if (!is.list(value) || is.null(value$gene) || is.null(value$metadata))
    stop("Output must contain gene and metadata")
  value$gene <- as_dgc(value$gene)
  cells <- colnames(value$gene)
  features <- rownames(value$gene)
  if (!length(cells) || !nrow(value$gene)) stop("Empty gene matrix")
  if (is.null(features) || length(features) != nrow(value$gene) ||
      anyNA(features) || any(!nzchar(features)))
    stop("Gene matrix lacks valid feature names")
  if (nrow(value$gene) > 100000L)
    stop("Gene matrix has more than 100,000 rows and is likely not an expression matrix")
  barcode_features <- grepl("(^|[_:-])[ACGTN]{12,}(-[0-9]+)?$", features,
                            ignore.case = TRUE)
  if (mean(barcode_features) > 0.50)
    stop("Gene matrix appears transposed: most feature names are cell barcodes")
  if (anyDuplicated(cells) || anyDuplicated(rownames(value$gene))) stop("Duplicate matrix names")
  if (any(!is.finite(value$gene@x))) stop("Gene matrix contains non-finite values")
  if (!identical(cells, as.character(value$metadata$cell_id)))
    stop("Metadata cell order mismatch")
  for (name in intersect(c("protein", "peak"), names(value))) {
    value[[name]] <- as_dgc(value[[name]])
    if (!identical(cells, colnames(value[[name]]))) stop(name, " cell order mismatch")
    if (any(!is.finite(value[[name]]@x))) stop(name, " contains non-finite values")
  }
  if (require_spatial) {
    if (is.null(value$coords)) stop("Spatial input lacks coordinates")
    if (!identical(cells, rownames(value$coords))) stop("Coordinate cell order mismatch")
    if (any(!is.finite(value$coords))) stop("Coordinates contain non-finite values")
  }
  value
}

rank_rds_files <- function(root, max_rds_gb = 6, accession = "") {
  files <- if (file.exists(root) && !dir.exists(root)) {
    if (grepl("\\.rds(\\.gz)?$", root, ignore.case = TRUE)) root else character()
  } else {
    list.files(root, pattern = "\\.rds(\\.gz)?$", recursive = TRUE,
               full.names = TRUE, ignore.case = TRUE)
  }
  size <- file.info(files)$size / 1024^3
  files <- files[is.finite(size) & size <= max_rds_gb]
  size <- size[is.finite(size) & size <= max_rds_gb]
  name <- basename(files)
  score <- 100 * grepl("WNN_SingleR_annotated|processed_final|celltype2", name,
                       ignore.case = TRUE) +
    70 * grepl("(^|_)sc1\\.rds|annotated|merged", name, ignore.case = TRUE) +
    50 * grepl("harmony|after_wnn|after_qc", name, ignore.case = TRUE) +
    20 * grepl("matrix|counts", name, ignore.case = TRUE) -
    60 * grepl("plot_data|embedding|metadata", name, ignore.case = TRUE) - size
  if (nzchar(accession)) {
    exact <- grepl(paste0("^", accession, "\\.rds(\\.gz)?$"), name,
                   ignore.case = TRUE)
    root_level <- !grepl("[/\\\\]samples[/\\\\]", files, ignore.case = TRUE)
    score <- score + 95 * exact + 25 * root_level
  }
  files[order(score, decreasing = TRUE)]
}

convert_rds_inputs <- function(root, row, max_cells = 1500L, max_rds_gb = 6) {
  files <- rank_rds_files(root, max_rds_gb, row$accession)
  if (!length(files)) stop("No RDS below memory safety threshold")
  aggregate_candidate <- !grepl("[/\\\\]samples[/\\\\]", files, ignore.case = TRUE) &
    grepl(paste0("^", row$accession, "\\.rds(\\.gz)?$|WNN_SingleR|processed_final|",
                 "merged|harmony|celltype2|(^|_)sc1\\.rds"),
          basename(files), ignore.case = TRUE)
  if (!any(aggregate_candidate) && length(files) > 1L) {
    answer <- convert_rds_sample_collection(files, row, max_cells)
    return(validate_standard_input(answer, row$omics == "Visium"))
  }
  errors <- character()
  for (file in files) {
    log_message("Trying RDS", file)
    value <- tryCatch({
      object <- read_rds_any(file)
      if (inherits(object, "Seurat")) {
        answer <- seurat_to_standard(object, row, max_cells)
      } else {
        raw <- matrix_object_to_raw(object)
        part <- qc_sample_raw(raw, tools::file_path_sans_ext(basename(file)),
                              max_cells, spatial = row$omics == "Visium",
                              search_root = root, data_path = file)
        answer <- raw_parts_to_standard(list(part), row, max_cells)
      }
      rm(object); invisible(gc())
      validate_standard_input(answer, row$omics == "Visium")
    }, error = function(e) e)
    if (!inherits(value, "error")) {
      value$provenance$source_file <- file
      return(value)
    }
    errors <- c(errors, paste(basename(file), conditionMessage(value), sep = ": "))
    log_message("RDS rejected:", tail(errors, 1L))
    invisible(gc())
  }
  stop(paste(c("All RDS candidates failed", errors), collapse = " | "))
}

extract_archives <- function(source_root, work_root) {
  archives <- list.files(source_root, pattern = "\\.(tar|tar\\.gz|tgz|zip)$",
                         recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(archives)) stop("No supported archive found")
  dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(archives)) {
    archive <- archives[[i]]
    destination <- file.path(work_root, sprintf("archive_%03d", i))
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    log_message("Extracting", archive)
    if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
      tryCatch(utils::unzip(archive, exdir = destination),
               error = function(e) log_message("Extraction failed:", conditionMessage(e)))
    } else {
      tryCatch(utils::untar(archive, exdir = destination),
               error = function(e) log_message("Extraction failed:", conditionMessage(e)))
    }
  }
  nested <- list.files(work_root, pattern = "\\.(tar|tar\\.gz|tgz|zip)$",
                       recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  for (i in seq_along(nested)) {
    destination <- file.path(dirname(nested[[i]]), paste0("nested_", i))
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    tryCatch(if (grepl("\\.zip$", nested[[i]], ignore.case = TRUE))
      utils::unzip(nested[[i]], exdir = destination) else
        utils::untar(nested[[i]], exdir = destination), error = function(e) NULL)
  }
  work_root
}

convert_dataset <- function(row, work_root, max_cells = 1500L, max_rds_gb = 6) {
  direct_rds <- tryCatch(convert_rds_inputs(row$source_root, row, max_cells, max_rds_gb),
                         error = function(e) e)
  if (!inherits(direct_rds, "error")) return(direct_rds)
  log_message("Direct RDS route failed:", conditionMessage(direct_rds))
  direct_matrix <- tryCatch(convert_matrix_inputs(row$source_root, row, max_cells),
                            error = function(e) e)
  if (!inherits(direct_matrix, "error")) return(direct_matrix)
  log_message("Direct matrix route failed:", conditionMessage(direct_matrix))
  direct_dense <- tryCatch(convert_dense_inputs(row$source_root, row, max_cells),
                           error = function(e) e)
  if (!inherits(direct_dense, "error")) return(direct_dense)
  log_message("Direct dense-table route failed:", conditionMessage(direct_dense))
  extracted <- extract_archives(row$source_root, work_root)
  archive_rds <- tryCatch(convert_rds_inputs(extracted, row, max_cells, max_rds_gb),
                         error = function(e) e)
  if (!inherits(archive_rds, "error")) return(archive_rds)
  log_message("Archive RDS route failed:", conditionMessage(archive_rds))
  archive_matrix <- tryCatch(convert_matrix_inputs(extracted, row, max_cells),
                             error = function(e) e)
  if (!inherits(archive_matrix, "error")) return(archive_matrix)
  archive_dense <- tryCatch(convert_dense_inputs(extracted, row, max_cells),
                            error = function(e) e)
  if (!inherits(archive_dense, "error")) return(archive_dense)
  stop("All routes failed. RDS: ", conditionMessage(direct_rds),
       "; matrix: ", conditionMessage(direct_matrix),
       "; dense: ", conditionMessage(direct_dense),
       "; archive RDS: ", conditionMessage(archive_rds),
       "; archive matrix: ", conditionMessage(archive_matrix),
       "; archive dense: ", conditionMessage(archive_dense))
}
