.libPaths(c("D:/R/R-4.3.1/library/seuratV5", .libPaths()))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
manifest_file <- if (length(args) >= 1L) args[[1L]] else
  "G:/multiomics_integration_data_corrected/disease_manifests/batch_manifest.csv"
destination_root <- if (length(args) >= 2L) args[[2L]] else
  "G:/multiomics_integration_data_corrected"

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_argument)) sub("^--file=", "", file_argument[[1L]]) else
  sys.frame(1)$ofile
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "multiomics_io.R"), encoding = "UTF-8")

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE,
                     fileEncoding = "UTF-8-BOM")
if (!"dataset_key" %in% names(manifest))
  manifest$dataset_key <- paste(manifest$accession, manifest$omics, sep = "__")
if (!"reuse_file" %in% names(manifest)) manifest$reuse_file <- ""
manifest <- manifest[nzchar(manifest$reuse_file) & file.exists(manifest$reuse_file), ,
                     drop = FALSE]

rename_standard_cells <- function(value, row) {
  old_cells <- colnames(value$gene)
  metadata <- value$metadata[old_cells, , drop = FALSE]
  original <- if ("original_cell_id" %in% names(metadata))
    as.character(metadata$original_cell_id) else old_cells
  original[is.na(original) | !nzchar(original)] <- old_cells[is.na(original) |
                                                               !nzchar(original)]
  new_cells <- paste0(row$dataset_key, "__", make.unique(original))
  colnames(value$gene) <- new_cells
  if (!is.null(value$protein)) colnames(value$protein) <- new_cells
  if (!is.null(value$peak)) colnames(value$peak) <- new_cells
  if (!is.null(value$coords)) rownames(value$coords) <- new_cells
  metadata$cell_id <- new_cells
  metadata$original_cell_id <- original
  metadata$dataset <- row$dataset_key
  metadata$accession <- row$accession
  metadata$disease <- row$disease
  metadata$omics <- row$omics
  rownames(metadata) <- new_cells
  value$metadata <- metadata
  value$provenance$original_reader <- value$provenance$reader %||% NA_character_
  value$provenance$reader <- "reused_clean_input"
  value$provenance$reused_from <- row$reuse_file
  value$provenance$dataset_key <- row$dataset_key
  value$provenance$accession <- row$accession
  value$provenance$disease <- row$disease
  value$provenance$omics <- row$omics
  value
}

records <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  row <- as.list(manifest[i, , drop = FALSE])
  dir.create(dirname(row$output_file), recursive = TRUE, showWarnings = FALSE)
  status <- tryCatch({
    if (file.exists(row$output_file)) {
      existing <- validate_standard_input(readRDS(row$output_file), row$omics == "Visium")
      if (all(existing$metadata$dataset == row$dataset_key) &&
          all(existing$metadata$omics == row$omics)) {
        rm(existing)
        "already_current"
      } else stop("existing output uses an obsolete dataset key")
    } else {
      value <- validate_standard_input(readRDS(row$reuse_file), row$omics == "Visium")
      value <- rename_standard_cells(value, row)
      value <- validate_standard_input(value, row$omics == "Visium")
      temporary <- paste0(row$output_file, ".partial")
      saveRDS(value, temporary, compress = TRUE)
      validate_standard_input(readRDS(temporary), row$omics == "Visium")
      if (!file.rename(temporary, row$output_file)) stop("could not finalize reused input")
      rm(value)
      "reused"
    }
  }, error = function(e) paste0("failed: ", conditionMessage(e)))
  records[[i]] <- data.frame(dataset_key = row$dataset_key,
                             accession = row$accession, omics = row$omics,
                             status = status, source = row$reuse_file,
                             output_file = row$output_file,
                             stringsAsFactors = FALSE)
  cat(sprintf("[%d/%d] %s %s\n", i, nrow(manifest), row$dataset_key, status))
  invisible(gc())
}

records <- if (length(records)) do.call(rbind, records) else data.frame()
dir.create(file.path(destination_root, "logs"), recursive = TRUE, showWarnings = FALSE)
write.csv(records, file.path(destination_root, "logs", "reused_input_status.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Reusable clean inputs processed:", nrow(records), "\n")
