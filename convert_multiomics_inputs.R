.libPaths(c("D:/R/R-4.3.1/library/seuratV5", .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(babelgene)
  library(GenomicRanges)
})
set.seed(2026)

args <- commandArgs(trailingOnly = TRUE)
manifest_file <- if (length(args) >= 1L) args[[1L]] else
  "G:/multiomics_integration_data/disease_manifests/batch_manifest.csv"
destination_root <- if (length(args) >= 2L) args[[2L]] else
  "G:/multiomics_integration_data"
mode <- if (length(args) >= 3L) tolower(args[[3L]]) else "all"
max_cells <- as.integer(Sys.getenv("MAX_CELLS_PER_DATASET", "1500"))
max_rds_gb <- as.numeric(Sys.getenv("MAX_RDS_GB", "6"))
selected_accessions <- strsplit(Sys.getenv("ACCESSIONS", ""), ",", fixed = TRUE)[[1L]]
selected_accessions <- trimws(selected_accessions[nzchar(selected_accessions)])
force_accessions <- strsplit(Sys.getenv("FORCE_ACCESSIONS", ""), ",", fixed = TRUE)[[1L]]
force_accessions <- trimws(force_accessions[nzchar(force_accessions)])

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_argument)) sub("^--file=", "", file_argument[[1L]]) else
  sys.frame(1)$ofile
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "multiomics_io.R"), encoding = "UTF-8")

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE,
                     fileEncoding = "UTF-8-BOM")
if (!"dataset_key" %in% names(manifest))
  manifest$dataset_key <- paste(manifest$accession, manifest$omics, sep = "__")
if (mode %in% c("ready", "convertible"))
  manifest <- manifest[manifest$source_status == mode, , drop = FALSE]
if (length(selected_accessions))
  manifest <- manifest[manifest$accession %in% selected_accessions, , drop = FALSE]
manifest <- manifest[order(factor(manifest$source_status, c("ready", "convertible")),
                           manifest$omics, manifest$accession), , drop = FALSE]

directories <- c(file.path(destination_root, "clean_inputs"),
                 file.path(destination_root, "logs", "datasets"),
                 file.path(destination_root, "work"))
invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))
status_file <- file.path(destination_root, "logs", "conversion_status.csv")
status <- if (file.exists(status_file))
  read.csv(status_file, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-BOM") else data.frame()

record_status <- function(record) {
  if (nrow(status)) {
    if (!"dataset_key" %in% names(status))
      status$dataset_key <<- paste(status$accession, status$omics, sep = "__")
    status <<- status[status$dataset_key != record$dataset_key, , drop = FALSE]
  }
  status <<- rbind(status, record)
  status <<- status[order(status$dataset_key), , drop = FALSE]
  write.csv(status, status_file, row.names = FALSE, fileEncoding = "UTF-8")
}

for (index in seq_len(nrow(manifest))) {
  row <- as.list(manifest[index, , drop = FALSE])
  output_file <- row$output_file
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  log_file <- row$log_file
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()

  existing_ok <- FALSE
  existing <- NULL
  if (file.exists(output_file) &&
      !row$accession %in% force_accessions && !row$dataset_key %in% force_accessions) {
    existing_ok <- tryCatch({
      existing <- validate_standard_input(readRDS(output_file), row$omics == "Visium")
      TRUE
    }, error = function(e) FALSE)
  }
  if (existing_ok) {
    log_message(sprintf("[%d/%d] %s already complete", index, nrow(manifest), row$accession))
    record_status(data.frame(
      dataset_key = row$dataset_key, accession = row$accession,
      disease = row$disease, omics = row$omics,
      status = "success", cells = ncol(existing$gene), genes = nrow(existing$gene),
      proteins = if (is.null(existing$protein)) 0L else nrow(existing$protein),
      peaks = if (is.null(existing$peak)) 0L else nrow(existing$peak),
      has_coordinates = !is.null(existing$coords),
      species = paste(unique(existing$metadata$species), collapse = ";"),
      output_file = output_file,
      source_reader = existing$provenance$reader %||% NA_character_,
      elapsed_minutes = 0,
      message = "validated existing clean input",
      updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      stringsAsFactors = FALSE
    ))
    rm(existing); invisible(gc())
    next
  }

  log_connection <- file(log_file, open = "at", encoding = "UTF-8")
  sink(log_connection, split = TRUE)
  work_root <- file.path(destination_root, "work", row$dataset_key)
  result <- tryCatch({
    log_message(sprintf("[%d/%d] START %s %s %s", index, nrow(manifest),
                        row$accession, row$omics, row$disease))
    log_message("Source:", row$source_root)
    value <- convert_dataset(row, work_root, max_cells = max_cells,
                             max_rds_gb = max_rds_gb)
    value <- validate_standard_input(value, row$omics == "Visium")
    value$provenance$accession <- row$accession
    value$provenance$dataset_key <- row$dataset_key
    value$provenance$disease <- row$disease
    value$provenance$omics <- row$omics
    value$provenance$source_root <- row$source_root
    value$provenance$created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
    value$provenance$max_cells <- max_cells
    temporary <- paste0(output_file, ".partial")
    saveRDS(value, temporary, compress = TRUE)
    check <- validate_standard_input(readRDS(temporary), row$omics == "Visium")
    if (file.exists(output_file)) unlink(output_file)
    if (!file.rename(temporary, output_file)) stop("Could not finalize output file")
    record <- data.frame(
      dataset_key = row$dataset_key, accession = row$accession,
      disease = row$disease, omics = row$omics,
      status = "success", cells = ncol(check$gene), genes = nrow(check$gene),
      proteins = if (is.null(check$protein)) 0L else nrow(check$protein),
      peaks = if (is.null(check$peak)) 0L else nrow(check$peak),
      has_coordinates = !is.null(check$coords),
      species = paste(unique(check$metadata$species), collapse = ";"),
      output_file = output_file,
      source_reader = check$provenance$reader %||% NA_character_,
      elapsed_minutes = round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2),
      message = "", updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      stringsAsFactors = FALSE
    )
    log_message("SUCCESS", row$accession, ncol(check$gene), "cells",
                nrow(check$gene), "genes")
    rm(value, check); invisible(gc())
    record
  }, error = function(e) {
    log_message("FAILED", row$accession, conditionMessage(e))
    data.frame(
      dataset_key = row$dataset_key, accession = row$accession,
      disease = row$disease, omics = row$omics,
      status = "failed", cells = NA_integer_, genes = NA_integer_,
      proteins = NA_integer_, peaks = NA_integer_, has_coordinates = NA,
      species = NA_character_, output_file = output_file,
      source_reader = NA_character_,
      elapsed_minutes = round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2),
      message = conditionMessage(e),
      updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      stringsAsFactors = FALSE
    )
  })
  sink()
  close(log_connection)
  record_status(result)
  if (dir.exists(work_root)) unlink(work_root, recursive = TRUE, force = TRUE)
  invisible(gc(full = TRUE))
}

success <- sum(status$status == "success" & status$dataset_key %in% manifest$dataset_key)
failed <- sum(status$status == "failed" & status$dataset_key %in% manifest$dataset_key)
cat(sprintf("Conversion pass complete: %d success, %d failed, %d selected datasets\n",
            success, failed, nrow(manifest)))
