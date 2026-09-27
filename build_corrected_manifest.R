.libPaths(c("D:/R/R-4.3.1/library/seuratV5", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
old_manifest_file <- if (length(args) >= 1L) args[[1L]] else
  "F:/老师/数据库/审稿意见/diseasedetail_copy_report/batch_manifest.csv"
disease_file <- if (length(args) >= 2L) args[[2L]] else
  "F:/xwechat_files/wxid_4kntt2elglj422_962c/msg/file/2026-09/diseasedetail.csv"
catalogue_file <- if (length(args) >= 3L) args[[3L]] else
  "F:/xwechat_files/wxid_4kntt2elglj422_962c/msg/file/2026-09/dataset_tools_new2(1)(1).csv"
destination_root <- if (length(args) >= 4L) args[[4L]] else
  "G:/multiomics_integration_data_corrected"
report_dir <- if (length(args) >= 5L) args[[5L]] else
  file.path(destination_root, "disease_manifests")
disease_scope_file <- if (length(args) >= 6L) args[[6L]] else ""

read_utf8_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-BOM")
}

canonical_dataset_id <- function(x) {
  value <- trimws(as.character(x))
  has_gse <- grepl("GSE[0-9]+", value, ignore.case = TRUE, perl = TRUE)
  value[has_gse] <- sub("^.*?(GSE[0-9]+).*$", "\\1",
                        toupper(value[has_gse]), perl = TRUE)
  value
}

normalize_existing <- function(path) {
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

old <- read_utf8_csv(old_manifest_file)
disease_detail <- read_utf8_csv(disease_file)
catalogue <- read_utf8_csv(catalogue_file)
scoped_diseases <- character()
if (nzchar(disease_scope_file) && file.exists(disease_scope_file)) {
  scope <- read_utf8_csv(disease_scope_file)
  if (!"disease" %in% names(scope))
    stop("Disease scope file must contain a disease column")
  scoped_diseases <- sort(unique(trimws(as.character(scope$disease))))
  scoped_diseases <- scoped_diseases[nzchar(scoped_diseases)]
}

disease_detail$dataset_id <- trimws(as.character(disease_detail$dataset_id))
disease_detail$disease <- trimws(as.character(disease_detail$disease))
disease_detail <- disease_detail[nzchar(disease_detail$dataset_id) &
                                   nzchar(disease_detail$disease), , drop = FALSE]
disease_detail$lookup_id <- canonical_dataset_id(disease_detail$dataset_id)
disease_exact_split <- split(disease_detail$disease, disease_detail$dataset_id)
disease_exact_map <- vapply(disease_exact_split,
                            function(x) sort(unique(x))[[1L]], character(1))
disease_fallback_split <- split(disease_detail$disease, disease_detail$lookup_id)
disease_fallback_map <- vapply(disease_fallback_split,
                               function(x) sort(unique(x))[[1L]], character(1))

catalogue$dataset_id <- trimws(as.character(catalogue$dataset_id))
catalogue$lookup_id <- canonical_dataset_id(catalogue$dataset_id)
summarize_species <- function(x) {
  value <- unique(tolower(trimws(x[!is.na(x) & nzchar(trimws(x))])))
  if (length(value) == 1L && value %in% c("human", "mouse")) value else "auto"
}
species_exact_split <- split(as.character(catalogue$species), catalogue$dataset_id)
species_exact_map <- vapply(species_exact_split, summarize_species, character(1))
species_fallback_split <- split(as.character(catalogue$species), catalogue$lookup_id)
species_fallback_map <- vapply(species_fallback_split, summarize_species, character(1))

lookup_disease <- function(accession) {
  value <- unname(disease_exact_map[accession])
  missing <- is.na(value)
  value[missing] <- unname(disease_fallback_map[
    canonical_dataset_id(accession[missing])
  ])
  value[is.na(value)] <- ""
  value
}

lookup_species <- function(accession) {
  value <- unname(species_exact_map[accession])
  missing <- is.na(value)
  value[missing] <- unname(species_fallback_map[
    canonical_dataset_id(accession[missing])
  ])
  value[is.na(value)] <- "auto"
  value
}

make_rows <- function(disease, accession, omics, source_status, source_reason,
                      source_root, species_hint, source_group) {
  data.frame(
    disease = as.character(disease),
    accession = as.character(accession),
    omics = as.character(omics),
    source_status = as.character(source_status),
    source_reason = as.character(source_reason),
    source_root = as.character(source_root),
    species_hint = as.character(species_hint),
    source_group = as.character(source_group),
    stringsAsFactors = FALSE
  )
}

# These modalities keep the exact source rows used by the previous audited run.
unchanged <- old[old$omics %in% c("Epigenomics", "Proteomics", "Visium"), , drop = FALSE]
unchanged_rows <- make_rows(
  unchanged$disease, unchanged$accession, unchanged$omics,
  unchanged$source_status, unchanged$source_reason, unchanged$source_root,
  unchanged$species_hint, paste0("unchanged_", unchanged$omics)
)

# The 39 previously selected I:/data2 rows are transcriptomics, not genomics.
data2 <- old[grepl("^I:/数据2/", old$source_root), , drop = FALSE]
if (nrow(data2) != 39L)
  stop("Expected 39 selected I:/数据2 datasets, found ", nrow(data2))
data2_rows <- make_rows(
  data2$disease, data2$accession, "Transcriptomics",
  data2$source_status, data2$source_reason, data2$source_root,
  data2$species_hint, "I_data2_selected_39"
)

sc_files <- sort(list.files("G:/SC", pattern = "^GSE[0-9]+\\.rds$",
                            full.names = TRUE, ignore.case = TRUE))
sc_accessions <- tools::file_path_sans_ext(basename(sc_files))
sc_rows <- make_rows(
  lookup_disease(sc_accessions), sc_accessions, "Transcriptomics",
  "ready", "Direct Seurat RDS in corrected G:/SC source",
  vapply(sc_files, normalize_existing, character(1)),
  lookup_species(sc_accessions), "G_SC"
)

genomics_dirs <- sort(list.dirs("I:/基因组/基因组数据", recursive = FALSE,
                                full.names = TRUE))
genomics_accessions <- basename(genomics_dirs)
genomics_rows <- make_rows(
  lookup_disease(genomics_accessions), genomics_accessions, "Genomics",
  "convertible", "Raw GEO archive or matrix in corrected genomics source",
  vapply(genomics_dirs, normalize_existing, character(1)),
  lookup_species(genomics_accessions), "I_genomics_corrected"
)

inventory <- rbind(unchanged_rows, data2_rows, sc_rows, genomics_rows)
inventory$dataset_key <- paste(inventory$accession, inventory$omics, sep = "__")
old_clean_omics <- ifelse(inventory$source_group == "I_data2_selected_39",
                          "Genomics", inventory$omics)
reuse_candidates <- file.path("G:/multiomics_integration_data/clean_inputs",
                              old_clean_omics, paste0(inventory$accession, ".rds"))
inventory$reuse_file <- rep("", nrow(inventory))
reuse_allowed <- inventory$source_group == "I_data2_selected_39" |
  grepl("^unchanged_", inventory$source_group)
reuse_exists <- reuse_allowed & file.exists(reuse_candidates)
inventory$reuse_file[reuse_exists] <- normalizePath(
  reuse_candidates[reuse_exists], winslash = "/", mustWork = TRUE
)
inventory$mapping_status <- ifelse(nzchar(inventory$disease),
                                   "disease_mapped", "disease_not_found")

mapped <- inventory[nzchar(inventory$disease), , drop = FALSE]
if (anyDuplicated(mapped$dataset_key)) {
  duplicate_keys <- unique(mapped$dataset_key[duplicated(mapped$dataset_key)])
  stop("Duplicate dataset_key values: ", paste(duplicate_keys, collapse = ", "))
}

omics_by_disease <- lapply(split(mapped$omics, mapped$disease), unique)
eligible_diseases <- names(Filter(function(x) length(x) >= 2L, omics_by_disease))
if (length(scoped_diseases))
  eligible_diseases <- intersect(eligible_diseases, scoped_diseases)
manifest <- mapped[mapped$disease %in% eligible_diseases, , drop = FALSE]
manifest <- manifest[order(manifest$disease, manifest$omics, manifest$accession), ]
manifest$output_file <- file.path(destination_root, "clean_inputs", manifest$omics,
                                  paste0(manifest$accession, ".rds"))
manifest$log_file <- file.path(destination_root, "logs", "datasets",
                               paste0(manifest$dataset_key, ".log"))

inventory$eligible_multiomics_disease <- inventory$disease %in% eligible_diseases
inventory$included_in_manifest <- inventory$dataset_key %in% manifest$dataset_key

disease_summary <- do.call(rbind, lapply(split(manifest, manifest$disease), function(x) {
  data.frame(
    disease = x$disease[[1L]],
    datasets = nrow(x),
    omics_count = length(unique(x$omics)),
    omics = paste(sort(unique(x$omics)), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
disease_summary <- disease_summary[order(disease_summary$disease), , drop = FALSE]

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, file.path(report_dir, "batch_manifest.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(inventory, file.path(report_dir, "corrected_source_inventory.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(inventory[inventory$mapping_status == "disease_not_found", , drop = FALSE],
          file.path(report_dir, "source_datasets_without_disease_mapping.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(inventory[nzchar(inventory$disease) &
                      !inventory$eligible_multiomics_disease, , drop = FALSE],
          file.path(report_dir, "mapped_but_single_omics_disease.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(disease_summary, file.path(report_dir, "corrected_disease_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("Physical source rows:", nrow(inventory), "\n")
cat("Disease-mapped candidates:", nrow(mapped), "\n")
cat("Eligible diseases:", length(eligible_diseases), "\n")
cat("Final manifest datasets:", nrow(manifest), "\n")
print(with(manifest, table(omics, source_group)))
