args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "G:/multiomics_integration_data"
manifest <- read.csv(file.path(root, "disease_manifests", "batch_manifest.csv"),
                     stringsAsFactors = FALSE)
if (!"dataset_key" %in% names(manifest))
  manifest$dataset_key <- paste(manifest$accession, manifest$omics, sep = "__")
status_file <- file.path(root, "logs", "conversion_status.csv")
status <- if (file.exists(status_file))
  read.csv(status_file, stringsAsFactors = FALSE) else
  data.frame(dataset_key = character(), status = character())
if (nrow(status) && !"dataset_key" %in% names(status))
  status$dataset_key <- paste(status$accession, status$omics, sep = "__")
joined <- merge(manifest, status, by = "dataset_key", all.x = TRUE,
                suffixes = c("", ".run"))
joined$conversion_status <- ifelse(is.na(joined$status), "pending", joined$status)
write.csv(joined, file.path(root, "logs", "dataset_progress_full.csv"), row.names = FALSE)

summary <- aggregate(accession ~ disease + omics + conversion_status, joined, length)
names(summary)[names(summary) == "accession"] <- "datasets"
write.csv(summary, file.path(root, "logs", "dataset_progress_summary.csv"), row.names = FALSE)

usable <- joined[joined$conversion_status == "success", ]
disease_summary <- do.call(rbind, lapply(split(usable, usable$disease), function(x) {
  data.frame(disease = x$disease[[1L]], successful_datasets = nrow(x),
             successful_omics = length(unique(x$omics)),
             integration_ready = length(unique(x$omics)) >= 2L,
             stringsAsFactors = FALSE)
}))
write.csv(disease_summary, file.path(root, "logs", "integration_readiness.csv"),
          row.names = FALSE)

known_failures <- data.frame(
  accession = c("GSE153056", "GSE179733", "GSE207006", "GSE245840",
                "GSE200115", "GSE215858", "GSE246762", "GSE185344"),
  blocking_reason = c(
    "Source gzip and GEO RAW tar are truncated or damaged",
    "Source gzip-wrapped RDS is truncated and cannot be deserialized",
    "Source 10X gzip/Matrix Market files are internally corrupted",
    "Source 10X gzip/Matrix Market files are truncated",
    "Expression matrix is present, but no true spot-coordinate table is supplied",
    "Expression matrix is present, but no true spot-coordinate table is supplied",
    "Expression matrix is present, but no true spot-coordinate table is supplied",
    "Double-gzip Seurat RDS exceeds the 5 GB memory safety limit"
  ),
  recommended_action = c(
    rep("Re-download the complete GEO supplementary files", 4L),
    rep("Obtain the matching tissue_positions/spot-coordinate file", 3L),
    "Run on a higher-memory machine or provide a smaller processed RDS"
  ),
  stringsAsFactors = FALSE
)
joined <- merge(joined, known_failures, by = "accession", all.x = TRUE)
joined$final_status <- joined$conversion_status
joined$final_status[joined$accession %in% known_failures$accession &
                      joined$conversion_status != "success"] <-
  "unavailable_known"
joined$usable_for_integration <- joined$conversion_status == "success"
joined$blocking_reason[is.na(joined$blocking_reason) &
                         joined$conversion_status != "success"] <-
  joined$message[is.na(joined$blocking_reason) &
                   joined$conversion_status != "success"]
joined$recommended_action[is.na(joined$recommended_action) &
                            joined$conversion_status != "success"] <-
  "Inspect the dataset log and repair or replace the source files"
write.csv(joined, file.path(root, "logs", "final_dataset_summary.csv"),
          row.names = FALSE)
write.csv(joined[joined$final_status != "success",
                 c("disease", "accession", "omics", "source_root",
                   "final_status",
                   "blocking_reason", "recommended_action")],
          file.path(root, "logs", "unavailable_datasets.csv"),
          row.names = FALSE)

integration_file <- file.path(root, "logs", "integration_status.csv")
if (file.exists(integration_file)) {
  integration <- read.csv(integration_file, stringsAsFactors = FALSE,
                          check.names = FALSE)
  slugify <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
  metrics <- lapply(integration$disease, function(disease) {
    path <- file.path(root, "results", slugify(disease), "stage_metrics.csv")
    if (!file.exists(path)) return(data.frame())
    value <- read.csv(path, stringsAsFactors = FALSE)
    final <- value[nrow(value), , drop = FALSE]
    data.frame(
      disease = disease,
      final_normalized_mixing = final$normalized_mixing,
      final_neighbor_celltype_consistency = final$neighbor_celltype_consistency,
      metrics_file = path,
      result_file = file.path(root, "results", slugify(disease),
                              "complete_results.rds"),
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, metrics)
  integration <- merge(integration, metrics, by = "disease", all.x = TRUE)
  write.csv(integration, file.path(root, "logs", "final_integration_summary.csv"),
            row.names = FALSE)
}

result_root <- file.path(root, "results")
list_result_files <- function(pattern) {
  if (!dir.exists(result_root)) return(character())
  list.files(result_root, pattern = pattern, recursive = TRUE,
             full.names = TRUE)
}
png_files <- list_result_files("\\.png$")
acceptance <- data.frame(
  metric = c(
    "manifest_datasets", "successful_clean_inputs", "known_unavailable_inputs",
    "target_diseases", "successful_integrations", "integrated_cells",
    "complete_result_rds", "nonempty_png", "umap_3d_png", "umap_3d_plot_data_csv"
  ),
  value = c(
    nrow(joined), sum(joined$final_status == "success"),
    sum(joined$final_status == "unavailable_known"),
    length(unique(joined$disease)),
    if (exists("integration")) sum(integration$status == "success") else 0L,
    if (exists("integration")) sum(integration$cells[integration$status == "success"],
                                    na.rm = TRUE) else 0L,
    length(list_result_files("^complete_results\\.rds$")),
    sum(file.info(png_files)$size > 0L, na.rm = TRUE),
    length(list_result_files("^UMAP_3D_.*\\.png$")),
    length(list_result_files("^UMAP_3D_plot_data\\.csv$"))
  ),
  stringsAsFactors = FALSE
)
write.csv(acceptance, file.path(root, "logs", "final_run_acceptance.csv"),
          row.names = FALSE)
print(table(joined$conversion_status, useNA = "ifany"))
print(table(joined$final_status, useNA = "ifany"))
cat("Integration-ready diseases:", sum(disease_summary$integration_ready), "\n")
