.libPaths(c("D:/R/R-4.3.1/library/seuratV5", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
report_dir <- if (length(args) >= 1L) args[[1L]] else
  "F:/老师/数据库/审稿意见/diseasedetail_copy_report"
metadata_file <- if (length(args) >= 2L) args[[2L]] else
  "F:/xwechat_files/wxid_4kntt2elglj422_962c/msg/file/2026-09/dataset_tools_new2(1)(1).csv"
destination_root <- if (length(args) >= 3L) args[[3L]] else
  "G:/multiomics_integration_data"
output_file <- if (length(args) >= 4L) args[[4L]] else
  file.path(report_dir, "batch_manifest.csv")

read_utf8_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-BOM")
}

details <- read_utf8_csv(file.path(report_dir, "integration_dataset_details.csv"))
catalogue <- read_utf8_csv(metadata_file)

species_by_accession <- split(as.character(catalogue$species), catalogue$dataset_id)
get_species <- function(accession, source_root) {
  values <- unique(species_by_accession[[accession]])
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 1L) return(tolower(values))
  path <- tolower(source_root)
  if (grepl("mouse|小鼠", path)) return("mouse")
  if (grepl("human|人类|重新跑人", path)) return("human")
  "auto"
}

status_map <- c("可直接读取" = "ready", "需先解压或转换" = "convertible")
manifest <- data.frame(
  disease = details[["疾病"]],
  accession = details[["GSE号"]],
  omics = details[["组学"]],
  source_status = unname(status_map[details[["处理状态"]]]),
  source_reason = details[["输入依据"]],
  source_root = normalizePath(details[["I盘源目录"]], winslash = "/",
                              mustWork = TRUE),
  species_hint = mapply(get_species, details[["GSE号"]], details[["I盘源目录"]],
                        USE.NAMES = FALSE),
  stringsAsFactors = FALSE
)
manifest$output_file <- file.path(destination_root, "clean_inputs", manifest$omics,
                                 paste0(manifest$accession, ".rds"))
manifest$log_file <- file.path(destination_root, "logs", "datasets",
                               paste0(manifest$accession, ".log"))
manifest <- manifest[order(manifest$disease, manifest$omics, manifest$accession), ]
if (anyDuplicated(manifest$accession))
  stop("The integration manifest must contain one row per accession")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, output_file, row.names = FALSE, fileEncoding = "UTF-8")
disease_directory <- file.path(dirname(output_file), "diseases")
dir.create(disease_directory, recursive = TRUE, showWarnings = FALSE)
slugify <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
for (disease in unique(manifest$disease)) {
  write.csv(manifest[manifest$disease == disease, , drop = FALSE],
            file.path(disease_directory, paste0(slugify(disease), ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}
cat(sprintf("Wrote %d datasets from %d diseases to %s\n",
            nrow(manifest), length(unique(manifest$disease)), output_file))
print(with(manifest, table(omics, source_status)))
