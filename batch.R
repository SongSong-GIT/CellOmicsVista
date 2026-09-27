# -*- coding: UTF-8 -*-
# ============================================================
# 单细胞数据按 disease_group 内部/外部组合与 Seurat V4 标准流程脚本
# 功能：
#   1. 匹配 input_dir 下 RDS 与样本信息表中的 dataset_id。
#   2. 读取每个 RDS，记录实际细胞数/基因数。
#   3. 所有疾病分组逻辑统一使用样本表 disease_group 列，不再使用 disease 列。
#   4. disease_in：
#      - Human / Mouse 分开；
#      - 在同一个 disease_group 内组织 dataset；
#      - 1 个 dataset 时也生成 1 个组合并执行分析；
#      - 2-3 个 dataset 时全部直接合并；
#      - >=4 个 dataset 时，使用 combn(..., 4) 生成所有不重复的 4-dataset 组合。
#   5. disease_out：
#      - Human / Mouse 分开；
#      - 按 disease_group 类型进行组合；
#      - 可用 disease_group 少于 6 种时，全部可用 dataset 直接合并；
#      - 可用 disease_group >=6 种时，构造 6-disease_group 组合空间，
#        固定随机种子后随机选择最多 20 个不重复组合运行；
#      - 每个 6-disease_group 组合合并这 6 个 disease_group 下的全部 dataset。
#   6. 每个最终组合最多保留 200,000 个细胞：
#      - 总细胞数 <=200,000 时不下采样；
#      - 总细胞数 >200,000 时，按各 dataset 原始细胞数比例分配配额，
#        单个 RDS 先下采样，再逐个 merge，降低内存峰值。
#   7. 合并/读取后执行 NormalizeData -> FindVariableFeatures -> ScaleData -> PCA ->
#      Harmony（可用时）-> 3D UMAP -> FindNeighbors -> FindClusters(resolution=0.1)。
#   8. 单 dataset disease_in 若 sample_id 唯一值 >=2，仍执行 Harmony；
#      若 sample_id 唯一值 <2，则自动使用 PCA 继续 UMAP/邻居/聚类，不因 Harmony 条件不足而停止。
#   9. 保存组合索引、组合文本记录、RDS、metadata CSV、3D UMAP PNG 和运行汇总。
#  10. human_disease_in / mouse_disease_in 额外保存 in_combination_diagnostics.csv，
#      其中单 dataset 的 expected_combination_count 记为 1。
#
# 说明：
#   - disease_in 的“同一疾病”定义完全由 disease_group 相同决定。
#   - disease_out 的“不同疾病类型”同样完全由 disease_group 决定。
#   - “排列组合”按无顺序组合处理，同一集合的不同排列不会重复生成。
#   - disease_out 随机最多 20 个组合使用固定随机种子，可重复复现。
#   - Human 与 Mouse 始终分别处理。
#   - 所有 UMAP 可视化均使用 UMAP_1 / UMAP_2 / UMAP_3 三维静态 PNG 展示。
#   - disease_group 仅 trim 首尾空格，不进行语义合并或大小写归一化。
# ============================================================

options(stringsAsFactors = FALSE)

# ============================================================
# 1. R 包检查与加载
# ============================================================

# 功能：检查并加载所需 R 包；若缺失则停止并给出安装提示。
check_and_load_packages <- function(pkgs){
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing_pkgs) > 0){
    stop(
      paste0(
        "[ERROR] 缺少 R 包：", paste(missing_pkgs, collapse = ", "),
        "\n请先安装，例如：install.packages(c('",
        paste(missing_pkgs, collapse = "','"), "'))"
      ),
      call. = FALSE
    )
  }
  for(pkg in pkgs){
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  cat("[CHECK] R 包已加载：", paste(pkgs, collapse = ", "), "\n")
}

required_pkgs <- c("Seurat", "SeuratObject", "harmony", "dplyr", "ggplot2", "Matrix", "scatterplot3d")
check_and_load_packages(required_pkgs)

seurat_version <- as.character(packageVersion("Seurat"))
cat("[CHECK] Seurat 版本：", seurat_version, "\n")
if(!grepl("^4\\.", seurat_version)){
  stop(
    paste0("[ERROR] 本脚本按 Seurat V4 编写，当前版本为 ", seurat_version,
           "。请在 Seurat V4 环境中运行。"),
    call. = FALSE
  )
}

# ============================================================
# 2. 路径与基础参数
# ============================================================

input_dir <- "/bioccNSW/donglonglong/donglonglong/DATABASE_SC_HAMARY/input/"
metadata_file <- "/bioccNSW/donglonglong/donglonglong/DATABASE_SC_HAMARY/disease_classified.csv"
output_dir <- file.path("/bioccNSW/donglonglong/donglonglong/DATABASE_SC_HAMARY/output/")

# 组合规则参数。
# disease_in：同一 disease_group 内，1 个 dataset 也产生组合；2-3 个全并；>=4 个时每 4 个组成一个组合。
in_combination_dataset_size <- 4L

# disease_out：按 6 个不同 disease_group 组成一个组合；随机运行最多 20 个组合。
out_combination_disease_size <- 6L
out_random_combination_n <- 20L

# 每个最终合并对象的细胞数上限。
max_merged_cells <- 200000L

# 当 disease 组合总空间较小时可直接枚举；超过该值时改用随机唯一抽样，
# 避免先生成海量 combn 结果占用内存。
combination_enumeration_threshold <- 100000L

# Seurat/Harmony 参数。
npcs <- 50L
dims_use <- 1:15
cluster_resolution <- 0.1
harmony_group_col <- "sample_id"

# 3D UMAP 静态 PNG 参数。
umap3d_angle <- 55
umap3d_png_width <- 1800L
umap3d_png_height <- 1400L
umap3d_png_res <- 160L

# 其他运行参数。
seed <- 42L
continue_on_error <- TRUE

set.seed(seed)

# 四个主输出目录，名称按需求固定。
category_dirs <- c(
  "human_disease_in",
  "human_disease_out",
  "mouse_disease_in",
  "mouse_disease_out"
)

# 四类组合的稳定索引前缀。
category_prefix <- c(
  human_disease_in = "HDI",
  human_disease_out = "HDO",
  mouse_disease_in = "MDI",
  mouse_disease_out = "MDO"
)

# ============================================================
# 3. 通用工具函数
# ============================================================

# 功能：将字符串转换为适合文件名的安全文本。
safe_filename <- function(x){
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "NA", x)
}

# 功能：判断 disease_group 是否为可用于组合的有效值。
is_valid_disease_group <- function(x){
  !is.na(x) & nzchar(trimws(as.character(x)))
}

# 功能：基础配置自检，避免路径、组合规模或维度参数错误。
validate_config <- function(){
  if(!dir.exists(input_dir)){
    stop(paste0("[ERROR] 输入目录不存在：", input_dir), call. = FALSE)
  }
  if(!file.exists(metadata_file)){
    stop(paste0("[ERROR] 样本信息表不存在：", metadata_file), call. = FALSE)
  }
  if(in_combination_dataset_size != 4L){
    stop("[ERROR] disease_in 当前要求 in_combination_dataset_size 固定为 4。", call. = FALSE)
  }
  if(out_combination_disease_size != 6L){
    stop("[ERROR] disease_out 当前要求 out_combination_disease_size 固定为 6。", call. = FALSE)
  }
  if(out_random_combination_n != 20L){
    stop("[ERROR] disease_out 当前要求随机组合数固定为 20。", call. = FALSE)
  }
  if(max_merged_cells <= 0){
    stop("[ERROR] max_merged_cells 必须为正整数。", call. = FALSE)
  }
  if(combination_enumeration_threshold < out_random_combination_n){
    stop("[ERROR] combination_enumeration_threshold 不能小于随机组合数。", call. = FALSE)
  }
  if(max(dims_use) > npcs){
    stop("[ERROR] dims_use 最大值不能超过 npcs。", call. = FALSE)
  }
  if(cluster_resolution != 0.1){
    warning("[WARN] 当前 resolution 不是需求指定的 0.1。", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for(d in category_dirs){
    dir.create(file.path(output_dir, d), recursive = TRUE, showWarnings = FALSE)
  }

  cat("[CHECK] 基础配置自检通过。\n")
  cat("[CHECK] disease_in：1 dataset 也产生组合；2-3 dataset 全并；>=4 dataset 时 C(n,4)。\n")
  cat("[CHECK] disease_out：<6 disease_group 全并；>=6 disease_group 时随机选择最多",
      out_random_combination_n, "个 6-disease_group 组合。\n")
  cat("[CHECK] 每个组合最大细胞数：", max_merged_cells, "\n")
}

# 功能：读取并检查样本信息表；组合字段统一使用 disease_group。
# 说明：disease 列即使存在也不参与本脚本的 disease_in / disease_out 分组。
read_and_validate_metadata <- function(file){
  cat("[READ] 样本信息表：", file, "\n")
  meta_df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  required_cols <- c("dataset_id", "species", "disease_group")
  missing_cols <- setdiff(required_cols, colnames(meta_df))
  if(length(missing_cols) > 0){
    stop(
      paste0(
        "[ERROR] 样本信息表缺少必要列：",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  meta_df$dataset_id <- trimws(as.character(meta_df$dataset_id))
  meta_df$species <- trimws(as.character(meta_df$species))
  meta_df$disease_group <- ifelse(
    is.na(meta_df$disease_group),
    NA_character_,
    trimws(as.character(meta_df$disease_group))
  )

  if(any(!nzchar(meta_df$dataset_id))){
    stop("[ERROR] 样本信息表存在空 dataset_id。", call. = FALSE)
  }

  if(anyDuplicated(meta_df$dataset_id)){
    dup_ids <- unique(meta_df$dataset_id[duplicated(meta_df$dataset_id)])
    stop(
      paste0(
        "[ERROR] 样本信息表 dataset_id 必须唯一；发现重复：",
        paste(head(dup_ids, 20), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unknown_species <- setdiff(unique(meta_df$species), c("Human", "Mouse"))
  unknown_species <- unknown_species[!is.na(unknown_species)]
  if(length(unknown_species) > 0){
    warning(
      paste0(
        "[WARN] 发现 Human/Mouse 之外的 species，将不进入四类组合：",
        paste(unknown_species, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  n_missing_group <- sum(!is_valid_disease_group(meta_df$disease_group))
  if(n_missing_group > 0){
    cat(
      "[WARN] disease_group 缺失/空值记录：", n_missing_group,
      "；这些 dataset 不参与 disease_in/out 组合。\n"
    )
  }

  cat("[DATA] 样本信息表：", nrow(meta_df), "行 x", ncol(meta_df), "列\n")
  cat(
    "[CHECK] 有效 disease_group 数：",
    length(unique(meta_df$disease_group[is_valid_disease_group(meta_df$disease_group)])),
    "\n"
  )
  meta_df
}

# 功能：优先根据 RDS 文件名匹配 dataset_id；若多个命中，选择最长 dataset_id 并报告歧义。
match_dataset_id_from_filename <- function(rds_file, dataset_ids){
  stem <- tools::file_path_sans_ext(basename(rds_file))
  stem_upper <- toupper(stem)
  ids_upper <- toupper(dataset_ids)

  exact_idx <- which(ids_upper == stem_upper)
  if(length(exact_idx) == 1){
    return(list(dataset_id = dataset_ids[exact_idx], method = "filename_exact", ambiguous = FALSE))
  }

  hit_idx <- which(vapply(ids_upper, function(id) grepl(id, stem_upper, fixed = TRUE), logical(1)))
  if(length(hit_idx) == 0){
    return(list(dataset_id = NA_character_, method = "unmatched", ambiguous = FALSE))
  }
  if(length(hit_idx) == 1){
    return(list(dataset_id = dataset_ids[hit_idx], method = "filename_contains", ambiguous = FALSE))
  }

  hit_ids <- dataset_ids[hit_idx]
  selected <- hit_ids[which.max(nchar(hit_ids))]
  list(dataset_id = selected, method = "filename_contains_longest", ambiguous = TRUE)
}

# 功能：文件名无法匹配时，从 Seurat metadata$dataset_id 或 project.name 尝试回退匹配。
match_dataset_id_from_object <- function(obj, dataset_ids){
  if("dataset_id" %in% colnames(obj@meta.data)){
    object_ids <- unique(trimws(as.character(obj@meta.data$dataset_id)))
    object_ids <- object_ids[!is.na(object_ids) & nzchar(object_ids)]
    hit <- intersect(object_ids, dataset_ids)
    if(length(hit) == 1){
      return(list(dataset_id = hit, method = "object_metadata"))
    }
  }
  project_name <- tryCatch(as.character(obj@project.name), error = function(e) NA_character_)
  if(length(project_name) == 1 && !is.na(project_name) && project_name %in% dataset_ids){
    return(list(dataset_id = project_name, method = "project_name"))
  }
  list(dataset_id = NA_character_, method = "unmatched")
}

# 功能：逐个读取 RDS，获取实际维度并建立 RDS 与样本表的匹配清单。
# disease_group 从样本信息表写入 inventory，后续所有组合均以该字段为准。
scan_rds_inventory <- function(rds_files, meta_df){
  inventory_list <- vector("list", length(rds_files))

  for(i in seq_along(rds_files)){
    rds_file <- rds_files[i]
    cat("[READ]", i, "/", length(rds_files), "|", basename(rds_file), "\n")

    file_match <- match_dataset_id_from_filename(rds_file, meta_df$dataset_id)
    obj <- tryCatch(readRDS(rds_file), error = function(e) e)

    if(inherits(obj, "error")){
      inventory_list[[i]] <- data.frame(
        rds_file = normalizePath(rds_file, winslash = "/", mustWork = FALSE),
        rds_name = basename(rds_file),
        dataset_id = file_match$dataset_id,
        match_method = file_match$method,
        match_ambiguous = file_match$ambiguous,
        object_class = "READ_ERROR",
        has_rna_assay = FALSE,
        n_cells = NA_integer_,
        n_features = NA_integer_,
        species = NA_character_,
        disease_group = NA_character_,
        dataset_name = NA_character_,
        status = paste0("READ_ERROR: ", conditionMessage(obj)),
        stringsAsFactors = FALSE
      )
      cat("[WARN] RDS 读取失败：", basename(rds_file), " | ", conditionMessage(obj), "\n")
      next
    }

    if(!inherits(obj, "Seurat")){
      inventory_list[[i]] <- data.frame(
        rds_file = normalizePath(rds_file, winslash = "/", mustWork = FALSE),
        rds_name = basename(rds_file),
        dataset_id = file_match$dataset_id,
        match_method = file_match$method,
        match_ambiguous = file_match$ambiguous,
        object_class = paste(class(obj), collapse = ";"),
        has_rna_assay = FALSE,
        n_cells = NA_integer_,
        n_features = NA_integer_,
        species = NA_character_,
        disease_group = NA_character_,
        dataset_name = NA_character_,
        status = "NOT_SEURAT",
        stringsAsFactors = FALSE
      )
      cat("[WARN] 非 Seurat 对象，跳过：", basename(rds_file), "\n")
      rm(obj)
      gc(verbose = FALSE)
      next
    }

    final_match <- file_match
    if(is.na(final_match$dataset_id)){
      final_match <- match_dataset_id_from_object(obj, meta_df$dataset_id)
      final_match$ambiguous <- FALSE
    }

    has_rna <- "RNA" %in% Assays(obj)
    n_cells <- ncol(obj)
    n_features <- if(has_rna) nrow(obj[["RNA"]]) else nrow(obj)

    if(is.na(final_match$dataset_id)){
      status <- "UNMATCHED_DATASET_ID"
      species <- NA_character_
      disease_group <- NA_character_
      dataset_name <- NA_character_
    } else {
      meta_row <- meta_df[
        match(final_match$dataset_id, meta_df$dataset_id),
        ,
        drop = FALSE
      ]
      species <- meta_row$species[1]
      disease_group <- meta_row$disease_group[1]
      dataset_name <- if("dataset_name" %in% colnames(meta_row)){
        meta_row$dataset_name[1]
      } else {
        NA_character_
      }
      status <- if(has_rna) "OK" else "NO_RNA_ASSAY"
    }

    inventory_list[[i]] <- data.frame(
      rds_file = normalizePath(rds_file, winslash = "/", mustWork = FALSE),
      rds_name = basename(rds_file),
      dataset_id = final_match$dataset_id,
      match_method = final_match$method,
      match_ambiguous = ifelse(
        is.null(final_match$ambiguous),
        FALSE,
        final_match$ambiguous
      ),
      object_class = paste(class(obj), collapse = ";"),
      has_rna_assay = has_rna,
      n_cells = as.integer(n_cells),
      n_features = as.integer(n_features),
      species = species,
      disease_group = disease_group,
      dataset_name = dataset_name,
      status = status,
      stringsAsFactors = FALSE
    )

    cat(
      "[DATA]", basename(rds_file),
      "| dataset_id =", ifelse(
        is.na(final_match$dataset_id),
        "UNMATCHED",
        final_match$dataset_id
      ),
      "|", n_cells, "cells x", n_features, "genes",
      "| species =", ifelse(is.na(species), "NA", species),
      "| disease_group =", ifelse(is.na(disease_group), "NA", disease_group),
      "\n"
    )

    rm(obj)
    gc(verbose = FALSE)
  }

  do.call(rbind, inventory_list)
}

# 功能：从扫描清单中提取可参与组合的数据，并检查同一 dataset_id 是否对应多个 RDS。
build_eligible_catalog <- function(inventory){
  catalog <- inventory[
    inventory$status == "OK" &
      inventory$species %in% c("Human", "Mouse") &
      is_valid_disease_group(inventory$disease_group) &
      !is.na(inventory$n_cells) &
      inventory$n_cells > 0,
    ,
    drop = FALSE
  ]

  if(nrow(catalog) == 0){
    return(catalog)
  }

  duplicated_ids <- unique(catalog$dataset_id[duplicated(catalog$dataset_id)])
  if(length(duplicated_ids) > 0){
    stop(
      paste0(
        "[ERROR] 同一 dataset_id 匹配到多个 RDS，无法安全决定哪个文件代表该 dataset：",
        paste(duplicated_ids, collapse = ", "),
        "。请先整理 input_dir，使每个 dataset_id 仅对应一个 RDS。"
      ),
      call. = FALSE
    )
  }

  catalog$catalog_row_id <- seq_len(nrow(catalog))
  rownames(catalog) <- NULL
  catalog
}

# 功能：为 disease_in 根据同一 disease_group 下的 dataset 数生成组合。
# 规则：
#   1) 1 个 dataset：也生成 1 个单 dataset 组合；
#   2) 2-3 个 dataset：全部 dataset 直接合并；
#   3) >=4 个 dataset：生成所有不重复的 4-dataset 组合。
generate_in_dataset_sets <- function(df){
  if(nrow(df) < 1){
    return(list())
  }

  df <- df[order(df$dataset_id), , drop = FALSE]
  n_dataset <- nrow(df)

  if(n_dataset < in_combination_dataset_size){
    return(list(df$catalog_row_id))
  }

  combn(
    df$catalog_row_id,
    m = in_combination_dataset_size,
    simplify = FALSE
  )
}

# 功能：从 n 个 disease_group 中随机获得最多 20 个不重复的 6-disease_group 组合。
# 当理论组合数较小时完整枚举后随机抽样；理论组合数很大时直接随机生成唯一组合。
sample_out_disease_group_sets <- function(disease_group_values, sample_seed){
  disease_group_values <- sort(unique(as.character(disease_group_values)))
  disease_group_values <- disease_group_values[
    !is.na(disease_group_values) &
      nzchar(trimws(disease_group_values))
  ]
  n_group <- length(disease_group_values)

  if(n_group < 2){
    return(list(
      sets = list(),
      total_possible = 0,
      sampling_method = "none"
    ))
  }

  # 少于 6 个 disease_group 时，直接合并该物种全部有效 disease_group。
  if(n_group < out_combination_disease_size){
    return(list(
      sets = list(disease_group_values),
      total_possible = 1,
      sampling_method = "merge_all_lt6_disease_groups"
    ))
  }

  total_possible <- choose(n_group, out_combination_disease_size)
  n_take <- min(as.numeric(out_random_combination_n), total_possible)
  n_take <- as.integer(n_take)

  set.seed(sample_seed)

  if(total_possible <= combination_enumeration_threshold){
    all_sets <- combn(
      disease_group_values,
      m = out_combination_disease_size,
      simplify = FALSE
    )

    if(length(all_sets) <= n_take){
      selected_sets <- all_sets
    } else {
      selected_idx <- sample(
        seq_along(all_sets),
        size = n_take,
        replace = FALSE
      )
      selected_sets <- all_sets[selected_idx]
    }

    return(list(
      sets = selected_sets,
      total_possible = total_possible,
      sampling_method = "enumerate_then_random_sample"
    ))
  }

  # 组合空间很大时，不完整生成 combn 结果，直接随机获得 20 个唯一组合。
  selected_sets <- list()
  selected_keys <- character()

  while(length(selected_sets) < n_take){
    candidate <- sort(sample(
      disease_group_values,
      size = out_combination_disease_size,
      replace = FALSE
    ))
    key <- paste(candidate, collapse = "\r")

    if(!key %in% selected_keys){
      selected_sets[[length(selected_sets) + 1L]] <- candidate
      selected_keys <- c(selected_keys, key)
    }
  }

  list(
    sets = selected_sets,
    total_possible = total_possible,
    sampling_method = "direct_random_unique_sample"
  )
}

# 功能：返回 disease_out 的可重复随机种子；Human / Mouse 使用不同 seed。
get_out_sampling_seed <- function(species){
  if(species == "Human"){
    return(seed + 10001L)
  }
  if(species == "Mouse"){
    return(seed + 20001L)
  }
  stop(paste0("[ERROR] 未知 species：", species), call. = FALSE)
}

# 功能：按指定类别生成组合。
# disease_in：同一 disease_group 内按 dataset 数量生成组合。
# disease_out：按 disease_group 类型组合；>=6 个 disease_group 时随机选择最多 20 个六组组合，
#              每个组合合并这 6 个 disease_group 下的全部 dataset。
build_category_combinations <- function(catalog, category){
  species_need <- if(grepl("^human_", category)) "Human" else "Mouse"
  mode_need <- if(grepl("_disease_in$", category)) "disease_in" else "disease_out"

  df <- catalog[catalog$species == species_need, , drop = FALSE]
  if(nrow(df) < 1){
    return(list())
  }

  combos <- list()

  if(mode_need == "disease_in"){
    group_names <- sort(unique(df$disease_group))

    for(group_name in group_names){
      group_df <- df[df$disease_group == group_name, , drop = FALSE]
      dataset_sets <- generate_in_dataset_sets(group_df)

      for(ids in dataset_sets){
        selected <- catalog[
          match(ids, catalog$catalog_row_id),
          ,
          drop = FALSE
        ]

        raw_cells <- sum(as.numeric(selected$n_cells))

        combination_rule <- if(nrow(group_df) == 1){
          "single_dataset"
        } else if(nrow(group_df) < in_combination_dataset_size){
          "merge_all_2to3_datasets"
        } else {
          "choose_4_datasets_of_n"
        }

        combos[[length(combos) + 1L]] <- list(
          catalog_row_ids = ids,
          species = species_need,
          combination_type = "disease_in",
          combination_rule = combination_rule,
          dataset_count = length(ids),
          disease_group_count = 1L,
          disease_group_info = group_name,
          merge_info = paste0(
            "same_disease_group=", group_name,
            ";datasets=", paste(selected$dataset_id, collapse = " | ")
          ),
          raw_cells = raw_cells,
          planned_final_cells = as.integer(
            min(raw_cells, as.numeric(max_merged_cells))
          ),
          max_cells = max_merged_cells,
          out_total_possible_combinations = NA_real_,
          out_sampling_method = NA_character_
        )
      }
    }

    return(combos)
  }

  # disease_out：至少需要两个不同 disease_group。
  disease_group_values <- sort(unique(df$disease_group))
  if(length(disease_group_values) < 2){
    return(list())
  }

  sampled <- sample_out_disease_group_sets(
    disease_group_values = disease_group_values,
    sample_seed = get_out_sampling_seed(species_need)
  )

  if(length(sampled$sets) == 0){
    return(list())
  }

  if(length(disease_group_values) >= out_combination_disease_size){
    cat(
      "[CHECK]", category,
      "| disease_group types =", length(disease_group_values),
      "| possible C(n,6) =", format(sampled$total_possible, scientific = FALSE),
      "| selected =", length(sampled$sets),
      "| method =", sampled$sampling_method,
      "\n"
    )

    if(sampled$total_possible < out_random_combination_n){
      cat(
        "[WARN]", category,
        "理论 6-disease_group 组合数少于", out_random_combination_n,
        "，将运行全部", length(sampled$sets), "个可用组合。\n"
      )
    }
  } else {
    cat(
      "[CHECK]", category,
      "| disease_group types =", length(disease_group_values),
      "<", out_combination_disease_size,
      "，直接合并该物种全部可用 dataset。\n"
    )
  }

  for(group_set in sampled$sets){
    selected <- df[
      df$disease_group %in% group_set,
      ,
      drop = FALSE
    ]
    selected <- selected[
      order(selected$disease_group, selected$dataset_id),
      ,
      drop = FALSE
    ]

    if(nrow(selected) < 2){
      next
    }

    selected_groups <- sort(unique(selected$disease_group))
    raw_cells <- sum(as.numeric(selected$n_cells))

    combination_rule <- if(
      length(disease_group_values) < out_combination_disease_size
    ){
      "merge_all_lt6_disease_groups"
    } else {
      "random_6_disease_group_combination"
    }

    combos[[length(combos) + 1L]] <- list(
      catalog_row_ids = selected$catalog_row_id,
      species = species_need,
      combination_type = "disease_out",
      combination_rule = combination_rule,
      dataset_count = nrow(selected),
      disease_group_count = length(selected_groups),
      disease_group_info = paste(selected_groups, collapse = " | "),
      merge_info = paste0(
        "disease_groups=", paste(selected_groups, collapse = " | "),
        ";datasets=", paste(selected$dataset_id, collapse = " | ")
      ),
      raw_cells = raw_cells,
      planned_final_cells = as.integer(
        min(raw_cells, as.numeric(max_merged_cells))
      ),
      max_cells = max_merged_cells,
      out_total_possible_combinations = as.numeric(sampled$total_possible),
      out_sampling_method = sampled$sampling_method
    )
  }

  combos
}

# 功能：为某一类别生成组合计划。
plan_category <- function(catalog, category){
  build_category_combinations(catalog, category)
}

# 功能：检查生成的组合计划是否符合 disease_group、单 dataset disease_in 和 20 万细胞上限规划。
validate_combination_plan <- function(combos, catalog, category){
  if(length(combos) == 0){
    cat("[CHECK]", category, "无可执行组合。\n")
    return(invisible(TRUE))
  }

  mode_need <- if(grepl("_disease_in$", category)) "disease_in" else "disease_out"

  if(mode_need == "disease_out" && length(combos) > out_random_combination_n){
    stop(
      paste0(
        "[ERROR] ", category,
        " disease_out 组合数超过随机上限 ",
        out_random_combination_n, "：", length(combos)
      ),
      call. = FALSE
    )
  }

  combo_keys <- character()

  for(combo in combos){
    selected <- catalog[
      match(combo$catalog_row_ids, catalog$catalog_row_id),
      ,
      drop = FALSE
    ]

    if(any(is.na(selected$catalog_row_id))){
      stop(
        paste0("[ERROR] ", category, " 组合索引包含不存在的 catalog_row_id。"),
        call. = FALSE
      )
    }

    selected_groups <- sort(unique(selected$disease_group))
    selected_datasets <- sort(unique(selected$dataset_id))

    if(combo$combination_type == "disease_in"){
      if(nrow(selected) < 1){
        stop(
          paste0("[ERROR] ", category, " disease_in 组合没有 dataset。"),
          call. = FALSE
        )
      }
      if(length(selected_groups) != 1){
        stop(
          paste0("[ERROR] ", category, " disease_in 计划中 disease_group 不唯一。"),
          call. = FALSE
        )
      }

      if(combo$combination_rule == "single_dataset" &&
         length(selected_datasets) != 1){
        stop(
          paste0("[ERROR] ", category, " single_dataset 组合不是 1 个 dataset。"),
          call. = FALSE
        )
      }

      if(combo$combination_rule == "merge_all_2to3_datasets" &&
         !(length(selected_datasets) %in% 2:3)){
        stop(
          paste0("[ERROR] ", category, " disease_in 直接合并组合不是 2-3 个 dataset。"),
          call. = FALSE
        )
      }

      if(combo$combination_rule == "choose_4_datasets_of_n" &&
         length(selected_datasets) != in_combination_dataset_size){
        stop(
          paste0("[ERROR] ", category, " disease_in 四数据集组合不是 4 个 dataset。"),
          call. = FALSE
        )
      }
    }

    if(combo$combination_type == "disease_out"){
      if(nrow(selected) < 2){
        stop(
          paste0("[ERROR] ", category, " disease_out 组合 dataset 数少于 2。"),
          call. = FALSE
        )
      }
      if(length(selected_groups) < 2){
        stop(
          paste0("[ERROR] ", category, " disease_out 计划少于 2 个 disease_group。"),
          call. = FALSE
        )
      }

      if(combo$combination_rule == "random_6_disease_group_combination" &&
         length(selected_groups) != out_combination_disease_size){
        stop(
          paste0(
            "[ERROR] ", category,
            " 随机 disease_out 组合不是 ",
            out_combination_disease_size,
            " 个 disease_group。"
          ),
          call. = FALSE
        )
      }

      if(combo$combination_rule == "merge_all_lt6_disease_groups" &&
         length(selected_groups) >= out_combination_disease_size){
        stop(
          paste0(
            "[ERROR] ", category,
            " merge_all_lt6_disease_groups 实际 disease_group 数 >= ",
            out_combination_disease_size,
            "。"
          ),
          call. = FALSE
        )
      }

      group_key <- paste(selected_groups, collapse = "\r")
      if(group_key %in% combo_keys){
        stop(
          paste0("[ERROR] ", category, " disease_out 出现重复 disease_group 组合。"),
          call. = FALSE
        )
      }
      combo_keys <- c(combo_keys, group_key)
    }

    raw_cells_check <- sum(as.numeric(selected$n_cells))
    if(!isTRUE(all.equal(as.numeric(combo$raw_cells), raw_cells_check))){
      stop(
        paste0("[ERROR] ", category, " raw_cells 与 catalog 不一致。"),
        call. = FALSE
      )
    }

    expected_planned <- min(
      raw_cells_check,
      as.numeric(max_merged_cells)
    )
    if(!isTRUE(all.equal(
      as.numeric(combo$planned_final_cells),
      expected_planned
    ))){
      stop(
        paste0("[ERROR] ", category, " planned_final_cells 计算错误。"),
        call. = FALSE
      )
    }
  }

  cat("[CHECK]", category, "组合计划校验通过：", length(combos), "个组合。\n")
  invisible(TRUE)
}

# 功能：为组合赋予稳定索引名称，后续所有保存文件均使用该 combo_id。
assign_combo_ids <- function(combos, category){
  if(length(combos) == 0){
    return(combos)
  }
  prefix <- category_prefix[[category]]
  for(i in seq_along(combos)){
    combos[[i]]$combo_id <- sprintf("%s_%04d", prefix, i)
    combos[[i]]$category <- category
  }
  combos
}

# 功能：按 dataset 原始细胞数比例分配最终细胞配额。
# 若总细胞数 <= target_cells，则返回原始细胞数；否则精确分配到 target_cells。
# 注意：所有中间算术统一使用 double/numeric，避免大细胞数相乘触发 R 32 位整数溢出。
allocate_dataset_cells <- function(counts, target_cells){
  counts_num <- as.numeric(counts)
  target_num <- as.numeric(target_cells)

  if(length(counts_num) == 0){
    stop("[ERROR] allocate_dataset_cells: counts 不能为空。", call. = FALSE)
  }
  if(anyNA(counts_num) || any(!is.finite(counts_num)) || any(counts_num <= 0)){
    stop(
      "[ERROR] allocate_dataset_cells: counts 必须全部为有限正数。",
      call. = FALSE
    )
  }
  if(length(target_num) != 1 ||
     is.na(target_num) ||
     !is.finite(target_num) ||
     target_num <= 0){
    stop(
      "[ERROR] allocate_dataset_cells: target_cells 必须是单个有限正数。",
      call. = FALSE
    )
  }

  # 细胞数必须是整数；显式检查，避免异常输入被静默截断。
  tol <- .Machine$double.eps^0.5
  if(any(abs(counts_num - round(counts_num)) > tol)){
    stop(
      "[ERROR] allocate_dataset_cells: counts 中存在非整数细胞数。",
      call. = FALSE
    )
  }
  if(abs(target_num - round(target_num)) > tol){
    stop(
      "[ERROR] allocate_dataset_cells: target_cells 必须为整数。",
      call. = FALSE
    )
  }

  counts_num <- round(counts_num)
  target_num <- round(target_num)
  total_cells <- sum(counts_num)

  if(!is.finite(total_cells)){
    stop(
      "[ERROR] allocate_dataset_cells: 总细胞数计算结果不是有限值。",
      call. = FALSE
    )
  }

  # 总细胞数不超过上限时，保持全部细胞，不下采样。
  if(total_cells <= target_num){
    return(as.integer(counts_num))
  }

  if(target_num < length(counts_num)){
    stop(
      "[ERROR] target_cells 小于 dataset 数量，无法保证每个 dataset 至少保留 1 个细胞。",
      call. = FALSE
    )
  }

  # 关键修复：
  # target_num/counts_num 均为 double，因此不会发生 32 位整数乘法溢出。
  desired <- target_num * counts_num / total_cells

  if(anyNA(desired) || any(!is.finite(desired))){
    stop(
      "[ERROR] allocate_dataset_cells: 比例配额计算产生 NA/Inf。",
      call. = FALSE
    )
  }

  allocation <- floor(desired)
  allocation[allocation < 1] <- 1
  allocation <- pmin(allocation, counts_num)

  # 将 floor 后剩余的细胞逐个补到“理论配额差值”最大的 dataset。
  while(sum(allocation) < target_num){
    candidates <- which(allocation < counts_num)
    if(length(candidates) == 0){
      break
    }

    score <- desired[candidates] - allocation[candidates]
    if(anyNA(score) || any(!is.finite(score))){
      stop(
        "[ERROR] allocate_dataset_cells: 增加配额时 score 出现 NA/Inf。",
        call. = FALSE
      )
    }

    chosen <- candidates[which.max(score)]
    allocation[chosen] <- allocation[chosen] + 1
  }

  # 若因每个 dataset 至少保留 1 个细胞造成超额，则逐个减回目标值。
  while(sum(allocation) > target_num){
    candidates <- which(allocation > 1)
    if(length(candidates) == 0){
      break
    }

    score <- allocation[candidates] - desired[candidates]
    if(anyNA(score) || any(!is.finite(score))){
      stop(
        "[ERROR] allocate_dataset_cells: 减少配额时 score 出现 NA/Inf。",
        call. = FALSE
      )
    }

    chosen <- candidates[which.max(score)]
    allocation[chosen] <- allocation[chosen] - 1
  }

  if(anyNA(allocation) || any(!is.finite(allocation))){
    stop(
      "[ERROR] allocate_dataset_cells: 最终 allocation 出现 NA/Inf。",
      call. = FALSE
    )
  }

  if(sum(allocation) != target_num){
    stop(
      paste0(
        "[ERROR] dataset 配额分配失败：",
        format(sum(allocation), scientific = FALSE),
        " != ",
        format(target_num, scientific = FALSE)
      ),
      call. = FALSE
    )
  }

  if(any(allocation > counts_num)){
    stop(
      "[ERROR] dataset 配额超过对应 dataset 原始细胞数。",
      call. = FALSE
    )
  }

  if(any(allocation < 1)){
    stop(
      "[ERROR] dataset 配额出现小于 1 的值。",
      call. = FALSE
    )
  }

  # 目标上限为 200000，安全低于 R integer 上限，最终转回 integer。
  as.integer(allocation)
}

# 功能：创建空索引模板，确保无有效组合时仍输出字段完整的 CSV。
empty_index_template <- function(){
  data.frame(
    combo_id = character(),
    category = character(),
    combination_type = character(),
    combination_rule = character(),
    species = character(),
    disease_group_info = character(),
    merge_info = character(),
    combination_dataset_count = integer(),
    combination_disease_group_count = integer(),
    out_total_possible_combinations = numeric(),
    out_sampling_method = character(),
    dataset_id = character(),
    disease_group = character(),
    dataset_n_cells = numeric(),
    dataset_planned_cells = integer(),
    dataset_downsample_required = logical(),
    combo_raw_cells = numeric(),
    planned_final_cells = integer(),
    max_cells = integer(),
    rds_file = character(),
    status = character(),
    final_cells = integer(),
    output_rds = character(),
    output_metadata_csv = character(),
    output_combination_txt = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

# 功能：将组合计划展开成长表索引；每行代表一个 combo_id 中的一个 dataset_id。
build_index_table <- function(combos, catalog){
  if(length(combos) == 0){
    return(empty_index_template())
  }

  index_rows <- list()

  for(combo in combos){
    selected <- catalog[
      match(combo$catalog_row_ids, catalog$catalog_row_id),
      ,
      drop = FALSE
    ]

    planned_cells <- allocate_dataset_cells(
      selected$n_cells,
      target_cells = combo$max_cells
    )
    planned_final <- sum(as.numeric(planned_cells))

    for(i in seq_len(nrow(selected))){
      index_rows[[length(index_rows) + 1L]] <- data.frame(
        combo_id = combo$combo_id,
        category = combo$category,
        combination_type = combo$combination_type,
        combination_rule = combo$combination_rule,
        species = combo$species,
        disease_group_info = combo$disease_group_info,
        merge_info = combo$merge_info,
        combination_dataset_count = nrow(selected),
        combination_disease_group_count = length(unique(selected$disease_group)),
        out_total_possible_combinations = ifelse(
          is.null(combo$out_total_possible_combinations),
          NA_real_,
          as.numeric(combo$out_total_possible_combinations)
        ),
        out_sampling_method = ifelse(
          is.null(combo$out_sampling_method),
          NA_character_,
          combo$out_sampling_method
        ),
        dataset_id = selected$dataset_id[i],
        disease_group = selected$disease_group[i],
        dataset_n_cells = as.numeric(selected$n_cells[i]),
        dataset_planned_cells = as.integer(planned_cells[i]),
        dataset_downsample_required = planned_cells[i] < selected$n_cells[i],
        combo_raw_cells = as.numeric(combo$raw_cells),
        planned_final_cells = as.integer(planned_final),
        max_cells = as.integer(combo$max_cells),
        rds_file = selected$rds_file[i],
        status = "PLANNED",
        final_cells = NA_integer_,
        output_rds = NA_character_,
        output_metadata_csv = NA_character_,
        output_combination_txt = NA_character_,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, index_rows)
}

# 功能：运行前内部逻辑自检，校验 disease_group、单 dataset disease_in、
#       随机 6-disease_group disease_out 和 20 万细胞配额。
run_internal_logic_check <- function(){
  # disease_in：1 个 dataset 必须产生 1 个组合。
  toy_in1 <- data.frame(
    dataset_id = "IN1_D1",
    disease_group = "Group_A",
    species = "Human",
    n_cells = 50000L,
    catalog_row_id = 1L,
    stringsAsFactors = FALSE
  )

  sets1 <- generate_in_dataset_sets(toy_in1)
  if(length(sets1) != 1 || length(sets1[[1]]) != 1){
    stop(
      "[ERROR] 内部自检失败：单 dataset disease_in 应生成 1 个组合。",
      call. = FALSE
    )
  }

  # disease_in：3 个 dataset 应直接形成一个组合。
  toy_in3 <- data.frame(
    dataset_id = paste0("IN3_D", 1:3),
    disease_group = rep("Group_B", 3),
    species = rep("Human", 3),
    n_cells = c(30000L, 40000L, 50000L),
    catalog_row_id = 11:13,
    stringsAsFactors = FALSE
  )

  sets3 <- generate_in_dataset_sets(toy_in3)
  if(length(sets3) != 1 || length(sets3[[1]]) != 3){
    stop(
      "[ERROR] 内部自检失败：3-dataset disease_in 应直接生成 1 个组合。",
      call. = FALSE
    )
  }

  # disease_in：5 个 dataset 应产生 C(5,4)=5 个组合。
  toy_in5 <- data.frame(
    dataset_id = paste0("IN5_D", 1:5),
    disease_group = rep("Group_C", 5),
    species = rep("Human", 5),
    n_cells = rep(30000L, 5),
    catalog_row_id = 21:25,
    stringsAsFactors = FALSE
  )

  sets5 <- generate_in_dataset_sets(toy_in5)
  if(length(sets5) != choose(5, 4) ||
     !all(vapply(sets5, length, integer(1)) == 4L)){
    stop(
      "[ERROR] 内部自检失败：5-dataset disease_in 的 C(5,4) 不正确。",
      call. = FALSE
    )
  }

  # disease_out：5 个 disease_group 时应直接形成一个组合。
  out5 <- sample_out_disease_group_sets(
    paste0("Group_", 1:5),
    sample_seed = seed + 1L
  )
  if(length(out5$sets) != 1 ||
     length(out5$sets[[1]]) != 5 ||
     out5$sampling_method != "merge_all_lt6_disease_groups"){
    stop(
      "[ERROR] 内部自检失败：<6 disease_group 的直接合并规则不正确。",
      call. = FALSE
    )
  }

  # disease_out：8 个 disease_group 时 C(8,6)=28，固定随机获得 20 个不重复组合。
  out8 <- sample_out_disease_group_sets(
    paste0("Group_", 1:8),
    sample_seed = seed + 2L
  )
  out8_keys <- vapply(
    out8$sets,
    function(x) paste(sort(x), collapse = "|"),
    character(1)
  )

  if(out8$total_possible != choose(8, 6) ||
     length(out8$sets) != out_random_combination_n ||
     any(vapply(out8$sets, length, integer(1)) != out_combination_disease_size) ||
     anyDuplicated(out8_keys)){
    stop(
      "[ERROR] 内部自检失败：随机 20 个 6-disease_group 组合不正确。",
      call. = FALSE
    )
  }

  out8_repeat <- sample_out_disease_group_sets(
    paste0("Group_", 1:8),
    sample_seed = seed + 2L
  )
  out8_repeat_keys <- vapply(
    out8_repeat$sets,
    function(x) paste(sort(x), collapse = "|"),
    character(1)
  )
  if(!identical(out8_keys, out8_repeat_keys)){
    stop(
      "[ERROR] 内部自检失败：disease_out 随机组合不能复现。",
      call. = FALSE
    )
  }

  # 验证 20 万细胞比例分配，并覆盖整数溢出场景。
  allocation <- allocate_dataset_cells(
    c(120000L, 90000L, 60000L, 40000L),
    target_cells = max_merged_cells
  )
  if(anyNA(allocation) ||
     sum(as.numeric(allocation)) != as.numeric(max_merged_cells) ||
     any(allocation > c(120000L, 90000L, 60000L, 40000L))){
    stop(
      "[ERROR] 内部自检失败：20 万细胞配额分配不正确。",
      call. = FALSE
    )
  }

  no_downsample <- allocate_dataset_cells(
    c(50000L, 40000L, 30000L),
    target_cells = max_merged_cells
  )
  if(!identical(no_downsample, c(50000L, 40000L, 30000L))){
    stop(
      "[ERROR] 内部自检失败：小于 20 万细胞时不应下采样。",
      call. = FALSE
    )
  }

  cat("[CHECK] disease_group、单 dataset、随机组合与细胞配额逻辑自检通过。\n")
}

# 功能：从原始对象中构造统一的 sample_id，优先使用 Sample/sample/sample_id/orig.ident。
derive_original_sample <- function(obj){
  meta <- obj@meta.data
  candidates <- c("Sample", "sample", "sample_id", "orig.ident")
  hit <- candidates[candidates %in% colnames(meta)]
  if(length(hit) == 0){
    return(rep("sample1", ncol(obj)))
  }
  value <- as.character(meta[[hit[1]]])
  value[is.na(value) | !nzchar(trimws(value))] <- "sample1"
  value
}

# 功能：读取单个组合成员，必要时在合并前按 dataset 配额下采样，
#       并写入统一的 dataset/species/disease_group/sample/组合 metadata。
prepare_seurat_for_merge <- function(sample_row, combo, planned_cells){
  obj <- readRDS(sample_row$rds_file)

  if(!inherits(obj, "Seurat")){
    stop(paste0("[ERROR] 非 Seurat 对象：", sample_row$rds_file), call. = FALSE)
  }
  if(!"RNA" %in% Assays(obj)){
    stop(paste0("[ERROR] 对象缺少 RNA assay：", sample_row$rds_file), call. = FALSE)
  }

  DefaultAssay(obj) <- "RNA"
  actual_cells <- ncol(obj)
  planned_cells <- as.integer(planned_cells)

  if(planned_cells <= 0 || planned_cells > actual_cells){
    stop(
      paste0(
        "[ERROR] ", sample_row$dataset_id,
        " 的 planned_cells 无效：", planned_cells,
        "；当前 RDS 实际细胞数：", actual_cells
      ),
      call. = FALSE
    )
  }

  if(planned_cells < actual_cells){
    keep_cells <- sample(
      colnames(obj),
      size = planned_cells,
      replace = FALSE
    )
    obj <- subset(obj, cells = keep_cells)
    cat(
      "[DATA]", combo$combo_id, "|", sample_row$dataset_id,
      "| pre-merge downsample =", actual_cells, "->", ncol(obj), "cells\n"
    )
  } else {
    cat(
      "[DATA]", combo$combo_id, "|", sample_row$dataset_id,
      "| no downsample =", actual_cells, "cells\n"
    )
  }

  original_sample <- derive_original_sample(obj)

  obj$dataset_id <- sample_row$dataset_id
  obj$species <- sample_row$species
  obj$disease_group <- sample_row$disease_group
  obj$sample_id <- paste(
    sample_row$dataset_id,
    original_sample,
    sep = "::"
  )

  obj$combination_id <- combo$combo_id
  obj$combination_type <- combo$combination_type
  obj$combination_rule <- combo$combination_rule
  obj$combination_disease_group_info <- combo$disease_group_info
  obj$combination_merge_info <- combo$merge_info
  obj$combination_species <- combo$species
  obj$combination_dataset_count <- combo$dataset_count
  obj$combination_disease_group_count <- combo$disease_group_count
  obj$combination_raw_cells <- combo$raw_cells
  obj$combination_planned_cells <- combo$planned_final_cells

  obj$out_total_possible_combinations <- ifelse(
    is.null(combo$out_total_possible_combinations),
    NA_real_,
    as.numeric(combo$out_total_possible_combinations)
  )
  obj$out_sampling_method <- ifelse(
    is.null(combo$out_sampling_method),
    NA_character_,
    combo$out_sampling_method
  )

  obj$dataset_raw_cells <- sample_row$n_cells
  obj$dataset_planned_cells <- planned_cells

  if(combo$combination_type == "disease_in"){
    obj$merge_label <- sample_row$dataset_id
  } else {
    obj$merge_label <- paste(
      sample_row$disease_group,
      sample_row$dataset_id,
      sep = " | "
    )
  }

  obj <- RenameCells(
    obj,
    add.cell.id = safe_filename(sample_row$dataset_id)
  )
  obj
}

# 功能：将 Seurat UMAP 的 3 个分量以静态 3D 散点图保存为 PNG。
# 参数：
#   data：已完成 3D UMAP 的 Seurat 对象；
#   group_col：用于着色/分组的 metadata 字段；
#   output_file：PNG 输出路径；
#   title_text：图标题。
# 说明：
#   - 使用 UMAP_1 / UMAP_2 / UMAP_3 三个轴；
#   - 不对细胞进行绘图下采样，最多展示当前组合中的全部 200,000 个细胞；
#   - 点大小会根据细胞数自动缩小，以降低大样本遮挡；
#   - ElbowPlot 属于 PCA/Harmony 维数诊断图，不是 UMAP 空间图，因此仍保持二维诊断形式。
save_umap_3d_plot <- function(data, group_col, output_file, title_text){
  if(!group_col %in% colnames(data@meta.data)){
    cat("[WARN] metadata 不存在", group_col, "，跳过 3D UMAP：",
        basename(output_file), "\n")
    return(invisible(FALSE))
  }

  if(!"umap" %in% names(data@reductions)){
    cat("[WARN] Seurat 对象不存在 umap reduction，跳过：",
        basename(output_file), "\n")
    return(invisible(FALSE))
  }

  umap_emb <- Embeddings(data[["umap"]])
  if(ncol(umap_emb) < 3){
    cat("[WARN] UMAP 少于 3 维，跳过 3D UMAP：",
        basename(output_file), "\n")
    return(invisible(FALSE))
  }

  cell_ids <- rownames(umap_emb)
  group_value <- as.character(data@meta.data[cell_ids, group_col])
  group_value[is.na(group_value) | !nzchar(trimws(group_value))] <- "NA"

  group_factor <- factor(group_value, levels = sort(unique(group_value)))
  group_levels <- levels(group_factor)
  n_groups <- length(group_levels)

  # 使用 base R 的 HCL 调色板；绘图点使用半透明颜色，图例使用不透明颜色。
  legend_colors <- grDevices::hcl.colors(
    max(n_groups, 1L),
    palette = "Dark 3"
  )
  point_colors <- grDevices::adjustcolor(
    legend_colors[as.integer(group_factor)],
    alpha.f = 0.55
  )

  n_cells <- nrow(umap_emb)
  point_cex <- if(n_cells >= 150000){
    0.16
  } else if(n_cells >= 100000){
    0.20
  } else if(n_cells >= 50000){
    0.28
  } else if(n_cells >= 20000){
    0.38
  } else {
    0.55
  }

  # 根据类别数动态调整图例布局。
  legend_ncol <- if(n_groups <= 20){
    1L
  } else if(n_groups <= 50){
    2L
  } else {
    3L
  }
  legend_cex <- if(n_groups <= 20){
    0.85
  } else if(n_groups <= 50){
    0.68
  } else {
    0.55
  }

  grDevices::png(
    filename = output_file,
    width = umap3d_png_width,
    height = umap3d_png_height,
    res = umap3d_png_res
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  # 左侧为 3D UMAP，右侧独立图例区域，避免图例覆盖数据点。
  graphics::layout(
    matrix(c(1, 2), nrow = 1),
    widths = c(4.5, 1.8)
  )

  graphics::par(mar = c(5, 5, 4, 2))
  scatterplot3d::scatterplot3d(
    x = umap_emb[, 1],
    y = umap_emb[, 2],
    z = umap_emb[, 3],
    color = point_colors,
    pch = 16,
    cex.symbols = point_cex,
    angle = umap3d_angle,
    grid = TRUE,
    box = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    zlab = "UMAP_3",
    main = title_text
  )

  graphics::par(mar = c(1, 1, 1, 1))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = group_levels,
    col = legend_colors,
    pch = 16,
    pt.cex = 0.9,
    cex = legend_cex,
    ncol = legend_ncol,
    bty = "n",
    title = group_col
  )

  cat("[SAVE]", output_file, "| 3D UMAP | cells =", n_cells,
      "| groups =", n_groups, "\n")
  invisible(TRUE)
}

# 功能：保存 Harmony ElbowPlot；若 Harmony reduction 无 stdev，则回退保存 PCA ElbowPlot 并记录警告。
save_elbow_plot <- function(data, output_file, combo_id){
  ok <- TRUE
  p <- tryCatch(
    ElbowPlot(data, ndims = npcs, reduction = "harmony"),
    error = function(e){
      ok <<- FALSE
      cat("[WARN]", combo_id, "Harmony ElbowPlot 失败，改用 PCA：", conditionMessage(e), "\n")
      ElbowPlot(data, ndims = npcs, reduction = "pca")
    }
  )
  ggsave(
    filename = output_file,
    plot = p,
    width = 8,
    height = 6,
    dpi = 150,
    units = "in"
  )
  cat("[SAVE]", output_file, if(ok) "" else " [PCA fallback]", "\n")
}

# 功能：写出每个组合的文本记录，包含 disease_group、组合规则和每个 dataset 的下采样配额。
write_combination_txt <- function(combo, selected, output_file){
  con <- file(output_file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  writeLines(c(
    paste0("combo_id\t", combo$combo_id),
    paste0("category\t", combo$category),
    paste0("combination_type\t", combo$combination_type),
    paste0("combination_rule\t", combo$combination_rule),
    paste0("species\t", combo$species),
    paste0("dataset_count\t", combo$dataset_count),
    paste0("disease_group_count\t", combo$disease_group_count),
    paste0(
      "out_total_possible_combinations\t",
      ifelse(
        is.null(combo$out_total_possible_combinations),
        "NA",
        format(combo$out_total_possible_combinations, scientific = FALSE)
      )
    ),
    paste0(
      "out_sampling_method\t",
      ifelse(
        is.null(combo$out_sampling_method),
        "NA",
        combo$out_sampling_method
      )
    ),
    paste0("disease_group_info\t", combo$disease_group_info),
    paste0("merge_info\t", combo$merge_info),
    paste0("raw_cells\t", combo$raw_cells),
    paste0("planned_final_cells\t", combo$planned_final_cells),
    paste0("max_cells\t", combo$max_cells),
    "",
    paste(
      "dataset_id",
      "disease_group",
      "species",
      "raw_cells",
      "planned_cells",
      "downsample_required",
      "rds_file",
      sep = "\t"
    )
  ), con = con)

  for(i in seq_len(nrow(selected))){
    writeLines(
      paste(
        selected$dataset_id[i],
        selected$disease_group[i],
        selected$species[i],
        selected$n_cells[i],
        selected$planned_cells[i],
        selected$planned_cells[i] < selected$n_cells[i],
        selected$rds_file[i],
        sep = "\t"
      ),
      con = con
    )
  }

  cat("[SAVE]", output_file, "\n")
}

# 功能：执行单个组合的完整 Seurat V4 + Harmony 分析，并保存所有要求结果。
# 关键规则：组合原始总细胞数 > 200,000 时，先按 dataset 比例下采样，再逐个合并。
process_one_combination <- function(combo, catalog, category_dir){
  selected <- catalog[
    match(combo$catalog_row_ids, catalog$catalog_row_id),
    ,
    drop = FALSE
  ]

  if(any(is.na(selected$catalog_row_id))){
    stop(paste0("[ERROR] ", combo$combo_id, " 存在无效 catalog_row_id。"), call. = FALSE)
  }

  selected$planned_cells <- allocate_dataset_cells(
    selected$n_cells,
    target_cells = combo$max_cells
  )
  expected_final_cells <- sum(selected$planned_cells)

  if(expected_final_cells != combo$planned_final_cells){
    stop(
      paste0(
        "[ERROR] ", combo$combo_id,
        " 规划细胞数不一致：", expected_final_cells,
        " vs ", combo$planned_final_cells
      ),
      call. = FALSE
    )
  }

  combo_dir <- file.path(category_dir, combo$combo_id)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)

  combination_txt <- file.path(combo_dir, paste0(combo$combo_id, "_combination.txt"))
  write_combination_txt(combo, selected, combination_txt)

  cat(
    "[STEP]", combo$combo_id,
    "|", combo$combination_type,
    "| rule =", combo$combination_rule,
    "| species =", combo$species,
    "| datasets =", nrow(selected),
    "| raw cells =", sum(selected$n_cells),
    "| planned cells =", expected_final_cells, "\n"
  )

  # 逐个读取、必要时下采样并顺序 merge，避免同时保留全部原始大对象。
  data <- NULL

  for(i in seq_len(nrow(selected))){
    cat(
      "[READ]", combo$combo_id,
      "|", selected$dataset_id[i],
      "|", basename(selected$rds_file[i]), "\n"
    )

    obj <- prepare_seurat_for_merge(
      sample_row = selected[i, , drop = FALSE],
      combo = combo,
      planned_cells = selected$planned_cells[i]
    )

    if(is.null(data)){
      data <- obj
    } else {
      data <- merge(data, y = obj)
    }

    rm(obj)
    gc(verbose = FALSE)
  }

  if(is.null(data)){
    stop(paste0("[ERROR] ", combo$combo_id, " 合并结果为空。"), call. = FALSE)
  }
  if(!"RNA" %in% Assays(data)){
    stop(paste0("[ERROR] ", combo$combo_id, " 合并后缺少 RNA assay。"), call. = FALSE)
  }

  DefaultAssay(data) <- "RNA"

  if(anyDuplicated(colnames(data))){
    stop(paste0("[ERROR] ", combo$combo_id, " 合并后出现重复 cell ID。"), call. = FALSE)
  }
  if(combo$combination_type == "disease_in"){
    if(length(unique(data$dataset_id)) < 1){
      stop(paste0("[ERROR] ", combo$combo_id, " disease_in 没有 dataset。"), call. = FALSE)
    }
    if(length(unique(data$disease_group)) != 1){
      stop(
        paste0("[ERROR] ", combo$combo_id, " disease_in 合并后 disease_group 不唯一。"),
        call. = FALSE
      )
    }
  }

  if(combo$combination_type == "disease_out"){
    if(length(unique(data$dataset_id)) < 2){
      stop(
        paste0("[ERROR] ", combo$combo_id, " disease_out 合并后 dataset_id 少于 2。"),
        call. = FALSE
      )
    }
    if(length(unique(data$disease_group)) < 2){
      stop(
        paste0("[ERROR] ", combo$combo_id, " disease_out 合并后 disease_group 少于 2 种。"),
        call. = FALSE
      )
    }
  }
  if(ncol(data) != expected_final_cells){
    stop(
      paste0(
        "[ERROR] ", combo$combo_id,
        " 合并后细胞数与规划不一致：", ncol(data),
        " vs ", expected_final_cells
      ),
      call. = FALSE
    )
  }
  if(ncol(data) > max_merged_cells){
    stop(
      paste0(
        "[ERROR] ", combo$combo_id,
        " 合并后细胞数超过上限 ", max_merged_cells,
        "：", ncol(data)
      ),
      call. = FALSE
    )
  }

  cat("[DATA]", combo$combo_id, "| merged =", ncol(data), "cells x", nrow(data), "genes\n")

  # ==========================================================
  # Seurat V4 标准流程：按需求固定的参数，其余使用函数默认值
  # ==========================================================
  cat("[STEP]", combo$combo_id, "| NormalizeData\n")
  data <- NormalizeData(data, verbose = FALSE)

  cat("[STEP]", combo$combo_id, "| FindVariableFeatures\n")
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

  cat("[STEP]", combo$combo_id, "| ScaleData\n")
  data <- ScaleData(data, verbose = FALSE)

  cat("[STEP]", combo$combo_id, "| RunPCA\n")
  data <- RunPCA(data, npcs = npcs, verbose = FALSE)

  if(ncol(Embeddings(data[["pca"]])) < max(dims_use)){
    stop(
      paste0(
        "[ERROR] ", combo$combo_id,
        " PCA 可用维度少于 dims_use 最大值 ", max(dims_use), "。"
      ),
      call. = FALSE
    )
  }

  n_harmony_groups <- length(unique(data@meta.data[[harmony_group_col]]))

  if(n_harmony_groups >= 2){
    cat(
      "[STEP]", combo$combo_id,
      "| RunHarmony by", harmony_group_col,
      "| groups =", n_harmony_groups, "\n"
    )
    data <- harmony::RunHarmony(
      data,
      reduction = "pca",
      group.by.vars = harmony_group_col,
      reduction.save = "harmony"
    )
    reduction_use <- "harmony"
  } else {
    # 单 dataset / 单 sample 组合没有可供 Harmony 校正的多批次信息。
    # 此时保留 PCA 作为后续 UMAP/邻居分析的输入，不中断整个组合。
    reduction_use <- "pca"
    cat(
      "[WARN]", combo$combo_id,
      "|", harmony_group_col, "唯一值 =", n_harmony_groups,
      "；跳过 Harmony，使用 PCA 继续分析。\n"
    )
  }

  elbow_file <- file.path(combo_dir, paste0(combo$combo_id, "_elbow.png"))
  save_elbow_plot(data, elbow_file, combo$combo_id)

  cat(
    "[STEP]", combo$combo_id,
    "| RunUMAP | reduction =", reduction_use, "\n"
  )
  data <- RunUMAP(
    data,
    reduction = reduction_use,
    dims = dims_use,
    n.components = 3
  )

  cat(
    "[STEP]", combo$combo_id,
    "| FindNeighbors | reduction =", reduction_use, "\n"
  )
  data <- FindNeighbors(
    data,
    reduction = reduction_use,
    dims = dims_use
  )

  # FindClusters 使用 FindNeighbors 已构建的默认 SNN graph；
  # reduction/dims 设置在 FindNeighbors，而不是 FindClusters。
  cat("[STEP]", combo$combo_id, "| FindClusters resolution =", cluster_resolution, "\n")
  data <- FindClusters(
    data,
    resolution = cluster_resolution
  )

  # 将 3D UMAP 坐标加入 metadata，再保存最终 RDS。
  umap_emb <- Embeddings(data[["umap"]])
  if(ncol(umap_emb) < 3){
    stop(paste0("[ERROR] ", combo$combo_id, " UMAP 不是 3 维。"), call. = FALSE)
  }

  data@meta.data$UMAP_1 <- umap_emb[, 1]
  data@meta.data$UMAP_2 <- umap_emb[, 2]
  data@meta.data$UMAP_3 <- umap_emb[, 3]

  if(!identical(rownames(data@meta.data), colnames(data))){
    stop(paste0("[ERROR] ", combo$combo_id, " metadata 行名与 cell ID 不一致。"), call. = FALSE)
  }

  rds_file <- file.path(combo_dir, paste0(combo$combo_id, ".rds"))
  saveRDS(data, rds_file)
  cat("[SAVE]", rds_file, "\n")

  metadata_csv <- file.path(combo_dir, paste0(combo$combo_id, "_metadata.csv"))
  write.csv(
    data@meta.data,
    metadata_csv,
    row.names = TRUE,
    fileEncoding = "UTF-8"
  )
  cat("[SAVE]", metadata_csv, "\n")

  # UMAP：两个细胞类型标签。
  save_umap_3d_plot(
    data,
    "celltype_label.main",
    file.path(combo_dir, paste0(combo$combo_id, "_umap3d_celltype_main.png")),
    paste0(combo$combo_id, " - celltype_label.main")
  )
  save_umap_3d_plot(
    data,
    "celltype_label.fine",
    file.path(combo_dir, paste0(combo$combo_id, "_umap3d_celltype_fine.png")),
    paste0(combo$combo_id, " - celltype_label.fine")
  )

  # UMAP：组合来源信息。
  save_umap_3d_plot(
    data,
    "dataset_id",
    file.path(combo_dir, paste0(combo$combo_id, "_umap3d_dataset_id.png")),
    paste0(combo$combo_id, " - dataset_id")
  )
  save_umap_3d_plot(
    data,
    "disease_group",
    file.path(combo_dir, paste0(combo$combo_id, "_umap3d_disease_group.png")),
    paste0(combo$combo_id, " - disease_group")
  )
  save_umap_3d_plot(
    data,
    "merge_label",
    file.path(combo_dir, paste0(combo$combo_id, "_umap3d_merge_info.png")),
    paste0(combo$combo_id, " - merge information")
  )

  # UMAP：Seurat 聚类结果，使用同一 3D UMAP 空间展示。
  save_umap_3d_plot(
    data,
    "seurat_clusters",
    file.path(combo_dir, paste0(combo$combo_id, "_umap3d_seurat_clusters.png")),
    paste0(combo$combo_id, " - seurat_clusters")
  )

  final_cells <- ncol(data)
  cat("[DONE]", combo$combo_id, "| final cells =", final_cells, "\n")

  rm(data)
  gc(verbose = FALSE)

  list(
    status = "DONE",
    final_cells = final_cells,
    output_rds = normalizePath(rds_file, winslash = "/", mustWork = FALSE),
    output_metadata_csv = normalizePath(metadata_csv, winslash = "/", mustWork = FALSE),
    output_combination_txt = normalizePath(combination_txt, winslash = "/", mustWork = FALSE),
    error_message = NA_character_
  )
}

# 功能：无有效组合时写出原因说明。
write_no_combination_note <- function(category, catalog, output_file){
  species_need <- if(grepl("^human_", category)) "Human" else "Mouse"
  mode_need <- if(grepl("_disease_in$", category)) "disease_in" else "disease_out"
  df <- catalog[catalog$species == species_need, , drop = FALSE]

  lines <- c(
    paste0("category\t", category),
    paste0("species\t", species_need),
    paste0("combination_type\t", mode_need),
    paste0("eligible_dataset_count\t", nrow(df)),
    paste0(
      "eligible_disease_group_count\t",
      length(unique(df$disease_group))
    ),
    paste0("in_combination_dataset_size\t", in_combination_dataset_size),
    paste0("out_combination_disease_size\t", out_combination_disease_size),
    paste0("out_random_combination_n\t", out_random_combination_n),
    paste0("max_merged_cells\t", max_merged_cells),
    "",
    "未生成有效组合。",
    "组合规则：",
    "1) disease_in：同一个 disease_group 即为同疾病组；",
    "2) disease_in：1 个 dataset 也生成 1 个组合；2-3 个全并；>=4 个生成所有 C(n,4)；",
    "3) disease_out：至少包含两个不同 disease_group；",
    "4) disease_out：<6 个 disease_group 时直接合并全部可用 dataset；",
    "5) disease_out：>=6 个 disease_group 时随机选择最多 20 个不重复 6-group 组合；",
    "6) 每个 6-group 组合合并相应 disease_group 下的全部 dataset；",
    "7) 每个组合总细胞数 >200000 时，在 merge 前按 dataset 比例下采样；",
    "8) Human 与 Mouse 始终分别处理。"
  )

  writeLines(lines, con = output_file, useBytes = TRUE)
  cat("[SAVE]", output_file, "\n")
}

# 功能：执行一个类别的全部组合，持续更新 combination_index.csv。
run_category <- function(category, combos, catalog){
  category_dir <- file.path(output_dir, category)
  dir.create(category_dir, recursive = TRUE, showWarnings = FALSE)

  combos <- assign_combo_ids(combos, category)
  index_df <- build_index_table(combos, catalog)
  index_file <- file.path(category_dir, "combination_index.csv")
  write.csv(index_df, index_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat("[SAVE]", index_file, "\n")

  if(length(combos) == 0){
    note_file <- file.path(category_dir, "NO_VALID_COMBINATION.txt")
    write_no_combination_note(category, catalog, note_file)
    return(list(combos = combos, index = index_df))
  }

  for(i in seq_along(combos)){
    combo <- combos[[i]]
    set.seed(seed + i)

    result <- tryCatch(
      process_one_combination(combo, catalog, category_dir),
      error = function(e){
        cat("[WARN]", combo$combo_id, "处理失败：", conditionMessage(e), "\n")
        list(
          status = "ERROR",
          final_cells = NA_integer_,
          output_rds = NA_character_,
          output_metadata_csv = NA_character_,
          output_combination_txt = normalizePath(
            file.path(category_dir, combo$combo_id, paste0(combo$combo_id, "_combination.txt")),
            winslash = "/",
            mustWork = FALSE
          ),
          error_message = conditionMessage(e)
        )
      }
    )

    idx <- which(index_df$combo_id == combo$combo_id)
    index_df$status[idx] <- result$status
    index_df$final_cells[idx] <- result$final_cells
    index_df$output_rds[idx] <- result$output_rds
    index_df$output_metadata_csv[idx] <- result$output_metadata_csv
    index_df$output_combination_txt[idx] <- result$output_combination_txt
    index_df$error_message[idx] <- result$error_message

    # 每完成一个组合就刷新索引，避免长任务中途停止后丢失进度记录。
    write.csv(index_df, index_file, row.names = FALSE, fileEncoding = "UTF-8")

    if(result$status == "ERROR" && !continue_on_error){
      stop(paste0("[ERROR] ", combo$combo_id, " 失败，continue_on_error=FALSE。"), call. = FALSE)
    }
  }

  list(combos = combos, index = index_df)
}

# 功能：为 human/mouse disease_in 保存每个 disease_group 的 dataset 数与理论组合数。
# 单 dataset disease_group 的 expected_combination_count 现在明确记为 1。
write_in_combination_diagnostics <- function(catalog){
  specs <- list(
    human_disease_in = "Human",
    mouse_disease_in = "Mouse"
  )

  for(category in names(specs)){
    species_need <- specs[[category]]
    df <- catalog[catalog$species == species_need, , drop = FALSE]
    group_names <- sort(unique(df$disease_group))

    if(length(group_names) == 0){
      diag_df <- data.frame(
        species = character(),
        disease_group = character(),
        dataset_count = integer(),
        expected_combination_count = numeric(),
        rule = character(),
        explanation = character(),
        stringsAsFactors = FALSE
      )
    } else {
      rows <- lapply(group_names, function(group_name){
        n_dataset <- sum(df$disease_group == group_name)

        if(n_dataset == 1){
          expected <- 1
          rule <- "single_dataset"
          explanation <- paste0(
            "同一 disease_group 只有 1 个 dataset；按当前要求仍生成 1 个组合并执行分析。"
          )
        } else if(n_dataset < in_combination_dataset_size){
          expected <- 1
          rule <- "merge_all_2to3_datasets"
          explanation <- paste0(
            "同一 disease_group 有 ", n_dataset,
            " 个 dataset，因此全部直接合并并产生 1 个组合。"
          )
        } else {
          expected <- choose(
            n_dataset,
            in_combination_dataset_size
          )
          rule <- "choose_4_datasets_of_n"
          explanation <- paste0(
            "同一 disease_group 有 ", n_dataset,
            " 个 dataset，理论组合数为 C(",
            n_dataset, ",4)=",
            format(expected, scientific = FALSE),
            "。"
          )
        }

        data.frame(
          species = species_need,
          disease_group = group_name,
          dataset_count = n_dataset,
          expected_combination_count = as.numeric(expected),
          rule = rule,
          explanation = explanation,
          stringsAsFactors = FALSE
        )
      })

      diag_df <- do.call(rbind, rows)
    }

    output_file <- file.path(
      output_dir,
      category,
      "in_combination_diagnostics.csv"
    )
    write.csv(
      diag_df,
      output_file,
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    cat("[SAVE]", output_file, "\n")
  }
}

# 功能：运行结束后检查 DONE 记录对应文件和 20 万细胞上限。
post_run_self_check <- function(category_results){
  n_done <- 0L
  n_error <- 0L

  for(category in names(category_results)){
    index_df <- category_results[[category]]$index
    if(nrow(index_df) == 0){
      next
    }

    combo_status <- unique(index_df[, c(
      "combo_id", "status", "final_cells",
      "output_rds", "output_metadata_csv"
    )])

    n_done <- n_done + sum(combo_status$status == "DONE")
    n_error <- n_error + sum(combo_status$status == "ERROR")

    done_rows <- combo_status[combo_status$status == "DONE", , drop = FALSE]
    if(nrow(done_rows) > 0){
      if(any(done_rows$final_cells > max_merged_cells, na.rm = TRUE)){
        bad <- done_rows$combo_id[
          !is.na(done_rows$final_cells) &
            done_rows$final_cells > max_merged_cells
        ]
        stop(
          paste0(
            "[ERROR] 运行后自检发现最终细胞数超过 ",
            max_merged_cells, "：", paste(bad, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      missing_rds <- done_rows$combo_id[!file.exists(done_rows$output_rds)]
      missing_csv <- done_rows$combo_id[!file.exists(done_rows$output_metadata_csv)]

      if(length(missing_rds) > 0){
        warning(
          paste0("[WARN] 自检发现 DONE 但 RDS 缺失：", paste(missing_rds, collapse = ", ")),
          call. = FALSE
        )
      }
      if(length(missing_csv) > 0){
        warning(
          paste0("[WARN] 自检发现 DONE 但 metadata CSV 缺失：", paste(missing_csv, collapse = ", ")),
          call. = FALSE
        )
      }
    }
  }

  cat("[CHECK] 运行后自检：DONE =", n_done, "| ERROR =", n_error, "\n")
  invisible(list(done = n_done, error = n_error))
}

# 功能：保存本次运行参数，保证分析可追踪和可复现。
write_run_config <- function(){
  config_file <- file.path(output_dir, "run_config.txt")

  lines <- c(
    paste0("input_dir\t", input_dir),
    paste0("metadata_file\t", metadata_file),
    paste0("output_dir\t", output_dir),
    paste0("grouping_column\tdisease_group"),
    paste0("in_combination_dataset_size\t", in_combination_dataset_size),
    paste0("out_combination_disease_size\t", out_combination_disease_size),
    paste0("out_random_combination_n\t", out_random_combination_n),
    paste0("combination_enumeration_threshold\t", combination_enumeration_threshold),
    paste0("human_out_sampling_seed\t", get_out_sampling_seed("Human")),
    paste0("mouse_out_sampling_seed\t", get_out_sampling_seed("Mouse")),
    paste0("max_merged_cells\t", max_merged_cells),
    paste0("npcs\t", npcs),
    paste0("dims_use\t", paste(dims_use, collapse = ",")),
    paste0("cluster_resolution\t", cluster_resolution),
    paste0("harmony_group_col\t", harmony_group_col),
    paste0("umap3d_angle\t", umap3d_angle),
    paste0("umap3d_png_width\t", umap3d_png_width),
    paste0("umap3d_png_height\t", umap3d_png_height),
    paste0("umap3d_png_res\t", umap3d_png_res),
    paste0("seed\t", seed),
    paste0("Seurat_version\t", seurat_version),
    paste0("harmony_version\t", as.character(packageVersion("harmony")))
  )

  writeLines(lines, con = config_file, useBytes = TRUE)
  cat("[SAVE]", config_file, "\n")
}

# ============================================================
# 4. 主程序：配置自检、匹配 RDS、读取维度
# ============================================================

cat("[START] Seurat V4 disease combination pipeline\n")
validate_config()
run_internal_logic_check()
write_run_config()

meta_df <- read_and_validate_metadata(metadata_file)

rds_files <- sort(list.files(
  input_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = FALSE
))
if(length(rds_files) == 0){
  stop(paste0("[ERROR] input_dir 下未找到 RDS：", input_dir), call. = FALSE)
}
cat("[CHECK] 发现 RDS 文件：", length(rds_files), "\n")

inventory <- scan_rds_inventory(rds_files, meta_df)
inventory_file <- file.path(output_dir, "sample_inventory.csv")
write.csv(inventory, inventory_file, row.names = FALSE, fileEncoding = "UTF-8")
cat("[SAVE]", inventory_file, "\n")

catalog <- build_eligible_catalog(inventory)
eligible_file <- file.path(output_dir, "eligible_sample_catalog.csv")
write.csv(catalog, eligible_file, row.names = FALSE, fileEncoding = "UTF-8")
cat("[SAVE]", eligible_file, "\n")

cat(
  "[CHECK] 可参与组合的数据：", nrow(catalog),
  "| Human =", sum(catalog$species == "Human"),
  "| Mouse =", sum(catalog$species == "Mouse"), "\n"
)

if(nrow(catalog) == 0){
  stop(
    "[ERROR] 没有可参与组合的 RDS。请检查 RDS 文件名/dataset_id、RNA assay、species 和 disease_group。",
    call. = FALSE
  )
}

# 保存 disease_in 组合数量诊断表，便于检查每个 disease_group 的 dataset 数量和理论组合数。
write_in_combination_diagnostics(catalog)

# ============================================================
# 5. 生成四类组合计划
# ============================================================

category_plans <- list()
for(category in category_dirs){
  cat("[STEP] 生成组合计划：", category, "\n")
  category_plans[[category]] <- plan_category(catalog, category)
  validate_combination_plan(category_plans[[category]], catalog, category)
  cat("[DATA]", category, "| combinations =", length(category_plans[[category]]), "\n")
}

# ============================================================
# 6. 执行四类组合并保存结果
# ============================================================

category_results <- list()
for(category in category_dirs){
  cat("[STEP] 执行类别：", category, "\n")
  category_results[[category]] <- run_category(
    category = category,
    combos = category_plans[[category]],
    catalog = catalog
  )
  gc(verbose = FALSE)
}

# ============================================================
# 7. 汇总索引与运行后自检
# ============================================================

all_index <- do.call(
  rbind,
  lapply(category_results, function(x) x$index)
)
if(is.null(all_index)){
  all_index <- empty_index_template()
}

summary_index_file <- file.path(output_dir, "all_combination_index.csv")
write.csv(all_index, summary_index_file, row.names = FALSE, fileEncoding = "UTF-8")
cat("[SAVE]", summary_index_file, "\n")

post_run_self_check(category_results)

session_file <- file.path(output_dir, "sessionInfo.txt")
capture.output(sessionInfo(), file = session_file)
cat("[SAVE]", session_file, "\n")

cat("[DONE] 所有流程结束。输出目录：", output_dir, "\n")
