if(T){
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(SPOTlight)
  library(scran)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(Seurat)
  library(ggplot2)
  library(fastCNV)
  library(irGSEA)
  library(ggpubr)
  library(ggsignif) # 画显著性横线
  library(fs)
  library(R.utils)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(tidyverse)
  library(stringr)
}

datadir = "D:/Desktop/mouse 7.14 - 副本/"
resdir = "D:/Desktop/mouse 7.14 - 副本/res2/"
###########sc###################
if(T){
  library(hdf5r)
  h5_file <- "F:/Download/BRCA_GSE136206_mouse_aPD1aCTLA4_expression.h5"
  mat <- Read10X_h5(filename = h5_file)
  sc <- CreateSeuratObject(
    counts = mat,
    project = "GSE136206",
    min.cells = 3,
    min.features = 200
  )
  meta = rio::import("F:/Download/BRCA_GSE136206_mouse_aPD1aCTLA4_CellMetainfo_table.tsv")
  meta$celltype = meta$`Celltype (major-lineage)`
  
  a = sc@meta.data
  a$Cell = rownames(a)
  a = left_join(a,meta)#keeps all observations in x.
  rownames(a) = a$Cell
  a = a[colnames(sc),]
  sc@meta.data = a
 
  table(sc$celltype)
  sce <- as.SingleCellExperiment(sc)
  sce <- logNormCounts(sce)
  genes <- !grepl(pattern = "^Rp[l|s]|^mt-", x = rownames(sce))
  dec <- modelGeneVar(sce, subset.row = genes)
  # plot(dec$mean, dec$total, xlab ="Mean log-expression", ylab ="Variance")
  #curve(metadata(dec)$trend(x), col ="blue", add =TRUE)
  hvg <- getTopHVGs(dec, n =3000)
  colLabels(sce) <- colData(sce)$celltype# 是否使用celltype要根据自己的数据集进行修改
  mgs <- scoreMarkers(sce, subset.row = genes)
  mgs_fil <- lapply(names(mgs),function(i) {
    x <- mgs[[i]]
    x <- x[x$mean.AUC >0.7, ]
    x <- x[order(x$mean.AUC, decreasing =TRUE), ]
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)
}

if(T){
  celltype_markers <- c(
    "Pecam1", "Mzb1", "Cd79a", "Col1a1","Itgax",
    "Krt8", "Cd68", "Cd3e","Foxp3","Mki67"
  )
  ##### 打分
  # 细胞周期基因
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
  
  # 其他状态基因集
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
  
  
  
  state_gene_sets_m <- lapply(state_gene_sets, function(genes) {
    str_to_title(genes)
  })
  
}
  
  
  
process_second_level <- function(raw_dir){
    files <- dir_ls(raw_dir, recurse = FALSE)
    for(f in files){
      ext <- path_ext(f)
      if(ext == "tar"){
        folder <- file.path(raw_dir, path_ext_remove(fname))
        dir_create(folder)
        try(untar(f, exdir = folder), silent = TRUE)
        next
      }
      if(ext=="gz"){
        if(grepl("(barcodes|features|matrix)\\.(tsv|mtx)\\.gz$", f)){
          message("Skip gunzip: ", f)
          next
        }
        outfile <- file.path(raw_dir, path_ext_remove(path_file(f)))
        try(gunzip(f, destname=outfile, overwrite=TRUE, remove=T), silent=TRUE)
      }
    }
    all_files <- dir_ls(raw_dir, type = "file", recurse = FALSE)
    for(f in all_files){
      fn <- path_file(f)
      sample_prefix <- str_extract(fn, "^GSM\\d+")
      # sample_prefix <-gsub(".rds","",basename(fn))
      if(!is.na(sample_prefix)){
        sample_dir <- file.path(raw_dir, sample_prefix)
        dir_create(sample_dir)
        new_path <- file.path(sample_dir, fn)
        file_move(f, new_path)
      }
    }
  }
  
  
  
organize_visium <- function(root){
    sample_dirs <- dir_ls(root, type = "directory", recurse = FALSE)
    for(sample in sample_dirs){
      files <- dir_ls(sample, recurse = FALSE)
      for(f in files){
        ext <- path_ext(f)
        if(T){
          if(ext == "tar"){
            tmp_dir = dir_create(file.path(sample, paste0("my",gsub(".tar","",basename(f)))))
            untar(f, exdir = tmp_dir)
            contents <- dir_ls(tmp_dir, recurse = TRUE, type = "file")
            if (length(contents) > 0) {
              file_move(contents, sample)
            }
          }
          if(ext == "zip"){
            tmp_dir <- dir_create(file.path(sample, paste0("my",gsub(".zip","",basename(f)))))
            unzip(f, exdir = tmp_dir)

            contents <- dir_ls(tmp_dir,recurse = TRUE, type = "file")
            if (length(contents) > 0) {
              file_move(contents, sample)
            }
          }
        }
      }
      ############解压样本文件夹下的tar并提出
      h5 <- dir_ls(sample, recurse = TRUE, regexp = "filtered_feature_bc_matrix\\.h5$")
      if(length(h5) == 1){
        file_move(h5, file.path(sample, "filtered_feature_bc_matrix.h5"))
      } 
      
      exp_dir <- file.path(sample, "filtered_feature_bc_matrix")
      if(!dir_exists(exp_dir)) dir.create(exp_dir)
      matrix_patterns <- c(
        "barcodes.tsv.gz",
        "features.tsv.gz",
        "matrix.mtx.gz"
      )
      for(target in matrix_patterns){
        found <- dir_ls(sample, recurse = TRUE, regexp = target)
        if(length(found) == 1){
          new_path <- file.path(exp_dir, target)
          file_move(found, new_path)
        }
      }
      
      spatial_dir <- file.path(sample, "spatial")
      if(!dir_exists(spatial_dir)) dir.create(spatial_dir)
      spatial_files <- c("scalefactors_json.json","tissue_hires_image.png",
                         "tissue_lowres_image.png","tissue_positions_list.csv",
                         "tissue_positions.csv")
      pos_csv <- file.path(spatial_dir, "tissue_positions.csv")
      pos_list_csv <- file.path(spatial_dir, "tissue_positions_list.csv")
      if (file.exists(pos_csv)) {
        file_move(pos_csv, pos_list_csv)
        message("已重命名: tissue_positions.csv -> tissue_positions_list.csv")
      }############重命名tissue_positions
      for(target in spatial_files){
        found <- dir_ls(sample, recurse = TRUE, regexp = target)
        if(length(found) == 1){
          target_name <- sub("\\\\", "", target)  # 去掉正则转义符
          new_path <- file.path(spatial_dir, target_name)
          file_move(found, new_path)
          message("Moved spatial file: ", target_name)
        }
      }
    }
  }

############################################begin################################
####################单套数据：解压RAW.tar到TMM文件夹下#################
GSE = "GSE303705" 
if(T){
  raw_tar <- dir_ls(paste0(datadir,GSE),
                    glob = "*_RAW.tar*", type = "file")
  out_dir <- file.path(paste0(datadir,GSE), "TMM")
  dir.create(out_dir,recursive = TRUE,showWarnings = FALSE)
  untar(raw_tar,exdir = out_dir)######解压后删除多余文件
  process_second_level(out_dir)###解压gz分样本建立子文件夹
  organize_visium(out_dir)###建立exp spatial子文件夹
  
  samples <- dir_ls(out_dir, type = "directory")#samples下包含三件套/h5  +   spatial
}


for(i in samples){
  GSM =basename(i)
  # if ("outs" %in% list.dirs(i, full.names = FALSE, recursive = FALSE)) { i <- file.path(i, "outs")}
  if ("filtered_feature_bc_matrix" %in% list.dirs(i, full.names = FALSE, recursive = FALSE)&&
      length(list.files(file.path(i, "filtered_feature_bc_matrix"))) > 0) {
    ##############三件套
    counts <- Read10X(data.dir = paste0(i,"/filtered_feature_bc_matrix"),gene.column = 1)
    ##############转ID
    if(grepl("^ENSMUSG", rownames(counts)[1])){
      ensembl_id <- rownames(counts)
      gene_map <- bitr(
        geneID = ensembl_id,
        fromType = "ENSEMBL",
        toType = "SYMBOL",
        OrgDb = org.Mm.eg.db,  # 换成鼠
        drop = FALSE
      )
      gene_map = na.omit(gene_map[!duplicated(gene_map$SYMBOL),])
      counts = counts[gene_map$ENSEMBL,]
      rownames(counts) = gene_map$SYMBOL
    }
    st <- CreateSeuratObject(counts = counts, assay = "Spatial")
    ##############处理tissue_positions
    pos_file = file.path(file.path(i, "spatial"), "tissue_positions_list.csv")
    if (file.exists(pos_file)) {
      if(grepl("^barcode", readLines( pos_file, n = 1))){
        pos <- read.csv(pos_file, header = TRUE)
        write.table(pos, pos_file, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
      }}
    ##############判断是否存在tissue_lowres_image.png
    if(file.exists(file.path(i, "tissue_lowres_image.png"))){
      image2 <- Read10X_Image(
        image.dir = file.path(i, "spatial"),
        image.name = "tissue_lowres_image.png",  # 明确指定低分辨率图
        filter.matrix = FALSE                   # 先关闭过滤，避免坐标丢失
      )
      DefaultAssay(st = image2) <- "Spatial"
      image2 <- image2[Cells(x = st)]
      st[["slice1"]] <- image2
    }else {
      GSM =basename(i)
      img = Read10X_Image(
        image.dir = file.path(i, "spatial"),
        image.name = "tissue_hires_image.png",  # 明确指定低分辨率图
        filter.matrix = FALSE                   # 先关闭过滤，避免坐标丢失
      )
      st <- CreateSeuratObject(counts = counts, assay = "Spatial")
      DefaultAssay(st = img) <- "Spatial"
      img <- img[Cells(x = st)]
      st[["slice1"]] <- img
      st@images$slice1@scale.factors$lowres = st@images$slice1@scale.factors$hires
    }
  }else if (file.exists(file.path(i, "filtered_feature_bc_matrix.h5"))) {##############标准
    if( file.exists(file.path(i, "spatial/tissue_lowres_image.png"))){
      ##############判断是否存在tissue_lowres_image.png  
      st <- Load10X_Spatial(data.dir = i,
                            filename = "filtered_feature_bc_matrix.h5",
                            image = NULL)
    } else {
      GSM =basename(i)
      img = Read10X_Image(
        image.dir = file.path(i, "spatial"),
        image.name = "tissue_hires_image.png",  # 明确指定低分辨率图
        filter.matrix = FALSE                   # 先关闭过滤，避免坐标丢失
      )
      st <- Load10X_Spatial(data.dir = i,filename = "filtered_feature_bc_matrix.h5",image = img)
      st@images$slice1@scale.factors$lowres = st@images$slice1@scale.factors$hires
    }
  }
  # SpatialFeaturePlot(st,features = "nFeature_Spatial",pt.size.factor = 2,alpha = 1)
  coords <- GetTissueCoordinates(st)
  # colnames(coords)[1:2] = c("x", "y")
  coords <- as.matrix(coords[, c("x", "y")])
  
  
  st <- st[, rownames(coords)]
  spatial <- as.SingleCellExperiment(st)
  spe <- SpatialExperiment(
    assays = assays(spatial),
    colData = colData(spatial),
    spatialCoords = coords
  )
  ################res
  res <- SPOTlight(
    x = sce,
    y = spe,
    groups = as.character(sce$celltype),
    mgs = mgs_df,
    hvg = hvg,
    weight_id ="mean.AUC",
    group_id ="cluster",
    gene_id ="gene")
  
  mat = res$mat
  coords = as.data.frame(coords)
  coords$sample = basename(i)
  coords$cancer = "mouse"
  
  irGSEAres <- irGSEA.score(
    object = st, assay = "Spatial",
    slot = "counts", 
    seeds = 123, 
    min.cells = 0, min.feature = 0,
    custom = T, ncores = 1 ,
    geneset = state_gene_sets_m,
    msigdb = FALSE,
    species = "Mus_musculus", geneid = "symbol",
    method = "AUCell",
    kcdf = 'Poisson',
    minGSSize = 1, maxGSSize = 2000)
  if ("layers" %in% slotNames(irGSEAres@assays$AUCell)) {
    # Seurat v5: Assay5 结构
    score <- data.frame(t(irGSEAres@assays$AUCell@layers$scale.data))
  } else {
    # Seurat v4: 普通 Assay 结构
    score <- data.frame(t(irGSEAres@assays$AUCell@scale.data))
  }
  rownames(score) <- colnames(st)
  colnames(score) <- names(state_gene_sets_m)
  
  

  st <- NormalizeData(st, normalization.method="LogNormalize")
  marker  = data.frame(t(st@assays$Spatial$data[celltype_markers,]))
  ####################整理表格
  res  =cbind(coords,mat)
  res$cellid = rownames(res)
  res  =cbind(res,marker)
  res  =cbind(res,score)
  
  res_long <- res %>%
    pivot_longer(
      cols = all_of(celltype_markers),
      names_to = "gene",
      values_to = "expr"
    )
  rm(st)
  gc()
  # GSM=basename(i)
  rio::export(res_long,paste0(resdir,GSE,"_",GSM,"_res.xlsx"))
}

######绘图...细胞注释######
values <- c(
  "B"              = "#D4FA38",
  "CD4Tconv"       = "#47A857",
  "CD8T"           = "#1DE07A",
  "CD8Tex"         = "#45C7A8",
  "DC"             = "#3CC0C0",
  "Endothelial"    = "#E8A080",
  "Fibroblasts"    = "#D988D9",
  "Malignant"      = "#F25A72",
  "Mono/Macro"     = "#66B2B2",
  "Neutrophils"    = "#3A96BB",
  "NK"             = "#8CBFAF",
  "Plasma"         = "#D6DB70",
  "Tprolif"        = "#52E01F"
)


st$celltype = colnames(mat)[max.col(mat, ties.method = "first")]
table(st$celltype)

p1 = SpatialDimPlot(st,group.by = "celltype",pt.size.factor = 4)+
  scale_fill_manual(values = values)

#####绘图...基因表达可视化####
p2 = SpatialFeaturePlot(st,features =celltype_markers,pt.size.factor = 2,
                   image.alpha = 0.5,ncol = 5)

res$celltype = st$celltype
res_long <- res %>%
  pivot_longer(
    cols = all_of(celltype_markers),
    names_to = "gene",
    values_to = "expr"
  )
p3 = ggplot(coords_long, aes(x = celltype, y = expression, fill = celltype)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = values) +
  facet_wrap(~gene, ncol = 5, scales = "free_y") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 1,
      hjust = 1,
      size = 9
    ),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(x = "Cell Type", y = "Expression")



st@meta.data = cbind(st@meta.data,score)
# 2. 自动找出所有 state_gene_sets 开头的模块分数
plot_features <- names(state_gene_sets)
plots <- lapply(names(state_gene_sets), function(gene) {
  SpatialFeaturePlot(
    st,
    features = gene,
    pt.size.factor = 2,
    image.alpha = 0.5
  ) +
    NoLegend() +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
})

# 拼接成一张大图（自动按列排布）
combined_plot <- wrap_plots(plots, ncol = 4)  # ncol 改列数


(p1+p2)/p3/combined_plot
