library(Seurat)
library(hdf5r)
library(tidyverse)
library(rio)
library(dplyr)
library(clustree)
library(ggplot2)
library(SingleR)
setwd("F:\\老师\\数据库")
da=readRDS("F:\\单细胞及服务器\\单细胞注释\\CHOL_GSE125449\\pbmc.rds")
dim(da)
head(da@meta.data)

SC_Autoanno_SingleR=function(seurat_obj,
                             module,
                             seurat_clusters){
  ## SingleR自动注释
  require(SingleR)
  require(Seurat)
  data_in <- GetAssayData(seurat_obj, slot="data")
  clusters <- seurat_obj@meta.data[,seurat_clusters]
  load(module)
  
  ###  自动注释
  for(methed in c("label.main","label.fine")){
    cell_singleR <- SingleR(test = data_in, ref = ref, labels = ref[[methed]],
                            method = "cluster", clusters = clusters,
                            assay.type.test = "logcounts", assay.type.ref = "logcounts")
    celltype = data.frame(ClusterID=rownames(cell_singleR),
                          celltype=cell_singleR$labels,
                          stringsAsFactors = FALSE)
    seurat_obj@meta.data[,paste0("celltype_",methed)] = "NA"
    for(i in 1:nrow(celltype)){
      seurat_obj@meta.data[which(seurat_obj@meta.data[,seurat_clusters] == celltype$ClusterID[i]),
                           paste0("celltype_",methed)] <- celltype$celltype[i]
    }
    message(paste0("SingleR ",methed," running OKK..."))
  }
  return(seurat_obj)
}
seurat_annotated <- SC_Autoanno_SingleR(
  seurat_obj = da,
  module = "F:\\单细胞及服务器\\单细胞注释\\celltype_anno\\celltype_anno\\singleR\\HumanPrimaryCellAtlas_hpca.se_human.RData",#换成自己的路径
  seurat_clusters = "RNA_snn_res.0.3"
)
DimPlot(seurat_annotated, group.by = "celltype_label.main", label = T)
saveRDS(da,"celltype2.rds")
#######################     细胞通讯             ##############################

library(CellChat)

# 提取数据
data <- GetAssayData(seurat_annotated, slot = "data")
meta <- seurat_annotated@meta.data
cellchat <- createCellChat(object = data, meta = meta, group.by = "celltype_label.main")

# 数据库
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# 分析
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, color.heatmap = "Reds")
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:7))#有多少细胞类型换为多少
###################################################
# # 细胞周期评分
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# seurat_annotated <- CellCycleScoring(
#   seurat_annotated,
#   s.features = s.genes,
#   g2m.features = g2m.genes
# )
# 
# # 免疫活性评分 (你可自定义基因集)
# immune_gene <- c("CD2","CD3D","CD3E","CD8A","GZMA","PRF1","IFNG")
# seurat_annotated <- AddModuleScore(
#   seurat_annotated,
#   features = list(immune_gene),
#   name = "ImmuneScore"
# )
# 
# # 干性评分
# stemness_gene <- c("POU5F1","SOX2","NANOG","KLF4","MYC")
# seurat_annotated <- AddModuleScore(
#   seurat_annotated,
#   features = list(stemness_gene),
#   name = "StemScore"
# )