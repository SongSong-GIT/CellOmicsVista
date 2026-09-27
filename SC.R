#############################    聚类画图            ############################################
SCrna_seq_analyse=function(seurat.path,
                           humanORmouse,
                           run_rm_CellCycle=FALSE,
                           nfeatures=3000,  ## FindVariableFeatures
                           pc.contribution=90,
                           my_ndims=FALSE,
                           ndims.1=NULL,  ## 参与降维聚类的PC数
                           ndims.2=NULL,  ## 参与去双细胞的PC数，默认和ndims.1一致
                           run_rm_doublet=TRUE,  ## 是否需要去双细胞
                           metadata_sample,  ## 样本列名
                           run_RunHarmony=TRUE,  ## 是否需要去批次
                           run_decontX=TRUE, ## 是否需要去环境污染细胞
                           decontX_cutoff=0.2, ## 去除环境污染阈值
                           run_FindClusters=TRUE ## FindClusters参数
){
  
  require(Seurat)
  require(dplyr)
  require(rio)
  require(DoubletFinder)
  require(decontX)
  require(clustree)
  
  
  rm_doublet <- function(input=NULL,dim.usage=30,poi.dfr=0.075,auto="true") {
    if(!requireNamespace("DoubletFinder", quietly = TRUE)){
      devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
    }
    library(DoubletFinder)
    if(class(input)=="character") EC <- readRDS(input)
    if(class(input)=="Seurat") EC <- input
    
    Find_doublet <- function(data){
      sweep.res.list <- paramSweep(data, PCs = 1:dim.usage, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      nExp_poi <- round(poi.dfr*ncol(data))
      p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
      data <- doubletFinder(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
      colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
      return(data)
    }
    
    EC<-Find_doublet(EC)
    EC<-subset(EC,subset=doublet_info=="Singlet")
    c <- grep("pANN_",colnames(EC@meta.data))
    EC@meta.data <- EC@meta.data[,-c]
    return(EC)
  }
  
  
  ## 参数调整
  if(class(seurat.path)=="character"){
    da=rio::import(seurat.path)
  }else{
    da=seurat.path
  }
  if(nfeatures==3000){
    nfeatures
  }else{
    nfeatures=nfeatures
  }
  if(my_ndims){
    if(!is.null(ndims.1)){
      ndims.1=30
    }else{
      ndims.1=ndims.1
    }
    if(!is.null(ndims.2)){
      ndims.2=30
    }else{
      ndims.2=ndims.2
    }
  }
  if(decontX_cutoff==0.2){
    decontX_cutoff=0.2
  }else{
    decontX_cutoff=decontX_cutoff
  }
  if(pc.contribution==90){
    pc.contribution=90
  }else{
    pc.contribution=pc.contribution
  }
  
  ## 分析流程
  da <- NormalizeData(da)
  da <- FindVariableFeatures(da, selection.method = "vst", nfeatures = nfeatures)
  if(run_rm_CellCycle){
    if(humanORmouse=="human"){
      s.genes=Seurat::cc.genes.updated.2019$s.genes
      g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
    }
    if(humanORmouse=="mouse"){
      s.genes=c("Mcm5","Pcna","Tyms","Fen1","Mcm7","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7",
                "Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Polr1b","Nasp","Rad51ap1","Gmnn","Wdr76",
                "Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45,Cdc6","Exo1",
                "Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Mrpl36","E2f8")
      g2m.genes=c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67",
                  "Tmpo","Cenpf","Tacc3","Pimreg","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
                  "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2",
                  "Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe",
                  "Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
    }
    da <- CellCycleScoring(da, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    da <- ScaleData(da, vars.to.regress = c("S.Score", "G2M.Score"))
  }else{
    da <- ScaleData(da)
  }
  da <- RunPCA(da)
  if(!my_ndims){
    #确定每个pc的百分比，主成分累计贡献大于90%，选择拐点
    pct<- da[["pca"]]@stdev/sum(da[["pca"]]@stdev)*100
    cumu<-cumsum(pct)
    # 计算每个PC的百分比
    ndims.1 <- min(which(cumu >pc.contribution & pct< 5)[1],
                   sort(which((pct[1:length(pct)-1]- pct[2:length(pct)])>0.1),
                        decreasing=T)[1]+1)
  }
  pca.ndims=ndims.1
  if(run_rm_doublet){
    da <- RunUMAP(da,dims = 1:ndims.1)
    ##### 去除双细胞
    da <- SplitObject(da,split.by = metadata_sample)
    da <- lapply(da,function(x){
      rm_doublet(input=x,dim.usage=ndims.1,poi.dfr=0.075,auto="true")
    })
    da <- merge(da[[1]],
                y = c(lapply(2:length(da),function(x)da[[x]]) %>% unlist()),
                project = "scRNA");dim(da)
    da <- NormalizeData(da)
    da <- FindVariableFeatures(da, selection.method = "vst", nfeatures = nfeatures)
    da <- ScaleData(da)
    da <- RunPCA(da)
  }
  if(run_RunHarmony){
    library(harmony)
    da <- RunHarmony(da,reduction = "pca",
                     group.by.vars = metadata_sample,
                     reduction.save = "harmony")
  }
  if(run_RunHarmony){
    if(!my_ndims){
      pct<- da[["harmony"]]@stdev/sum(da[["harmony"]]@stdev)*100
      cumu<-cumsum(pct)
      # 计算每个PC的百分比
      ndims.2<- min(which(cumu >pc.contribution & pct< 5)[1],sort(which((pct[1:length(pct)-1]- pct[2:length(pct)])>0.1),decreasing=T)[1]+1)
    }else{
      ndims.2=ndims.2
    }
    da <- FindNeighbors(da, dims = 1:ndims.2,reduction = "harmony")
  }else{
    if(!my_ndims){
      pct<- da[["harmony"]]@stdev/sum(da[["harmony"]]@stdev)*100
      cumu<-cumsum(pct)
      # 计算每个PC的百分比
      ndims.2<- min(which(cumu >pc.contribution & pct< 5)[1],sort(which((pct[1:length(pct)-1]- pct[2:length(pct)])>0.1),decreasing=T)[1]+1)
    }else{
      ndims.2=ndims.2
    }
    da <- FindNeighbors(da, dims = 1:ndims.2)
  }
  Harmony_ndims=ndims.2
  if(isTRUE(run_FindClusters)){
    da <- FindClusters(da, resolution = 0.3)
    da <- FindClusters(da, resolution = 0.4)
    da <- FindClusters(da, resolution = 0.5)
    da <- FindClusters(da, resolution = 0.6)
    da <- FindClusters(da, resolution = 0.7)
    da <- FindClusters(da, resolution = 0.8)
    da <- FindClusters(da, resolution = 0.9)
    da <- FindClusters(da, resolution = 1)
    da <- FindClusters(da, resolution = 1.2)
    
    # 聚类树
    clus.tree.out<-clustree(da)+
      theme(legend.position = "bottom") +
      scale_color_brewer(palette = "Set1") +
      scale_edge_color_continuous(low = "grey80", high = "red")
  }
  if(!isTRUE(run_FindClusters)){
    da <- FindClusters(da, resolution = run_FindClusters)
    clus.tree.out =NULL
  }
  if(run_RunHarmony){
    da <- RunUMAP(da, dims = 1:ndims.2,reduction = "harmony",reduction.name = "umap")
    da <- RunTSNE(da, dims = 1:ndims.2,reduction = "harmony",reduction.name = "tsne")
  }else{
    da <- RunUMAP(da, dims = 1:ndims.1)
    da <- RunTSNE(da, dims = 1:ndims.1)
  }
  return(list(result=da,pca.ndims=pca.ndims,Harmony_ndims=Harmony_ndims,clus.tree.out=clus.tree.out))
}

gc()
#调用函数
results <- SCrna_seq_analyse(
  seurat.path =scRNA_ha,       # 数据路径或Seurat对象
  humanORmouse = "human",          # 设置物种   human
  run_rm_CellCycle = F,         # 是否去除细胞周期影响
  nfeatures = 3000,                # 高变基因数量
  pc.contribution = 90,            # 主成分累计贡献阈值
  my_ndims = T,                # 自动确定PC数量
  ndims.1=30,                    # 参与降维聚类的PC数
  ndims.2=30,                    #参与去双细胞的PC数，默认和ndims.1一致
  run_rm_doublet = F,           # 是否去除双细胞
  metadata_sample = "orig.ident",  # 样本名称
  run_RunHarmony = F,           # 是否去批次
  run_decontX = F,              # 是否去环境污染
  decontX_cutoff = 0.2,            # 环境污染阈值
  run_FindClusters = TRUE          # 聚类分辨率
)

# 提取结果
sce <- results$result
clus_tree_plot <- results$clus.tree.out
if (!is.null(clus_tree_plot)) {
  print(clus_tree_plot)
}
head(sce@meta.data)

DimPlot(sce,group.by = "RNA_snn_res.0.5",label = T)
DimPlot(seurat_object,group.by = "orig.ident",label = T)
markers <- c("EPCAM","SFTPC","CC10",#Epithelial
             "PECAM1","CDH5","CD31",#Endothelial
             "PDGFRA","COL1A1","TNC",#Fibroblast
             "CD68","CD64","CD11C",#Macrophag
             "CD3D","CD3E",  #T cells
             "CD79A","MS4A1","CD79B",  #B/plasma cells
             "S100A8","S100A9","LY6G" ,#Neutrophil
             "CD14","FCN1" ,"CCR2","LY6C","CD11B" , #Monocyte
             "TPSAB1","MS4A2", #Mast cell
             "FCGR2B","CYBB","CD53","CTSS","FCER1G" ,#Myeloid cell
             "KRT5","KRT17",#Basal cell	
             "NKG7", #nk
             "PDGFRB","RGS5","ACTA2" #MCs
)
da <- seurat_object
DotPlot(da, features = markers,
        assay='RNA' ,group.by = 'RNA_snn_res.0.6' )  + coord_flip()

da@meta.data$celltype <- NA
da@meta.data$celltype[da@meta.data$RNA_snn_res.0.6 %in% c(3,5,6,10,14)] <- 'Mural cells'
da@meta.data$celltype[da@meta.data$RNA_snn_res.0.6 %in% c(2,4,11,9,12)] <- 'Endothelial'
da@meta.data$celltype[da@meta.data$RNA_snn_res.0.6 %in% c(15,16,17,8)] <- 'Myeloid cell'
da@meta.data$celltype[da@meta.data$RNA_snn_res.0.6 %in% c(0,1,13,7,19)] <- 'Fibroblast'
da@meta.data$celltype[da@meta.data$RNA_snn_res.0.6 %in% c(18)] <- 'Mast cell'
DimPlot(da, group.by = "celltype", label = T)


