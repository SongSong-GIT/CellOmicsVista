# 加载包
library(Seurat)
library(tidyverse)
library(ggpubr)
library(scales)

# 读取数据
sc_obj <- readRDS("F:\\老师\\数据库\\celltype2.rds")

# 标准化
# sc_obj <- NormalizeData(sc_obj)
# sc_obj <- ScaleData(sc_obj, features = rownames(sc_obj))

# ===================== 1. 各类状态评分（人类基因） =====================
# 细胞周期基因集
s_genes <- c(
  "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM3",
  "CDT1","GMNN","CDC45","MCM6","DHFR","CDK1","TYMP","MCM7","RRM2","RNASEH2A",
  "H2AFZ","HIST1H2AC","HIST1H2BC","HIST1H2BD","HIST1H2BE","HIST1H2BF",
  "HIST1H2BG","HIST1H2BH","HIST1H2BK","HIST2H2AC"
)
g2m_genes <- c(
  "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
  "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","CKAP2L",
  "CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B",
  "HJURP","CDCA3","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2",
  "DLGAP5","CDCA2","CDCA8","ECT2","KIF23","CENPE","PLK1","AURKA","PSRC1","ANLN"
)
sc_obj <- CellCycleScoring(sc_obj, s.features = s_genes, g2m.features = g2m_genes)

# 干性评分（40个）
stemness_genes <- c(
  "POU5F1","SOX2","NANOG","KLF4","MYC","LIN28A","TDGF1","DPPA5","DPPA4","DPPA2",
  "ZFP42","ESRRB","UTF1","FOXD3","GDF3","FGF4","LEFTY1","LEFTY2","SALL4","ZIC3",
  "MYCN","BMP4","FBXO15","TBX3","TFCP2L1","GRB2","STAT3","LIFR","IL6ST","ETV4",
  "ETV5","FGF5","OTX2","SOX11","MYBL2","PCNA","MKI67","TOP2A","CDK1","AURKB"
)
sc_obj <- AddModuleScore(sc_obj, features = list(stemness = stemness_genes), name = "StemnessScore")

# 免疫活性评分（50个）
immune_genes <- c(
  "CD3D","CD3E","CD3G","CD4","CD8A","CD8B","IFNG","TNF","IL2","IL4","IL5","IL6",
  "IL10","IL12A","IL12B","IL17A","IL21","FOXP3","GZMB","PRF1","NKG7","KLRB1",
  "PTPRC","ITGAM","ITGAX","CD14","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD",
  "HLA-DRA","HLA-DRB1","HLA-DPB1","CD274","PDCD1","CTLA4","LAG3","TIGIT","NFKB1",
  "REL","STAT1","STAT4","STAT6","TBX21","GATA3","RORC","BCL6","ICOS","CD28"
)
sc_obj <- AddModuleScore(sc_obj, features = list(immune = immune_genes), name = "ImmuneScore")

# ===================== 2. 绘图函数 =====================
# 状态-基因表达相关性散点图
plot_correlation <- function(seurat_obj, score_col, gene, title){
  plot_data <- data.frame(
    score = seurat_obj@meta.data[[score_col]],
    expression = as.numeric(seurat_obj@assays$RNA@data[gene, ])
  )
  p <- ggplot(plot_data, aes(x=score, y=expression)) +
    geom_point(size=1, alpha=0.6, color="#3498db") +
    geom_smooth(method="lm", color="#e74c3c", se=TRUE) +
    stat_cor(method="pearson", size=4.5) +
    labs(x="State Score", y=paste0(gene," Expression"), title=title) +
    theme_bw(base_size = 11) + 
    theme(plot.title=element_text(hjust=0.5))
  return(p)
}

# 分组状态评分箱线图
plot_group_boxplot <- function(seurat_obj, score_col, group_col="celltype_label.main", title){
  plot_data <- seurat_obj@meta.data[, c(group_col, score_col)]
  colnames(plot_data) <- c("Group","Score")
  p <- ggplot(plot_data, aes(x=Group, y=Score, fill=Group)) +
    geom_boxplot(alpha=0.7, outlier.size=0.3) + 
    stat_compare_means(method="kruskal.test", size=3.5) +
    scale_fill_brewer(palette="Set2") + 
    labs(x="Cell Type", y="State Score", title=title) +
    theme_bw(base_size = 11) + 
    theme(plot.title=element_text(hjust=0.5), legend.position="none",
          axis.text.x=element_text(angle=45,hjust=1))
  return(p)
}

# ===================== 3. 批量出全量图（补齐细胞周期） =====================
# 散点图：S期、G2M期、干性、免疫
p_s_corr  <- plot_correlation(sc_obj, "S.Score", "PCNA", "S期评分与PCNA表达相关性")
p_g2m_corr<- plot_correlation(sc_obj, "G2M.Score", "MKI67", "G2M期评分与MKI67表达相关性")
p_stem_corr<-plot_correlation(sc_obj, "StemnessScore1", "SOX2", "干性评分与SOX2表达相关性")
p_immune_corr<-plot_correlation(sc_obj, "ImmuneScore1", "CD3D", "免疫活性评分与CD3D表达相关性")

# 箱线图：S期、G2M期、干性、免疫
p_s_box   <- plot_group_boxplot(sc_obj, "S.Score", title = "各细胞类型S期评分分布")
p_g2m_box <- plot_group_boxplot(sc_obj, "G2M.Score", title = "各细胞类型G2M期评分分布")
p_stem_box<- plot_group_boxplot(sc_obj, "StemnessScore1", title = "各细胞类型干性评分分布")
p_immune_box<-plot_group_boxplot(sc_obj, "ImmuneScore1", title = "各细胞类型免疫活性评分分布")

# 统一拼图排版 4行2列
all_plot <- ggarrange(
  p_s_corr, p_g2m_corr,
  p_stem_corr, p_immune_corr,
  p_s_box, p_g2m_box,
  p_stem_box, p_immune_box,
  ncol = 2, nrow = 4
)

# 保存高清图
ggsave("细胞全状态关联分析图.png", all_plot, width = 14, height = 16, dpi = 300)
ggsave("细胞全状态关联分析图.pdf", all_plot, width = 14, height = 16)

# 导出所有评分表格
all_score <- sc_obj@meta.data[,c("S.Score","G2M.Score","StemnessScore1","ImmuneScore1")]
write.csv(all_score, "细胞四大状态评分表.csv", row.names = T)

cat("✅ 已补齐细胞周期S期、G2M期所有图表，全套8张图全部生成完毕！\n")