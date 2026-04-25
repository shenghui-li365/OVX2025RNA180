library(DESeq2)
library(tidyverse)
library(DESeq2)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(readxl)

# 001.ShamvsOVX 差异分析 -------------------------------------------------------------------
setwd("/share/home/shli24/2025OVX/result/") #mamba activate r_httpgd_env
color_ovx_e2 <- c(
  "Sham"   = "#009E73",
  "OVX"    = "#D55E00",
  "NS" = "#7f7f7f"#,
  #"OVX+E2" = "#0072B2"  # 深蓝（干预 / 恢复）
)
i = "Liver"
aging_hallmarkers_table <- read_xlsx("/share/home/shli24/2025OVX/data/Aging_hallmarkers_GSEA.xlsx")
DE <- read.delim(paste0("DE_result/DEgenename_ShamvsOVX_",i,".txt"))
a <- intersect(toupper(DE$gene_name),aging_hallmarkers_table$Symbol)
DE_sub <- DE[toupper(DE$gene_name) %in% a, ]
marker <- rbind((DE_sub %>%
                   arrange(desc(log2FoldChange)) %>%
                   slice_head(n = 10)),
                DE_sub %>%
                  arrange(log2FoldChange) %>%
                  slice_head(n = 10))
marker <- marker[!is.na(marker$gene_name),]
load("/share/home/shli24/2025OVX/result/DE_result/DEnumber.RData")
res <- read.delim(paste0("DE_result/DE2_ShamvsOVX_",i,".txt"))
res <- res[!is.na(res$direction),]
res$direction <- factor(res$direction ,levels = c('NS','Sham','OVX'))
# table(res$direction, useNA = "ifany")

#    NS  Sham   OVX  <NA>
# 16035   842   819    28
P1 <- ggplot(res,aes(x = log2FoldChange,y = -log10(pvalue))) +
  geom_point(aes(color = direction),size=4) +
  scale_color_manual(name = '',
                     # color or three types
                     values = color_ovx_e2,
                     # legend labels
                     labels = c('OVX'=paste0('OVX (n=',DEnumber_pvalue_ShamvsOVX[2,i],')'),
                               'NS'= paste0('NS (n=',DEnumber_pvalue_ShamvsOVX[1,i],')'),
                               'Sham'=paste0('Sham (n=',DEnumber_pvalue_ShamvsOVX[3,i],')'))) +
    scale_fill_manual(name = '',guide = "none",
                     # color or three types
                     values = color_ovx_e2) +                             
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',linewidth = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',linewidth = 0.8) +
  ggtitle(paste0("Sham vs OVX ",i))+theme_classic()+
  # ylim(-2,7)+
  # xlim(-10,12)+
  geom_text_repel(data = marker,aes(x = log2FoldChange,y = -log10(pvalue),label = marker$gene_name),
                  force=20,color="black",size=4,point.padding = 0.5,hjust = 0.5,
                  arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                  segment.color="black",segment.size=0.2,nudge_y=1)+
  geom_point(data = marker,aes(x = log2FoldChange,y = -log10(pvalue),fill = direction),
             color = "black",shape = 21,size=4)+ 
  theme(
    legend.position = "top",#c(0.4,0.8),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "#ffffff"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black",
                               vjust = 0.5, hjust = 0.5),
    legend.text = element_text(size = 12, color = "black"),
    axis.line.x = element_line(color = "black", size = 0.6),
    axis.line.y = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6) ,axis.ticks.length = unit(2, "mm")     ## 加粗坐标轴刻度线
  )
ggsave(P1, filename =paste0("r_plot/002_DEpvalue_ShamvsOVX_",i,"2.pdf"),width = 4.1, height = 4.6, useDingbats=FALSE)

a <- intersect(toupper(DE$gene_name),aging_hallmarkers_table$Symbol)
DE_sub <- DE[toupper(DE$gene_name) %in% a, ]
DE_sub <- DE_sub[!is.na(DE_sub$gene_name),]
genes <- rbind((DE_sub %>%
                   arrange(desc(log2FoldChange)) %>%
                   slice_head(n = 20)),
                DE_sub %>%
                  arrange(log2FoldChange) %>%
                  slice_head(n = 20))


AA <- c("Hnf4a","Ppara","Foxa1","Foxa2","Foxa3","Alb","Apoa1","Apob","Cyp3a11","Cyp2e1",
"Pck1","G6pc","Fasn","Acaca","Scd1",
"Hoxa1","Hoxa5","Hoxa9","Sox2","Pax6","Gata4","Tbx3", "Ccnb1","Cdk1","Mki67","Pcna","E2f1"
 )
f <- intersect(DE$gene_name,AA)
DE_sub2 <- DE[(DE$gene_name) %in% f, ]
[1] "Ppara"  "Fasn"   "Cyp2e1" "Foxa1"  "Scd1"


# expression heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)

# ================================
# 1. 基因分模块（最终版本）
# ================================
load("/share/home/shli24/2025OVX/result/OVXRNA2.RData")
library(DESeq2)
i=tissue_name= "Liver"
  A_info <- info[info$Group != "OVX+E2" & info$Tissue == tissue_name, ]
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
  rownames(A_info) <- A_info$Name
genename <- read.delim("gene_name.txt")
dds <- DESeqDataSetFromMatrix(countData = A_count,
                             colData = A_info,
                             design = ~ Group)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
module_list <- list(
  Cholesterol = c(
    "Hmgcs1","Hmgcr","Mvk","Pmvk","Mvd","Idi1","Fdps","Sqle","Nsdhl"
  ),
  Lipid_metabolism = c(
    "Fasn","Scd1","Acly","Acss2","Acsl3","Ppara","Foxa1","Insig1","Ppard","Fabp5"
  ),
  Lipid_storage = c(
    "Pnpla3","Mogat1","Mgll","Fitm1","Apoa4","Elovl3"
  ),
  Stress = c(
    "Txnip","Cyp1b1","Cyp2e1","Aldh1a3","Degs2","Gal3st1","Ptges"
  ),
  DNA_inflammation = c(
    "Rad51b","Ung","Mmp7","Thbs4","Inhbe","Rufy4","Tspyl5","Ntrk1"
  )
)
# 建立映射（gene_name → gene_id）
gene_symbols <- unlist(module_list)

# ================================
# 2. 合并 gene + module 映射
# ================================

gene_module <- data.frame(
  gene = unlist(module_list),
  module = rep(names(module_list), lengths(module_list))
)

# 保证顺序（热图按模块排序）
gene_order <- gene_module$gene
gene_ids <- genename[match(gene_order,genename$gene_name),2]


# 如果 mat 没有这些基因顺序
mat_z <- t(scale(t(mat)))
gene_ids2 <- match(gene_ids,rownames(mat_z))
mat2 <- mat_z[gene_ids2, ]
rownames(mat2) <- gene_order 
# 按模块排序
gene_module <- gene_module[match(rownames(mat2), gene_order), ]

# ================================
# 3. 定义模块颜色
# ================================

module_colors <- c(
  Cholesterol       = "#1f78b4",
  Lipid_metabolism  = "#33a02c",
  Lipid_storage     = "#ff7f00",
  Stress            = "#e31a1c",
  DNA_inflammation  = "#6a3d9a"
)

# ================================
# 4. 行注释（module annotation）
# ================================

row_ha <- rowAnnotation(
  Module = gene_module$module,
  col = list(Module = module_colors),
  show_annotation_name = TRUE,
  annotation_legend_param = list(title = "Gene Module")
)

# ================================
# 5. 列注释（Sham / OVX）
# ================================

# 假设列名包含 Sham / OVX 信息
group <- ifelse(grepl("S", colnames(mat2)), "Sham",
         ifelse(grepl("T", colnames(mat2)), "OVX", "Other"))
colnames(mat2) <- c("Sham3","Sham1","Sham2","Sham4","OVX1","OVX2","OVX3")
col_ha <- HeatmapAnnotation(
  Group = group,
  col = list(
    Group = c(
      Sham = "#009E73",
      OVX  = "#D55E00",
      Other = "#999999"
    )
  ),
  annotation_legend_param = list(title = "Group")
)

# ================================
# 6. 颜色函数（表达量）
# ================================

col_fun <- colorRamp2(
  c(min(mat2), 0, max(mat2)),
  c("#2166ac", "white", "#b2182b")
)

# ================================
# 7. 绘制热图
# ================================
pdf("/share/home/shli24/2025OVX/result/r_plot/004_Liver_heatmap.pdf",width = 6.2 ,height = 6.4)
ht <- Heatmap(
  mat2,
  name = "Expression",
  
  col = col_fun,
  
  row_split = gene_module$module,
  
  left_annotation = row_ha,
  top_annotation  = col_ha,
  
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  row_names_gp = gpar(fontface = "italic", fontsize = 10),#bold.italic
  column_names_gp = gpar(fontsize = 10),
  
  row_title_rot = 0,
  row_title_gp = gpar(fontface = "bold"),
  
  column_title = "Samples",
  
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)
dev.off()
library(matrixStats)

# 计算基因方差
gene_var <- rowVars(mat2)

# 取top variable genes
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:30]
mat_top <- mat2[top_genes, ]
col_fun <- colorRamp2(
  c(min(mat_top), 0, max(mat_top)),
  c("#2166ac", "white", "#b2182b")
)

