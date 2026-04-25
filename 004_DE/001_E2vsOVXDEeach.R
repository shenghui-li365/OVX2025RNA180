# 001.ShamvsOVX 差异分析 -------------------------------------------------------------------
setwd("/share/home/shli24/2025OVX/result/") #mamba activate r_httpgd_env
color_ovx_e2 <- c(
  "Sham"   = "#009E73",
  "OVX"    = "#D55E00",
  "NS" = "#7f7f7f",
  "OVX+E2" = "#0072B2"  # 深蓝（干预 / 恢复）
)
load("/share/home/shli24/2025OVX/result/OVXRNA.RData")
colnames(count)==colnames(tpm)
colnames(tpm)
colnames(count)
info$Name
library(DESeq2)
library(tidyverse)
unique(info$Tissue)
colnames(info)
info$Name
i=tissue_name = "Spleen"
library(DESeq2)
library(dplyr)
#定义函数
#001 define function---------------------------------------------------------
run_DESeq2_by_tissue <- function(tissue_name, info, count, tpm) {
  # 1. 提取当前组织的数据
  DE_info   <- info
  DE_count  <- count
  # 2. 对齐样本名
  rownames(DE_info) <- DE_info$Name
  colnames(DE_count) <- DE_info$Name[colnames(DE_count) %in% DE_info$Name]
  DE_count <- round(DE_count)
  DE_count <- DE_count[rowSums(DE_count) >= 2, ]
  DE_tpm   <- tpm[, DE_info$Name]
  
  # 3. 构建 DESeq2 对象
  conditions <- data.frame(conditions = factor(DE_info$Group))
  rownames(conditions) <- colnames(DE_count)
  dds <- DESeqDataSetFromMatrix(countData = DE_count,
                                colData = conditions,
                                design = ~ conditions)
  dds <- DESeq(dds)
  
  # 4. 提取结果
  res <- results(dds, contrast = list(c("conditions_OVX.E2_vs_OVX"))) %>%
    as.data.frame()
  
  # 5. 添加各种方向标签（基于 pvalue 和 padj，不同 log2FC 阈值）
  # 注意：先基于 pvalue，后基于 padj，注意变量名不能重复
  res <- res %>%
    mutate(
      # 基于 pvalue
      direction = if_else(pvalue > 0.05, "NS",
                          if_else(log2FoldChange > 0, "OVX+E2",
                                  if_else(log2FoldChange < 0, "OVX", "NS"))),
      direction_5 = if_else(pvalue > 0.05, "NS",
                            if_else(log2FoldChange > 0.5, "OVX+E2",
                                    if_else(log2FoldChange < -0.5, "OVX", "NS"))),
      direction1 = if_else(pvalue > 0.05, "NS",
                           if_else(log2FoldChange > 1, "OVX+E2",
                                   if_else(log2FoldChange < -1, "OVX", "NS"))),
      direction1_5 = if_else(pvalue > 0.05, "NS",
                             if_else(log2FoldChange > 1.5, "OVX+E2",
                                     if_else(log2FoldChange < -1.5, "OVX", "NS"))),
      direction2 = if_else(pvalue > 0.05, "NS",
                           if_else(log2FoldChange > 2, "OVX+E2",
                                   if_else(log2FoldChange < -2, "OVX", "NS"))),
      # 基于 padj
      Direction = if_else(padj > 0.05, "NS",
                          if_else(log2FoldChange > 0, "OVX+E2",
                                  if_else(log2FoldChange < 0, "OVX", "NS"))),
      Direction_5 = if_else(padj > 0.05, "NS",
                            if_else(log2FoldChange > 0.5, "OVX+E2",
                                    if_else(log2FoldChange < -0.5, "OVX", "NS"))),
      Direction1 = if_else(padj > 0.05, "NS",
                           if_else(log2FoldChange > 1, "OVX+E2",
                                   if_else(log2FoldChange < -1, "OVX", "NS"))),
      Direction1_5 = if_else(padj > 0.05, "NS",
                             if_else(log2FoldChange > 1.5, "OVX+E2",
                                     if_else(log2FoldChange < -1.5, "OVX", "NS")))
    )
  
  # 6. 可选：输出每个方向的计数（或直接返回 res）
  cat("\n=== Tissue:", tissue_name, "===\n")
  cat("基于 pvalue (direction):\n")
  print(table(res$direction))
  cat("基于 pvalue 且 log2FC > 0.5:\n")
  print(table(res$direction_5))
  cat("基于 pvalue 且 log2FC > 1:\n")
  print(table(res$direction1))
  cat("基于 pvalue 且 log2FC > 1.5:\n")
  print(table(res$direction1_5))
  cat("基于 pvalue 且 log2FC > 2:\n")
  print(table(res$direction2))
  cat("基于 padj (Direction):\n")
  print(table(res$Direction))
  cat("基于 padj 且 log2FC > 0.5:\n")
  print(table(res$Direction_5))
  cat("基于 padj 且 log2FC > 1:\n")
  print(table(res$Direction1))
  cat("基于 padj 且 log2FC > 1.5:\n")
  print(table(res$Direction1_5))
  
  # 返回结果对象，以便后续保存或进一步分析
  return(res)
}
#002 "Liver"---------------------------------------------------------
i=tissue_name= "Liver"
  A_info <- info[info$Group != "Sham" & info$Tissue == tissue_name, ][c(-3,-6),] #T52
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
colnames()
#002 "Kidney"---------------------------------------------------------
i= tissue_name ="Kidney"
  A_info <- info[info$Group != "Sham" & info$Tissue == tissue_name, ][c(-1,-8),] #S72 T79
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Name
#003 "Heart"---------------------------------------------------------
i= tissue_name ="Heart"
  A_info <- info[info$Group != "Sham" & info$Tissue == tissue_name, ][c(-3),] #S85
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Name
#004 "Lung"---------------------------------------------------------
i= tissue_name ="Lung"
  A_info <- info[info$Group != "Sham" & info$Tissue == tissue_name, ][c(-4,-8),] #S72 T79
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Name
#005 "Spleen"---------------------------------------------------------
i= tissue_name ="Spleen"
  A_info <- info[info$Group != "Sham" & info$Tissue == tissue_name, ][c(-6),] #T52
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Name