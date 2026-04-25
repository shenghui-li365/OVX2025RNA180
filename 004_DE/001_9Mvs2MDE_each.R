# 001.2Mvs9M 差异分析 -------------------------------------------------------------------
#module load arm/r/4.4.1
#R
setwd("/share/home/shli24/2025OVX/result/") #mamba activate r_httpgd_env
color_age <- c(
  "NS" = "#7f7f7f",
  "2M" = "#009E73",   # young
  "9M" = "#CC79A7"    # aged（muted purple）
)
load("/share/home/shli24/2025OVX/result/agingRNA.RData")
colnames(count)==colnames(tpm)
colnames(tpm)
colnames(count)
info$Name
library(DESeq2)
library(tidyverse)
unique(info$Tissue)
colnames(info)
info$Name
i = "Spleen"
library(DESeq2)
library(dplyr)
library(ggrepel)
aa <- read.delim("DE_result/DE_ShamvsOVX_Liver.txt")
base_df <- as.data.frame(table(aa$direction))
colnames(base_df) <- c("direction", "Count")
  DEnumber_pvalue_2Mvs9M = base_df
  DEnumber_pvalue_52Mvs9M = base_df
  DEnumber_pvalue1_2Mvs9M = base_df
  DEnumber_pvalue1_52Mvs9M = base_df
  DEnumber_pvalue2_2Mvs9M = base_df
  DEnumber_padj_2Mvs9M = base_df
  DEnumber_padj_52Mvs9M = base_df
  DEnumber_padj1_2Mvs9M = base_df
  DEnumber_padj1_52Mvs9M = base_df

#定义函数
#000 define function---------------------------------------------------------
run_DESeq2_by_tissue <- function(tissue_name, info, count, tpm) {
  # 1. 提取当前组织的数据
  DE_info   <- info
  DE_count  <- count
  # 2. 对齐样本名
  rownames(DE_info) <- DE_info$Name
  colnames(DE_count) <- DE_info$Name[colnames(DE_count) %in% DE_info$Name]
  DE_count <- round(DE_count)
  keep <- rowSums(DE_count >= 10) >= 2
  DE_count <- DE_count[keep, ]
  #DE_count <- DE_count[rowSums(DE_count) >= 2, ]
  DE_tpm   <- tpm[, DE_info$Name]
  
  # 3. 构建 DESeq2 对象
  conditions <- data.frame(conditions = factor(DE_info$Group))
  rownames(conditions) <- colnames(DE_count)
  dds <- DESeqDataSetFromMatrix(countData = DE_count,
                                colData = conditions,
                                design = ~ conditions)
  dds <- DESeq(dds)
  
  # 4. 提取结果
  res <- results(dds, contrast = list(c("conditions_9M_vs_2M"))) %>%
    as.data.frame()
  
  # 5. 添加各种方向标签（基于 pvalue 和 padj，不同 log2FC 阈值）
  # 注意：先基于 pvalue，后基于 padj，注意变量名不能重复
  res <- res %>%
    mutate(
      # 基于 pvalue
      direction = if_else(pvalue > 0.05, "NS",
                          if_else(log2FoldChange > 0, "9M",
                                  if_else(log2FoldChange < 0, "2M", "NS"))),
      direction_5 = if_else(pvalue > 0.05, "NS",
                            if_else(log2FoldChange > 0.5, "9M",
                                    if_else(log2FoldChange < -0.5, "2M", "NS"))),
      direction1 = if_else(pvalue > 0.05, "NS",
                           if_else(log2FoldChange > 1, "9M",
                                   if_else(log2FoldChange < -1, "2M", "NS"))),
      direction1_5 = if_else(pvalue > 0.05, "NS",
                             if_else(log2FoldChange > 1.5, "9M",
                                     if_else(log2FoldChange < -1.5, "2M", "NS"))),
      direction2 = if_else(pvalue > 0.05, "NS",
                           if_else(log2FoldChange > 2, "9M",
                                   if_else(log2FoldChange < -2, "2M", "NS"))),
      # 基于 padj
      Direction = if_else(padj > 0.05, "NS",
                          if_else(log2FoldChange > 0, "9M",
                                  if_else(log2FoldChange < 0, "2M", "NS"))),
      Direction_5 = if_else(padj > 0.05, "NS",
                            if_else(log2FoldChange > 0.5, "9M",
                                    if_else(log2FoldChange < -0.5, "2M", "NS"))),
      Direction1 = if_else(padj > 0.05, "NS",
                           if_else(log2FoldChange > 1, "9M",
                                   if_else(log2FoldChange < -1, "2M", "NS"))),
      Direction1_5 = if_else(padj > 0.05, "NS",
                             if_else(log2FoldChange > 1.5, "9M",
                                     if_else(log2FoldChange < -1.5, "2M", "NS")))
    )
   A <- table(res$direction)
   A <- A[c("NS", "9M","2M")]
   A[is.na(A)] <- 0
   # 转成列写入
   DEnumber_pvalue_2Mvs9M[, i] <- A

   A <- table(res$direction_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue_52Mvs9M[, i] <- A

   A <- table(res$direction1)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue1_2Mvs9M[, i] <- A
   A <- table(res$direction1_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue1_52Mvs9M[, i] <- A

  A <- table(res$direction2)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue2_2Mvs9M[, i] <- A
A <- table(res$Direction)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj_2Mvs9M[, i] <- A

   A <- table(res$Direction_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj_52Mvs9M[, i] <- A

  A <- table(res$Direction1)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj1_2Mvs9M[, i] <- A
  A <- table(res$Direction1_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj1_52Mvs9M[, i] <- A
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
#001 "Liver"---------------------------------------------------------
i=tissue_name= "Liver"
  A_info <- info[info$Tissue == tissue_name, ][c(-3,-7),]#O3 Y5
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Sample
#002 "Lung"---------------------------------------------------------
i=tissue_name= "Lung" #可加测Y5
  A_info <- info[info$Tissue == tissue_name, ][c(-4,-7),]#O4 Y9
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Sample
#003 "Kidney"---------------------------------------------------------
i=tissue_name= "Kidney"
  A_info <- info[info$Tissue == tissue_name, ][c(-3,-6),]#O3 Y4
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Sample
#004 "Heart"---------------------------------------------------------
i=tissue_name= "Heart"
  A_info <- info[info$Tissue == tissue_name, ][c(-2,-7),]#O3 Y5
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Sample
#005 "Spleen"---------------------------------------------------------
i=tissue_name= "Spleen"
  A_info <- info[info$Tissue == tissue_name, ][c(-4),]# Y3
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Sample
info$SampleName[c(9:11,17,24,30,32,38)]
info <-info[c(-9:-11,-17,-24,-30,-32,-38),]
count <- count[, info$SampleName]
tpm   <- tpm[, info$SampleName]
save(info,count,tpm,file = "/share/home/shli24/2025OVX/result/agingRNA2.RData")
