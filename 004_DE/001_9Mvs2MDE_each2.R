# 001.2Mvs9M 差异分析 -------------------------------------------------------------------
#module load arm/r/4.4.1
#R
setwd("/share/home/shli24/2025OVX/result/") #mamba activate r_httpgd_env
color_age <- c(
  "NS" = "#7f7f7f",
  "2M" = "#009E73",   # young
  "9M" = "#CC79A7"    # aged（muted purple）
)
load("/share/home/shli24/2025OVX/result/agingRNA2.RData")
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
base_df$direction <- c("NS", "9M","2M")
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
  A_info <- info[info$Tissue == tissue_name, ]
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)
A_info$Sample
#002 "batch"---------------------------------------------------------
genename <- read.delim("gene_name.txt")

DEresult <- list()
DEgenenameresult <- list()

n = 2
for (i in unique(info$Tissue)){
  A_info <- info[info$Tissue == i, ]
  A_count <- count[, A_info$Name]
  A_tpm   <- tpm[, A_info$Name]
  res   <- run_DESeq2_by_tissue(i, A_info, A_count, A_tpm)

write.table(res,paste0("DE_result/DE2_9Mvs2M_",i,".txt"),sep = "\t",quote = F,row.names = T,col.names = T)
res$gene_id <-rownames(res)
DEresult[[i]] <- res
res <-res[!is.na(res$direction),]
DE <- res[which(res$direction!="NS"),]
DE <-merge(DE,genename,by="gene_id", all.x = FALSE)
DEgenenameresult[[i]] <- DE
marker <- rbind((DE %>%
                   arrange(desc(log2FoldChange)) %>%
                   slice_head(n = 10)),
                DE %>%
                  arrange(log2FoldChange) %>%
                  slice_head(n = 10))
marker <- marker[!is.na(marker$gene_name),]

P1 <- ggplot(res,aes(x = log2FoldChange,y = -log10(pvalue))) +
  geom_point(aes(color = direction),size=4) +
  scale_color_manual(name = '',
                     # color or three types
                     values = color_age,
                     # legend labels
                     label = c('9M'=paste0('9M (n=',DEnumber_pvalue_2Mvs9M[2,n],')'),
                               'NS'= paste0('NS (n=',DEnumber_pvalue_2Mvs9M[1,n],')'),
                               '2M'=paste0('2M (n=',DEnumber_pvalue_2Mvs9M[3,n],')'))) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',linewidth = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',linewidth = 0.8) +
  ggtitle(paste0("2M vs 9M ",i))+theme_classic()+
  # ylim(-2,7)+
  # xlim(-10,12)+
  geom_text_repel(data = marker,aes(x = log2FoldChange,y = -log10(pvalue),label = marker$gene_name),
                  force=20,color="black",size=4,point.padding = 0.5,hjust = 0.5,
                  arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                  segment.color="black",segment.size=0.2,nudge_y=1)+
  geom_point(data = marker,aes(x = log2FoldChange,y = -log10(pvalue)),
             color = "Black",shape = 21,size=4)+
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
write.table(DE,paste0("DE_result/DEgenename2_9Mvs2M_",i,".txt"),sep = "\t",quote = F,row.names = T,col.names = T)
ggsave(P1, filename =paste0("r_plot/002_DEpvalue_9Mvs2M_",i,".pdf"),width = 4.1, height = 4.6, useDingbats=FALSE)



A <- table(res$direction)
# 强制包含所有类别（没有的自动补0）
 A <- A[c("NS", "9M","2M")]
 A[is.na(A)] <- 0
# 转成列写入
DEnumber_pvalue_2Mvs9M[, n] <- A

colnames(DEnumber_pvalue_2Mvs9M)[n] <- i
 
  ##direction_5
  A <- table(res$direction_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue_52Mvs9M[, n] <- A
  colnames(DEnumber_pvalue_52Mvs9M)[n] <- i
  ##direction1
   A <- table(res$direction1)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue1_2Mvs9M[, n] <- A
  colnames(DEnumber_pvalue1_2Mvs9M)[n] <- i

  ##direction1_5
   A <- table(res$direction1_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue1_52Mvs9M[, n] <- A
  colnames(DEnumber_pvalue1_52Mvs9M)[n] <- i
  
  ##direction2
  A <- table(res$direction2)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue2_2Mvs9M[, n] <- A
  colnames(DEnumber_pvalue2_2Mvs9M)[n] <- i
  ##Direction
  A <- table(res$Direction)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj_2Mvs9M[, n] <- A
  colnames(DEnumber_padj_2Mvs9M)[n] <- i
  
  ##Direction_5
   A <- table(res$Direction_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj_52Mvs9M[, n] <- A
  colnames(DEnumber_padj_52Mvs9M)[n] <- i
  ##Direction1
  A <- table(res$Direction1)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj1_2Mvs9M[, n] <- A
  colnames(DEnumber_padj1_2Mvs9M)[n] <- i
  ##Direction1_5
  A <- table(res$Direction1_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "9M","2M")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj1_52Mvs9M[, n] <- A
  colnames(DEnumber_padj1_52Mvs9M)[n] <- i
  n=n+1
}
save(  DEnumber_pvalue_2Mvs9M ,DEnumber_pvalue_52Mvs9M ,
  DEnumber_pvalue1_2Mvs9M ,DEnumber_pvalue2_2Mvs9M ,
  DEnumber_padj_2Mvs9M ,DEnumber_padj_52Mvs9M ,
  DEnumber_padj1_2Mvs9M ,DEnumber_padj1_52Mvs9M ,
  file = "/share/home/shli24/2025OVX/result/DE_result/aging_DEnumber.RData")
save(  DEresult,DEgenenameresult,
  file = "/share/home/shli24/2025OVX/result/DE_result/aging_DEresult.RData")
#003 "DEnumber plot"---------------------------------------------------------
## 3.1 DEnumber plot---------------------------------------------------------------------
#module load arm/r/4.4.1
library(tidyverse)
library(ggplot2)
library(ggrepel)
color_age <- c(
  "NS" = "#7f7f7f",
  "2M" = "#009E73",   # young
  "9M" = "#CC79A7"    # aged（muted purple）
)
load("/share/home/shli24/2025OVX/result/DE_result/aging_DEnumber.RData")
A <- as.data.frame(t(DEnumber_pvalue_2Mvs9M))

df_long <- DEnumber_pvalue_2Mvs9M %>%
  filter(direction != "NS") %>%     # 不画 NS
  pivot_longer(
    cols = -direction,
    names_to = "Tissue",
    values_to = "Count"
  )
df_long <- df_long %>%
  mutate(
    Count = ifelse(direction == "2M", -Count, Count),
    Group = recode(direction,
                   "9M" = "9M",
                   "2M" = "2M")
  )
df_long$Label <- abs(df_long$Count)
df_long$Tissue <- factor(df_long$Tissue,levels = c("Liver","Lung","Kidney","Spleen","Heart"))
df_long$Position <- (df_long$Count)*1.2
#df_long$Position[c(3,8)] <- c(20,-20)
#df_long$Position[c(4,9)] <- c(50,-60)
p1 <- ggplot(df_long, aes(x = Tissue, y = Count, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, size = 0.8) +
  scale_fill_manual(values = color_age) +
  geom_text(aes(x = Tissue,y = Position, label = Label), color = "black",fontface  = "bold") +
  theme_bw() +
  labs(
    y = "Number of DEGs",
    x = ""
  ) +
  #theme_update(text = element_text(family = "Arial")) +
  theme(
   
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    # ⭐ 关键：统一边框粗细
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # ❌ 删除 axis.line
    axis.line = element_blank(),
    axis.text.y = element_text(size = 12, color = "black",face = "bold"),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    legend.text = element_text( color = "black",face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
ggsave(p1, filename ="/share/home/shli24/2025OVX/result/r_plot/002_DEnumber_2Mvs9M_.pdf",width = 4.1, height = 4.6)

# 计算每个组织的总 DEG 数（用于标注）
total_counts <- df_long %>%
  group_by(Tissue) %>%
  summarise(Total = sum(Count))

# 设置组织顺序
df_long$Tissue <- factor(df_long$Tissue, levels = c("Liver", "Lung", "Spleen", "Kidney", "Heart"))
total_counts$Tissue <- factor(total_counts$Tissue, levels = c("Liver", "Lung", "Spleen", "Kidney", "Heart"))

# 绘制堆叠柱状图
p1 <- ggplot(df_long, aes(x = Tissue, y = Count, fill = direction)) +
  geom_col(position = "stack", width = 0.7) +          # 堆叠柱状图
  geom_hline(yintercept = 0, linewidth = 0.8) +
  scale_fill_manual(values = color_age) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +   # 👈 y轴从0开始，顶部留5%空间
  # 添加总数标签（在柱顶上方）
  geom_text(data = total_counts,
            aes(x = Tissue, y = Total + 50, label = Total),   # +5 使标签略微高于柱顶
            color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  
  theme_bw() +
  labs(
    y = "Number of DEGs",
    x = "",
    fill = "Direction"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    legend.text = element_text(color = "black", face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )

# 保存图形
ggsave(p1, filename = "/share/home/shli24/2025OVX/result/r_plot/002_DEnumber_2Mvs9M_stacked.pdf",
       width = 4.2, height = 4.3, device = cairo_pdf)   # 推荐使用 Cairo PDF 避免字体问题
## 3.2 DEnumber5 plot---------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggrepel)
color_age <- c(
  "2M"   = "#009E73",
  "9M"    = "#D55E00",
  "NS" = "#7f7f7f",
  "9M+E2" = "#0072B2"  # 深蓝（干预 / 恢复）
)
load("/share/home/shli24/2025OVX/result/DE_result/DEnumber.RData")
A <- as.data.frame(t(DEnumber_pvalue_52Mvs9M))

df_long <- DEnumber_pvalue_52Mvs9M %>%
  filter(direction != "NS") %>%     # 不画 NS
  pivot_longer(
    cols = -direction,
    names_to = "Tissue",
    values_to = "Count"
  )
df_long <- df_long %>%
  mutate(
    Count = ifelse(direction == "2M", -Count, Count),
    Group = recode(direction,
                   "9M" = "9M",
                   "2M" = "2M")
  )
df_long$Label <- abs(df_long$Count)
df_long$Tissue <- factor(df_long$Tissue,levels = c("Liver","Lung","Spleen","Kidney","Heart"))
df_long$Position <- (df_long$Count)*1.2
#df_long$Position[c(3,8)] <- c(20,-20)
#df_long$Position[c(4,9)] <- c(50,-60)
p1 <- ggplot(df_long, aes(x = Tissue, y = Count, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, size = 0.8) +
  scale_fill_manual(values = color_age) +
  geom_text(aes(x = Tissue,y = Position, label = Label), color = "black",fontface  = "bold") +
  theme_bw() +
  labs(
    y = "Number of DEGs",
    x = ""
  ) +
  #theme_update(text = element_text(family = "Arial")) +
  theme(
   
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    # ⭐ 关键：统一边框粗细
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # ❌ 删除 axis.line
    axis.line = element_blank(),
    axis.text.y = element_text(size = 12, color = "black",face = "bold"),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    legend.text = element_text( color = "black",face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
ggsave(p1, filename ="/share/home/shli24/2025OVX/result/r_plot/002_DEnumber_52Mvs9M_.pdf",width = 4.1, height = 4.6)

# 计算每个组织的总 DEG 数（用于标注）
total_counts <- df_long %>%
  group_by(Tissue) %>%
  summarise(Total = sum(Count))

# 设置组织顺序
df_long$Tissue <- factor(df_long$Tissue, levels = c("Liver", "Lung", "Spleen", "Kidney", "Heart"))
total_counts$Tissue <- factor(total_counts$Tissue, levels = c("Liver", "Lung", "Spleen", "Kidney", "Heart"))

# 绘制堆叠柱状图
p1 <- ggplot(df_long, aes(x = Tissue, y = Count, fill = direction)) +
  geom_col(position = "stack", width = 0.7) +          # 堆叠柱状图
  geom_hline(yintercept = 0, linewidth = 0.8) +
  scale_fill_manual(values = color_age) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +   # 👈 y轴从0开始，顶部留5%空间
  # 添加总数标签（在柱顶上方）
  geom_text(data = total_counts,
            aes(x = Tissue, y = Total + 50, label = Total),   # +5 使标签略微高于柱顶
            color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  
  theme_bw() +
  labs(
    y = "Number of DEGs",
    x = "",
    fill = "Direction"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    legend.text = element_text(color = "black", face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )

# 保存图形
ggsave(p1, filename = "/share/home/shli24/2025OVX/result/r_plot/002_DEnumber_52Mvs9M_stacked.pdf",
       width = 4.2, height = 4.3, device = cairo_pdf)   # 推荐使用 Cairo PDF 避免字体问题


