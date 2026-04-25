setwd("/share/home/shli24/2025OVX/result/")
load("agingRNA.RData")
organ_colors <- c(
  Heart = "#E64B35",
  Liver = "#4DBBD5",
  Kidney = "#00A087",
  Lung = "#3C5488",
  Spleen = "#F39B7F"
)
library(matrixStats)
library(ggplot2)
library(Rtsne)
library(FactoMineR)
library(factoextra)
# 转换为matrix
tpm <- as.matrix(tpm)

# log2转换（非常重要）
tpm_log <- log2(tpm + 1)
var_genes <- order(rowVars(tpm_log), decreasing = TRUE)[1:2000]
expr_var <- tpm_log[var_genes, ]
expr_t <- t(expr_var)
meta <- info
pca_res <- prcomp(expr_t, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Group = meta$Group,
  Tissue = meta$Tissue
)
p1 <- ggplot(pca_df, aes(PC1, PC2, color=Tissue, shape=Group)) +
  geom_point(size=4, alpha=0.9) +
  #stat_ellipse(level=0.95, linewidth=0.8) +
  theme_classic(base_size = 14, base_family = "sans") +
  scale_fill_manual(values = organ_colors) +
  scale_color_manual(values = organ_colors) +
  labs(title="aging model",
       x=paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
       y=paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")) +
  theme(
    legend.position = "right",
    #legend.direction = "horizontal",
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "#ffffff"),
    # 文字样式（你已有的）
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10,  face = "bold", color = "black"),
    legend.text = element_text(size = 10, face = "bold", color = "black"),
    legend.key.size = unit(0.6, "cm"),legend.spacing.y = unit(0.2, "cm"),
    # ========== 新增：加粗坐标轴和刻度线 ==========
    axis.line = element_line(color = "black", linewidth = 1.1),      # 坐标轴线加粗
    axis.ticks = element_line(color = "black", linewidth = 1.1), # 刻度线加粗
    axis.ticks.length = unit(1.2, "mm")                               # 刻度线长度（可选）
  )

p1
ggsave(plot = p1,width = 4.1,height = 3.5,device = cairo_pdf,
  file = "/share/home/shli24/2025OVX/result/r_plot/000_aging_PCA.pdf"
)

set.seed(123)
tsne_res <- Rtsne(expr_t,
                  dims = 2,
                  perplexity = 3,
                  max_iter = 1000,
                  verbose = TRUE)

tsne_df <- data.frame(
  tSNE1 = tsne_res$Y[,1],
  tSNE2 = tsne_res$Y[,2],
  Group = meta$Group,
  Tissue = meta$Tissue
)
p2 <- ggplot(tsne_df, aes(tSNE1, tSNE2, color=Tissue, shape=Group)) +
  geom_point(size=4, alpha=0.9) +
  theme_classic(base_size = 14, base_family = "sans") +
  scale_fill_manual(values = organ_colors) +
  scale_color_manual(values = organ_colors) +
  labs(title="aging model") +
  theme(
    legend.position = "right",
    #legend.direction = "horizontal",
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "#ffffff"),
    # 文字样式（你已有的）
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10,  face = "bold", color = "black"),
    legend.text = element_text(size = 10, face = "bold", color = "black"),
    legend.key.size = unit(0.6, "cm"),legend.spacing.y = unit(0.2, "cm"),
    # ========== 新增：加粗坐标轴和刻度线 ==========
    axis.line = element_line(color = "black", linewidth = 1.1),      # 坐标轴线加粗
    axis.ticks = element_line(color = "black", linewidth = 1.1), # 刻度线加粗
    axis.ticks.length = unit(1.2, "mm")                               # 刻度线长度（可选）
  )
p2
ggsave(plot = p2,width = 4.1,height = 3.5,device = cairo_pdf,
  file = "/share/home/shli24/2025OVX/result/r_plot/000_aging_tSNE.pdf"
)
