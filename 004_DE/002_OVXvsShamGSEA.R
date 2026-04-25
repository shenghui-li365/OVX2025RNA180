
#module load arm/r/4.4.1
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)   # 人类，小鼠请换成 org.Mm.eg.db
library(openxlsx)
library(dplyr)
library(readxl)
library(circlize)
library(ComplexHeatmap)
color_ovx_e2 <- c(
  "Sham"   = "#009E73",
  "OVX"    = "#D55E00",
  "NS" = "#7f7f7f",
  "OVX+E2" = "#0072B2"  # 深蓝（干预 / 恢复）
)
load("/share/home/shli24/2025OVX/result/DE_result/DEresult.RData")
tissue <- c("Heart","Liver","Spleen","Lung","Kidney")
genename <- read.delim("/share/home/shli24/2025OVX/result/gene_name.txt")


aging_hallmarkers_table <- read_xlsx("/share/home/shli24/2025OVX/data/Aging_hallmarkers_GSEA.xlsx")
setwd("/share/home/shli24/2025OVX/result")
i = "Liver"
All_GSEA_results <- list()
pdf("r_plot/002_ShamvsOVX_GSEAall.pdf",width = 7,height = 5)
for(i in names(DEresult)){
  DE_results <- read.delim(paste0("DE_result/DE2_ShamvsOVX_",i,".txt"), row.names=1)
  if( "gene_id" %in% colnames(DE_results) ==0 ){DE_results$gene_id = rownames(DE_results)}
  DE_results <-merge(DE_results,genename,by="gene_id", all.x = FALSE)
  DE_results$Symbol_toupp=toupper(DE_results$gene_name)
  DE_results<-DE_results[duplicated(DE_results$Symbol_toupp)==FALSE,]
  nrow(DE_results)
  #geneList <- DE_results$log2FoldChange
  geneList <- DE_results$stat
  names(geneList) <- DE_results$Symbol_toupp
  geneList <- sort(geneList, decreasing = T)
  geneList <- geneList[geneList != 0]
  
  set.seed(10000)
  GSEA_results <- GSEA(geneList, TERM2GENE=aging_hallmarkers_table, verbose=FALSE, 
                       minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  results_table<-GSEA_results@result
  results_table$Condition <- i
  All_GSEA_results[[i]] <- results_table
  p <-ridgeplot(GSEA_results,
            showCategory = 20,
            fill = "pvalue", #填充色 "pvalue", "p.adjust", "qvalue" 
            core_enrichment = TRUE,#是否只使用 core_enriched gene
            label_format = 30,
            orderBy = "NES",
            decreasing = FALSE
  )+
    theme_bw()+
    ggtitle(i)+
    theme(plot.title = element_text(size = 25),
          axis.text = element_text(color = 'black',size=20),
          axis.title.x  = element_blank(),
          axis.title = element_text(color = 'black'),
          axis.title.y   = element_text(color = 'black',size=20),
          legend.text = element_text(size=15),
          legend.title = element_text(size=20)
    )
  print(p)
  rm(results_table)
  rm(GSEA_results)
}


data <- do.call(rbind, All_GSEA_results)
unique(data$Description)
data$Description[which(data$Description == "Inflammatory response")] <- "Inflammatory"

mat <- data %>%
  dplyr::select(Description, Condition, NES) %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = NES
  ) %>%
  tibble::column_to_rownames("Description") %>%
  as.matrix()
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

Heatmap(
  mat,
  name = " ",
  col = col_fun,
  cluster_columns = T,   # 固定顺序
  cluster_rows = FALSE,
  column_names_rot = 45,
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.circle(
      x = x,
      y = y,
      r = min(w, h) * 0.42,
      gp = gpar(fill = fill, col = NA)
    )
  }
)
mat1 <- mat[c(2,8:10),]
Heatmap(
  mat1,
  name = " ",
  col = col_fun,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  column_names_rot = 45,
  
  column_names_gp = gpar(
    fontfamily = "Arial",
    fontface   = "bold",
    fontsize   = 12
  ),
  
  row_names_gp = gpar(
    fontfamily = "Arial",
    fontface   = "bold",
    fontsize   = 12
  ),
  
  heatmap_legend_param = list(
    labels_gp = gpar(fontfamily = "Arial", fontface = "bold"),
    title_gp  = gpar(fontfamily = "Arial", fontface = "bold")
  ),
  
  rect_gp = gpar(type = "none"),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.circle(
      x = x,
      y = y,
      r = min(w, h) * 0.42,
      gp = gpar(fill = fill, col = NA)
    )
  }
)
dev.off()
save(  DEresult,DEgenenameresult,
  file = "/share/home/shli24/2025OVX/result/DE_result/DEresult.RData")