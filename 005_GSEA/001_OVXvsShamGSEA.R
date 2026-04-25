 #module load arm/r/4.4.1
 #R
 library(clusterProfiler)
 library(readxl)
 library(tidyverse)
 library(ComplexHeatmap)
setwd("/share/home/shli24/2025OVX/result/DE_result")
tissue <- c("Heart","Liver","Spleen","Lung","Kidney")
i = "Heart"
aging_hallmarkers_table<-read_xlsx("/share/home/shli24/2025OVX/data/Aging_hallmarkers_GSEA.xlsx")
All_GSEA_results <- list()

 genename <- read.delim("/share/home/shli24/2025OVX/result/gene_name.txt")
for (i in tissue) {
  
  DE_results <- read.delim(paste0("DE2_ShamvsOVX_",i,".txt"), row.names=1)
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
  
  set.seed(123)
  GSEA_results <- GSEA(geneList, TERM2GENE=aging_hallmarkers_table, verbose=FALSE, 
                       minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  results_table<-GSEA_results@result
  results_table$Condition <- i
  All_GSEA_results[[i]] <- results_table
  
  
  p <-ridgeplot(GSEA_results,
                showCategory = 20,
                fill = "pvalue", #濉厖鑹?"pvalue", "p.adjust", "qvalue" 
                core_enrichment = TRUE,#鏄惁鍙娇鐢?core_enriched gene
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


pdf("/share/home/shli24/2025OVX/result/r_plot/003_GSEA_stat_ShamvsOVX_all.pdf",width = 7,height = 5)

# 濡傛灉姣忎釜 list 鍏冪礌閮芥槸鏁版嵁妗?鐭╅樀
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
  library(circlize)
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

Heatmap(
  mat,
  name = " ",
  col = col_fun,
  cluster_columns = T,   # 鍥哄畾椤哄簭
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
    fontfamily = "sans",
    fontface   = "bold",
    fontsize   = 12
  ),
  
  row_names_gp = gpar(
    fontfamily = "sans",
    fontface   = "bold",
    fontsize   = 12
  ),
  
  heatmap_legend_param = list(
    labels_gp = gpar(fontfamily = "sans", fontface = "bold"),
    title_gp  = gpar(fontfamily = "sans", fontface = "bold")
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


