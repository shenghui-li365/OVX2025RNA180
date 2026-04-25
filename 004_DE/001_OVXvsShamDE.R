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
i = "Spleen"
A_info <- info[(info$Group != "OVX+E2")&(info$Tissue == i),]
A_count<- count[,A_info$Name]
A_tpm<- tpm[,A_info$Name]

DE_info <- A_info[A_info$Tissue ==i,]
DE_count <- A_count[,DE_info$Name]

rownames(DE_count) 
rownames(DE_info) <- DE_info$Name
colnames(DE_count) <- DE_info$Name[which(colnames(DE_count)%in% DE_info$Name)]
DE_count = round(DE_count)
DE_count = DE_count[rowSums(DE_count)>=2,]
DE_tpm <- A_tpm[,DE_info$Name]
library(ggrepel)
rm(conditions)
conditions = data.frame(conditions=factor(DE_info$Group))
rownames(conditions) = colnames( DE_count)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData =  DE_count ,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
resultsNames(dds)
res = results(dds,contrast=list(c("conditions_Sham_vs_OVX")))%>%as.data.frame()
res <- mutate(res,direction = if_else(
  pvalue >0.05,"NS",if_else(
    log2FoldChange> 0,'Sham',if_else(
      log2FoldChange < 0,'OVX',"NS"
    )
  )
))
table(res$direction)
DEnumber_pvalue_ShamvsOVX <-as.data.frame(table(res$direction))
colnames(DEnumber_pvalue_ShamvsOVX) <- c("direction",i)

res <- mutate(res,direction_5 = if_else(
  pvalue >0.05,"NS",if_else(
    log2FoldChange> 0.5,'Sham',if_else(
      log2FoldChange < -0.5,'OVX',"NS"
    )
  )
))
table(res$direction_5)
DEnumber_pvalue_5ShamvsOVX <-as.data.frame(table(res$direction_5))
colnames(DEnumber_pvalue_5ShamvsOVX) <- c("direction",i)

res <- mutate(res,direction1 = if_else(
  pvalue >0.05,"NS",if_else(
    log2FoldChange> 1,'Sham',if_else(
      log2FoldChange < -1,'OVX',"NS"
    )
  )
))
table(res$direction1)
DEnumber_pvalue1_ShamvsOVX <-as.data.frame(table(res$direction1))
colnames(DEnumber_pvalue1_ShamvsOVX) <- c("direction",i)

res <- mutate(res,direction1_5 = if_else(
  pvalue >0.05,"NS",if_else(
    log2FoldChange> 1.5,'Sham',if_else(
      log2FoldChange < -1.5,'OVX',"NS"
    )
  )
))
table(res$direction1_5)
DEnumber_pvalue1_5ShamvsOVX <-as.data.frame(table(res$direction1_5))
colnames(DEnumber_pvalue1_5ShamvsOVX) <- c("direction",i)

res <- mutate(res,direction2 = if_else(
  pvalue >0.05,"NS",if_else(
    log2FoldChange> 2,'Sham',if_else(
      log2FoldChange < -2,'OVX',"NS"
    )
  )
))
table(res$direction2)
DEnumber_pvalue2_ShamvsOVX <-as.data.frame(table(res$direction2))
colnames(DEnumber_pvalue2_ShamvsOVX) <- c("direction",i)

res <- mutate(res,Direction = if_else(
  padj >0.05,"NS",if_else(
    log2FoldChange> 0,'Sham',if_else(
      log2FoldChange < 0,'OVX',"NS"
    )
  )
))
table(res$Direction)
DEnumber_padj_ShamvsOVX <-as.data.frame(table(res$Direction))
colnames(DEnumber_padj_ShamvsOVX) <- c("Direction",i)

res <- mutate(res,Direction_5 = if_else(
  padj >0.05,"NS",if_else(
    log2FoldChange> 0.5,'Sham',if_else(
      log2FoldChange < -0.5,'OVX',"NS"
    )
  )
))
table(res$Direction_5)
DEnumber_padj_5ShamvsOVX <-as.data.frame(table(res$direction))
colnames(DEnumber_padj_5ShamvsOVX) <- c("Direction",i)

res <- mutate(res,Direction1 = if_else(
  padj >0.05,"NS",if_else(
    log2FoldChange> 1,'Sham',if_else(
      log2FoldChange < -1,'OVX',"NS"
    )
  )
))
table(res$Direction1)
DEnumber_padj1_ShamvsOVX <-as.data.frame(table(res$Direction1))
colnames(DEnumber_padj1_ShamvsOVX) <- c("Direction",i)

res <- mutate(res,Direction1_5 = if_else(
  padj >0.05,"NS",if_else(
    log2FoldChange> 1.5,'Sham',if_else(
      log2FoldChange < -1.5,'OVX',"NS"
    )
  )
))
table(res$Direction1_5)
DEnumber_padj1_5ShamvsOVX <-as.data.frame(table(res$Direction1_5))
colnames(DEnumber_padj1_5ShamvsOVX) <- c("Direction",i)

genename <- read.delim("gene_name.txt")
write.table(res,paste0("DE_result/DE_ShamvsOVX_",i,".txt"),sep = "\t",quote = F,row.names = T,col.names = T)
res$gene_id <-rownames(res)
res <-res[!is.na(res$direction),]
DE <- res[which(res$direction!="NS"),]
DE <-merge(DE,genename,by="gene_id", all.x = FALSE)

marker <- rbind((DE %>%
                   arrange(desc(log2FoldChange)) %>%
                   slice_head(n = 10)),
                DE %>%
                  arrange(log2FoldChange) %>%
                  slice_head(n = 10))
marker <- marker[!is.na(marker$gene_name),]
n=2
P1 <- ggplot(res,aes(x = log2FoldChange,y = -log10(pvalue))) +
  geom_point(aes(color = direction),size=4) +
  scale_color_manual(name = '',
                     # color or three types
                     values = color_ovx_e2,
                     # legend labels
                     label = c('OVX'=paste0('OVX (n=',DEnumber_pvalue_ShamvsOVX[2,n],')'),
                               'NS'= paste0('NS (n=',DEnumber_pvalue_ShamvsOVX[1,n],')'),
                               'Sham'=paste0('Sham (n=',DEnumber_pvalue_ShamvsOVX[3,n],')'))) +
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
write.table(DE,paste0("DE_result/DEgenename_ShamvsOVX_",i,".txt"),sep = "\t",quote = F,row.names = T,col.names = T)
ggsave(P1, filename =paste0("r_plot/002_DEpvalue_ShamvsOVX_",i,".pdf"),width = 4.1, height = 4.6, useDingbats=FALSE)
n = 2
i = "Heart"
genename <- read.delim("gene_name.txt")
for (i in unique(info$Tissue)) {
  A_info <- info[(info$Group != "OVX+E2")&(info$Tissue == i),]
  A_count<- count[,A_info$Name]
  A_tpm<- tpm[,A_info$Name]
  
  DE_info <- A_info[A_info$Tissue ==i,]
  DE_count <- A_count[,DE_info$Name]
  colnames(DE_info)
  rownames(DE_count) 
  rownames(DE_info) <- DE_info$Name
  colnames(DE_count) <- DE_info$Name[which(colnames(DE_count)%in% DE_info$Name)]
  DE_count = round(DE_count)
  DE_count = DE_count[rowSums(DE_count)>=2,]
  DE_tpm <- A_tpm[,DE_info$Name]
  rm(conditions)
  conditions = data.frame(conditions=factor(DE_info$Group))
  rownames(conditions) = colnames( DE_count)
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData =  DE_count ,
    colData = conditions,
    design = ~ conditions)
  dds = DESeq(ddsFullCountTable)
  resultsNames(dds)
  res = results(dds,contrast=list(c("conditions_Sham_vs_OVX")))%>%as.data.frame()
  ##direction
  res <- mutate(res,direction = if_else(
    pvalue >0.05,"NS",if_else(
      log2FoldChange> 0,'Sham',if_else(
        log2FoldChange < 0,'OVX',"NS"
      )
    )
  ))
  A <- table(res$direction)
# 强制包含所有类别（没有的自动补0）
 A <- A[c("NS", "OVX","Sham")]
 A[is.na(A)] <- 0
# 转成列写入
DEnumber_pvalue_ShamvsOVX[, n] <- A
# 设置列名
colnames(DEnumber_pvalue_ShamvsOVX)[n] <- i
 
  ##direction_5
  res <- mutate(res,direction_5 = if_else(
    pvalue >0.05,"NS",if_else(
      log2FoldChange> 0.5,'Sham',if_else(
        log2FoldChange < -0.5,'OVX',"NS"
      )
    )
  ))
  A <- table(res$direction_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue_5ShamvsOVX[, n] <- A
  colnames(DEnumber_pvalue_5ShamvsOVX)[n] <- i
  ##direction1
  res <- mutate(res,direction1 = if_else(
    pvalue >0.05,"NS",if_else(
      log2FoldChange> 1,'Sham',if_else(
        log2FoldChange < -1,'OVX',"NS"
      )
    )
  ))
  A <- table(res$direction1)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue1_ShamvsOVX[, n] <- A
  colnames(DEnumber_pvalue1_ShamvsOVX)[n] <- i
  ##direction1_5
  res <- mutate(res,direction1_5 = if_else(
    pvalue >0.05,"NS",if_else(
      log2FoldChange> 1.5,'Sham',if_else(
        log2FoldChange < -1.5,'OVX',"NS"
      )
    )
  ))
  A <- table(res$direction1_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue1_5ShamvsOVX[, n] <- A
  colnames(DEnumber_pvalue1_5ShamvsOVX)[n] <- i
  
  ##direction2
  res <- mutate(res,direction1_5 = if_else(
    pvalue >0.05,"NS",if_else(
      log2FoldChange> 2,'Sham',if_else(
        log2FoldChange < -2,'OVX',"NS"
      )
    )
  ))
  A <- table(res$direction2)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_pvalue2_ShamvsOVX[, n] <- A
  colnames(DEnumber_pvalue2_ShamvsOVX)[n] <- i
  ##Direction
  res <- mutate(res,Direction = if_else(
    padj >0.05,"NS",if_else(
      log2FoldChange> 0,'Sham',if_else(
        log2FoldChange < 0,'OVX',"NS"
      )
    )
  ))
  A <- table(res$Direction)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj_ShamvsOVX[, n] <- A
  colnames(DEnumber_padj_ShamvsOVX)[n] <- i
  
  ##Direction_5
  res <- mutate(res,Direction_5 = if_else(
    padj >0.05,"NS",if_else(
      log2FoldChange> 0.5,'Sham',if_else(
        log2FoldChange < -0.5,'OVX',"NS"
      )
    )
  ))
  table(res$Direction_5)
  A <- table(res$Direction_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj_5ShamvsOVX[, n] <- A
  colnames(DEnumber_padj_5ShamvsOVX)[n] <- i
  ##Direction1
  res <- mutate(res,Direction1 = if_else(
    padj >0.05,"NS",if_else(
      log2FoldChange> 1,'Sham',if_else(
        log2FoldChange < -1,'OVX',"NS"
      )
    )
  ))
  A <- table(res$Direction1)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj1_ShamvsOVX[, n] <- A
  colnames(DEnumber_padj1_ShamvsOVX)[n] <- i
  ##Direction1_5
  res <- mutate(res,Direction1_5 = if_else(
    padj >0.05,"NS",if_else(
      log2FoldChange> 1.5,'Sham',if_else(
        log2FoldChange < -1.5,'OVX',"NS"
      )
    )
  ))
  table(res$Direction_5)
  A <- table(res$Direction1_5)
# 强制包含所有类别（没有的自动补0）
  A <- A[c("NS", "OVX","Sham")]
  A[is.na(A)] <- 0
# 转成列写入
  DEnumber_padj1_5ShamvsOVX[, n] <- A
  colnames(DEnumber_padj1_5ShamvsOVX)[n] <- i
  write.table(res,paste0("DE_result/DE_ShamvsOVX_",i,".txt"),sep = "\t",quote = F,row.names = T,col.names = T)
  res$gene_id <-rownames(res)
  res <-res[!is.na(res$direction),]
  DE <- res[which(res$direction!="NS"),]
  DE <-merge(DE,genename,by="gene_id", all.x = FALSE)
  
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
                       values = color_ovx_e2,
                       # legend labels
                       label = c('OVX'=paste0('OVX (n=',DEnumber_pvalue_ShamvsOVX[2,n],')'),
                                 'NS'= paste0('NS (n=',DEnumber_pvalue_ShamvsOVX[1,n],')'),
                                 'Sham'=paste0('Sham (n=',DEnumber_pvalue_ShamvsOVX[3,n],')'))) +
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
  write.table(DE,paste0("DE_result/DEgenename_ShamvsOVX_",i,".txt"),sep = "\t",quote = F,row.names = T,col.names = T)
  ggsave(P1, filename =paste0("r_plot/002_DEpvalue_ShamvsOVX_",i,".pdf"),width = 4.1, height = 4.6, useDingbats=FALSE)
  n=n+1
}
save(DEnumber_pvalue_ShamvsOVX,DEnumber_pvalue_5ShamvsOVX,DEnumber_pvalue1_ShamvsOVX,
     DEnumber_pvalue1_5ShamvsOVX,DEnumber_pvalue2_ShamvsOVX,
     DEnumber_padj_ShamvsOVX,DEnumber_padj_5ShamvsOVX,DEnumber_padj1_ShamvsOVX,
     DEnumber_padj1_5ShamvsOVX,
     file = "DE_result/DEnumber_ShamvsOVX.RData"
)
DEnumber_pvalue_ShamvsOVX[,-3]
DEnumber_pvalue_5ShamvsOVX[,-3]
DEnumber_pvalue1_ShamvsOVX[,-3]
DEnumber_pvalue1_5ShamvsOVX[,-3]
DEnumber_pvalue2_ShamvsOVX[,-3]
DEnumber_padj_ShamvsOVX[,-3]
DEnumber_padj_5ShamvsOVX[,-3]
DEnumber_padj1_ShamvsOVX[,-3]
DEnumber_padj1_5ShamvsOVX[,-3]
