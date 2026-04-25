#module load arm/r/4.4.1
library(clusterProfiler)
library(org.Mm.eg.db)
setwd("/share/home/shli24/2025OVX/result/")
## 0. 安装并加载所需包 ---------------------------------------------
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler","org.Hs.eg.db","openxlsx"))
library(clusterProfiler)
library(org.Mm.eg.db)   # 人类，小鼠请换成 org.Mm.eg.db
library(openxlsx)
library(dplyr)
GO_top10 <- list()
GO <- list()  # 需要添加GO列表初始化
KEGG_result <- list()  # 需要添加KEGG_result列表初始化
KEGG_top10 <- list()  # 需要添加KEGG_top10列表初始化
file = list.files("/share/home/shli24/2025OVX/result/DE_result",pattern = "DEgenename2_")
group <- c("direction","direction_5","direction1","direction1_5","direction2","Direction","Direction_5","Direction1","Direction1_5")
m = "direction"
i = "DEgenename2_ShamvsOVX_Heart.txt"
for (i in file){
    a <- gsub("DEgenename2_","",gsub(".txt","",i))
    DE_results <- read.delim(paste0("DE_result/",i))
    for (m in group){
        DE_results <- DE_results[DE_results[,m]!="NS",]
        DE_results <- DE_results[!is.na(DE_results$gene_name), ]
        all_genes <- DE_results$gene_name
        ## 1. 统一把 Symbol 转成 ENTREZID -----------------------------------
    all_genes <- DE_results$gene_name
    name_ID <- bitr(all_genes,
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Mm.eg.db)
    
    # GO分析
    ego_BP <- enrichGO(gene = name_ID$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = 'ENTREZID',
                       ont = "BP",
                       pAdjustMethod = "BH",
                       minGSSize = 5,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
    if (!is.null(ego_BP) && nrow(ego_BP) > 0) {
    ego_BP <- as.data.frame(ego_BP) %>%
        arrange(pvalue) %>%
        mutate(
            BgTotal = as.numeric(sub(".*/", "", BgRatio)),
            RichFactor = Count / BgTotal
        )
    f <- paste(a, m, sep = "_")
    ego_BP$cell <-  f
    GO[[f]] <- ego_BP
    GO_top10[[f]] <- head(ego_BP, 10)
    write.table(ego_BP, paste0("GO_result/", f, "_all", ".csv"), sep = ",", quote = F, row.names = F, col.names = T)
    rm(ego_BP)
            }else {
            # 如果KEGG结果为空，创建空数据框并跳过后续操作
            warning(paste("No GO results for",f))

        }
        # KEGG分析
    KEGG <- enrichKEGG(gene = name_ID$ENTREZID,
                       organism = 'mmu',
                       keyType = 'kegg',
                       pvalueCutoff = 0.05,
                       pAdjustMethod = 'hochberg')
    
    # 检查KEGG结果是否为空
    if (!is.null(KEGG) && nrow(KEGG) > 0) {
        KEGG <- as.data.frame(KEGG) %>%
            arrange(pvalue) %>%
            mutate(
                BgTotal = as.numeric(sub(".*/", "", BgRatio)),
                RichFactor = Count / BgTotal
            )
        KEGG$cell <- f
        KEGG_result[[f]] <- KEGG
        KEGG_top10[[f]] <- head(KEGG, 10)
        write.table(KEGG, paste0("KEGG_result/", f, "_all", ".csv"), sep = ",", quote = F, row.names = F, col.names = T)
    } else {
        # 如果KEGG结果为空，创建空数据框并跳过后续操作
        warning(paste("No KEGG results for", f))
    }
    rm(KEGG)
    for (n in unique(DE_results[,m])){
        gene <- DE_results[which(DE_results[,m] == n),]
        name_ID <- bitr(gene$gene_name,  # 修改这里：使用gene而不是rownames(gene)
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Mm.eg.db)
                       # GO分析
        ego_BP <- enrichGO(gene = name_ID$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           keyType = 'ENTREZID',
                           ont = "BP",
                           pAdjustMethod = "BH",
                           minGSSize = 5,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
        if (!is.null(ego_BP) && nrow(ego_BP) > 0) {
        ego_BP <- as.data.frame(ego_BP) %>%
            arrange(pvalue) %>%
            mutate(
                BgTotal = as.numeric(sub(".*/", "", BgRatio)),
                RichFactor = Count / BgTotal
            )
        g <- paste(a, m,n,sep = "_")
        ego_BP$cell <- g
        GO[[g]] <- ego_BP
        GO_top10[[g]] <- head(ego_BP, 10)
        write.table(ego_BP, paste0("GO_result/", g,".csv"), sep = ",", quote = F, row.names = F, col.names = T)
        rm(ego_BP)
        }else {
            # 如果KEGG结果为空，创建空数据框并跳过后续操作
            warning(paste("No GO results for", g))

        }
        # KEGG分析
        KEGG <- enrichKEGG(gene = name_ID$ENTREZID,
                           organism = 'mmu',
                           keyType = 'kegg',
                           pvalueCutoff = 0.05,
                           pAdjustMethod = 'hochberg')
        
        # 检查KEGG结果是否为空
        if (!is.null(KEGG) && nrow(KEGG) > 0) {
            KEGG <- as.data.frame(KEGG) %>%
                arrange(pvalue) %>%
                mutate(
                    BgTotal = as.numeric(sub(".*/", "", BgRatio)),
                    RichFactor = Count / BgTotal
                )
            KEGG$cell <-  g
            KEGG_result[[ g]] <- KEGG
            KEGG_top10[[ g]] <- head(KEGG, 10)
            write.table(KEGG, paste0("KEGG_result/", g, ".csv"), sep = ",", quote = F, row.names = F, col.names = T)
        } else {
            # 如果KEGG结果为空，创建空数据框并跳过后续操作
            warning(paste("No KEGG results for",g))

        }
        rm(KEGG)
    }
    }
    }

save(KEGG_result,GO,KEGG_top10,GO_top10,file = "All_TOP10_GO_KEGG.Rdata")
#load("All_TOP10_GO_KEGG.Rdata")
library(ggplot2)
pdf("GO_result/All_TOP10_GO.pdf",width = 8 ,height=6.1)
for(i in names(GO_top10)){
    bb <- GO_top10[[i]]
    bb$description <- sub(" - Mus musculus \\(house mouse\\)", "", bb$Description)
    bb$description <- factor(bb$description,levels = bb$description)
    p <- ggplot(bb, aes(x=-log(pvalue),y=description,fill = RichFactor)
) +
  geom_bar(position="dodge",stat="identity",show.legend = T) +
  scale_fill_gradient2(
  low  = "#2166ac",   # 深蓝
  mid  = "white",     # 中间白色
  high = "#b2182b",   # 深红
  midpoint = 0        # 对称中心
)+
  theme_classic()+
  ggtitle(i)+
  ylab("")+xlab("-(log Pvalue)")+
  # scale_y_discrete(limits=unique(BP_up$condition))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL))+	
  #theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 12,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))
        print(p)
        rm(bb)
}
dev.off()

pdf("KEGG_result/All_TOP10_KEGG.pdf",width = 8 ,height=6.1)
for(i in names(KEGG_top10)){
    bb <- KEGG_top10[[i]]
    bb$description <- sub(" - Mus musculus \\(house mouse\\)", "", bb$Description)
    bb$description <- factor(bb$description,levels = bb$description)
    p <- ggplot(bb, aes(x=-log(pvalue),y=description,fill = RichFactor)
) +
  geom_bar(position="dodge",stat="identity",show.legend = T) +
  scale_fill_gradient2(
  low  = "#2166ac",   # 深蓝
  mid  = "white",     # 中间白色
  high = "#b2182b",   # 深红
  midpoint = 0        # 对称中心
)+
  theme_classic()+
  ggtitle(i)+
  ylab("")+xlab("-(log Pvalue)")+
  # scale_y_discrete(limits=unique(BP_up$condition))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL))+	
  #theme(axis.title =element_text(size = 18,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 12,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",
        legend.title = element_text(size=14,face = "bold"),
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))
        print(p)
        rm(bb)
}
dev.off()
    
    
    