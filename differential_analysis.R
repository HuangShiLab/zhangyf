rm(list = ls())

library(crossRanger)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(biomformat)
library(ggrepel)
library(DESeq2)

table <- read.table("../06_GeneExpression/all.fpkm.xls", sep = "\t", header = T, row.names = 1)
meta <- read.table("meta.txt", sep = "\t", header = T, row.names = 1)

# fpkm_columns <- grep("_FPKM$", colnames(table), value = TRUE)
# fpkm <- table[, fpkm_columns]
# colnames(fpkm) <- gsub("_FPKM", "", colnames(fpkm))
fpkm_columns <- grep("_Count$", colnames(table), value = TRUE)
fpkm <- table[, fpkm_columns]
colnames(fpkm) <- gsub("_Count", "", colnames(fpkm))
fpkm <- subset(fpkm, rowSums(fpkm) != 0)
fpkm <- t(fpkm)
meta <- subset(meta, rownames(meta) %in% rownames(fpkm))
identical(rownames(fpkm), rownames(meta))

df <- melt(fpkm)
colnames(df) <- c("Samples", "Gene", "FPKM")
A_P_df <- subset(df, Samples %in% c("A1", "A2", "A3", "P1", "P2", "P3"))
A_P_df$Samples <- factor(A_P_df$Samples, levels = c("A1", "A2", "A3", "P1", "P2", "P3"))

Volcano_plot <- function(data) {
  p <- ggplot(data, aes(x = mean_logfc, y = `-log10 (qvalue)`, colour=Enr2)) + #x、y轴取值限制，颜色根据"Sig"
    geom_point(alpha=0.65, size=2) +  #点的透明度、大小
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
    xlim(c(-10, 10)) +  #调整点的颜色和x轴的取值范围
    geom_vline(xintercept= c(-1, 1), lty=4, col="black", lwd=0.8) + #添加x轴辅助线
    geom_hline(yintercept = -log10(0.05), lty=4, col="black", lwd=0.8) +  #添加y轴辅助线
    labs(x="log2FC", y="-log10 (qvalue)") +  #x、y轴标签
    theme_bw() + # 主题，help(theme)查找其他个性化设置
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right", 
          legend.title = element_blank())
  
  return (p)
}

differential_analysis2 <- function(data, meta, group1, group2) {
  x <- subset(data, meta[rownames(data), "group"] %in% c(group1, group2))
  x <- t(x)
  x <- as.data.frame(apply(x, c(1,2), as.integer))
  y <- subset(meta, group %in% c(group1, group2))
  y <- factor(y$group, levels = c(group1, group2))
  y <- data.frame(group = y)
  rownames(y) <- colnames(x)
  dds <- DESeqDataSetFromMatrix(countData = x, colData = y, design = ~group)
  dds1 <- DESeq(dds)
  res <- results(dds1, contrast = c("group", group2, group1))
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  #log2FC≥1 & padj<0.01 标识 Up，代表显著上调的基因
  #log2FC≤-1 & padj<0.01 标识 Down，代表显著下调的基因
  #其余标识 NoSignificant，代表非差异的基因
  res1[which(res1$log2FoldChange > 1 & res1$padj < 0.05),'sig'] <- 'Up'
  res1[which(res1$log2FoldChange < -1 & res1$padj < 0.05),'sig'] <- 'Down'
  res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'NoSignificant'
  
  write.table(res1, paste0(group2, "_VS_", group1, ".DESeq2.txt"), col.names = NA, sep = '\t', quote = FALSE)
  
  p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(size = 1) +  #绘制散点图
    scale_color_manual(values = c('darkred', 'gray', 'darkgreen'), limits = c('Up', 'NoSignificant', 'Down')) +  #自定义点的颜色
    labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = paste0(group2, "_VS_", group1), color = '') +  #坐标轴标题
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
          # panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
    geom_hline(yintercept = 2, lty = 3, color = 'black') +
    xlim(-10, 10) + ylim(0, 40)  #定义刻度边界
  
  ggsave(paste0(group2, "_VS_", group1, ".Volcano_plot.pdf"), p, width = 5, height = 4)
  
  return (res1)
}

P_A_test_output <- differential_analysis2(fpkm, meta, "P", "A")

result <- matrix(, nrow = 2, ncol = 1, dimnames = list(c("Down", "Up"), c("A_VS_P")))
result <- as.data.frame(result)
P_A_test_output <- table(P_A_test_output$sig)
result["Down", "A_VS_P"] <- P_A_test_output["Down"]
result["Up", "A_VS_P"] <- P_A_test_output["Up"]

result <- as.matrix(result)
pdf("Bar graph of genes significantly up- or down-regulation between groups.pdf", width=4, height=5)
p <- barplot(result,
        col = c("darkgreen", "darkred"),
        xlab = "Group",
        ylab = "Number",
        beside = TRUE)
text(p, result + 20, format(result), xpd = TRUE, col = "black")
legend("topright",
       legend = c("Down", "Up"),
       fill = c("darkgreen", "darkred"),)
dev.off()



