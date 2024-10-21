library(reshape2)
library(ggplot2)
library(RColorBrewer)
table <-reshape2table <- read.table("../06_GeneExpression/all.fpkm_by_company.xls", sep = "\t", header = T, row.names = 1)
meta <- read.table("meta.txt", sep = "\t", header = T, row.names = 1)

fpkm_columns <- grep("_FPKM$", names(table), value = TRUE)
fpkm <- table[, fpkm_columns]
colnames(fpkm) <- gsub("_FPKM", "", colnames(fpkm))
fpkm <- t(fpkm)
meta <- subset(meta, rownames(meta) %in% rownames(fpkm))
identical(rownames(fpkm), rownames(meta))

df <- melt(fpkm)
colnames(df) <- c("Samples", "Gene", "FPKM")
matched_index <- match(df$Samples, rownames(meta))
df$group <- meta[matched_index, "group"]
boxplot <- ggplot(df, aes(x = Samples, y = log10(FPKM))) +
  geom_boxplot(outliers = T, outlier.shape = NA) +
  scale_colour_manual(values = c("darkred", "darkgreen")) +
  geom_jitter(aes(color = group), shape = 1, size = 0.05, width = 0.15) +
  xlab("Samples") + 
  ylab("log10(FPKM)") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
boxplot
ggsave("FPKM_boxplot.pdf", plot = boxplot, width = 4, height = 3)
