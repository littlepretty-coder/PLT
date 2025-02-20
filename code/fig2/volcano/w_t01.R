library(tidyr)
library(nnet)
library(ggVolcano)
#rm(list=ls())
tax = read.delim('./data/tax.txt')
rownames(tax) = tax$Species
rownames(tax) = gsub('-','_',rownames(tax))
patterns <- c("\\[", "\\]", "\\(", "\\)", "\\.", "\\'", "\\/", "\\:", "\\=", "\\+")
rownames(tax) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(tax) )

data0 = read.delim('./fig2/new/真菌数据.txt')
data1 = read.delim('./fig2/new/细菌数据.txt')
#data = data[,-2]
data = rbind(data0,data1)
data$group = c(rep('fungi',nrow(data0)),rep('bacteria',nrow(data1)))
# w_t01= data0[,c(1,10,14,22)]
# w_t12= data0[,c(1,10,14,22)]
# w_t23= data0[,c(1,10,14,22)]
# wo_t01= data0[,c(1,10,14,22)]
# wo_t13= data0[,c(1,10,14,22)]
w_t01 = data[,c(1,10,11,17,26)]
w_t01 = merge(w_t01,tax[,6:7],by.x = 'key',by.y = 'Species')
w_t01 = w_t01[w_t01$Genus %in% names,]
mapping <- list(
  bac_1 = bac_1,
  bac_2 = bac_2,
  bac_3 = bac_3,
  bac_4 = bac_4,
  fungi_1 = fungi_1,
  fungi_2 = fungi_2,
  fungi_3 = fungi_3,
  fungi_4 = fungi_4
)
w_t01$Module <- 0
for (name in names(mapping)) {
  w_t01$Module[w_t01$Genus %in% mapping[[name]]] <- name
}
w_t01 <- w_t01 %>%
  mutate(Color = case_when(
    Module == "bac_1" ~ "#A51C36",
    Module == "bac_2" ~ "#7ABBDB",
    Module == "bac_3" ~ "#84BA42",
    Module == "bac_4" ~ "#682487",
    Module == 'fungi_1' ~ '#DBB428',
    Module == 'fungi_2' ~ '#D4562E',
    Module == 'fungi_3' ~ '#FFFFB3',
    Module == 'fungi_4' ~ '#FCCDE5'
  ))


w_t01 = w_t01[rowSums(w_t01[,2:3])>0,]
w_t01$fdr <- p.adjust(w_t01$p_wt01, method = "fdr")
#拉普拉斯变换 最小非零值的1/10 最小非零值为0.75 这里取0.1
min <- min(unlist(w_t01[,2:3][w_t01[,2:3] != 0]))
w_t01$log2fc = log2((w_t01$median_wt1+min)/(w_t01$median_wt0+min))
data <- add_regulate(w_t01, log2FC_name = "log2fc",
                     fdr_name = "fdr",log2FC = 1, fdr = 0.01)

# plot
gradual_volcano(data, x = "log2FoldChange", y = "padj",
                label = "key", label_number = 0, output = FALSE)

data$significance <- ifelse(
  data$padj < 0.05 & data$log2FoldChange > 1, "Up",  # 显著上调
  ifelse(data$padj < 0.05 & data$log2FoldChange < -1, "Down", "Not Significant")  # 显著下调或不显著
)

# 转换P值为-log10(pvalue)
data$negLogP <- -log10(data$padj)

# 绘制火山图
# 假设你的数据框包含列：log2FoldChange, negLogP, significance, name
library(ggplot2)
library(ggrepel)

# 按条件选择左上角和右上角的点
# top_left <- data[data$log2FoldChange < -1 & data$negLogP > -log10(0.05), ] %>%
#   dplyr::arrange(-negLogP) %>%  # 按 negLogP 降序排列
#   head(2)                      # 取前 5 个点
# 
# top_right <- data[data$log2FoldChange > 1 & data$negLogP > -log10(0.05), ] %>%
#   dplyr::arrange(-negLogP) %>%
#   head(2)

microbes <- c(
  "Megamonas_funiformis",
  "Mediterraneibacter_gnavus",
  "Alistipes_onderdonkii",
  "Faecalibacterium_prausnitzii",
  "Bacteroides_ovatus",
  "Klebsiella_aerogenes",
  "Bifidobacterium_breve",
  "Hungatella_hathewayi",
  "Lacticaseibacillus_paracasei",
  "Streptococcus_lactarius",
  "Enterococcus_faecalis",
  "Bifidobacterium_pseudocatenulatum",
  "Aspergillus_nomiae",
  "Neurospora_tetrasperma"
)
key16 = data[(data$key %in% microbes)&(data$padj < 0.05),]
# 将两个选择的点合并
top_points <- key16
data[data[, 12] == "Not Significant", 8] <- "grey"
# 绘制火山图并添加标注
ggplot(data, aes(x = log2FoldChange, y = negLogP, color = Color,shape = group)) +
  geom_point(alpha = 0.8, size = 4) +
  xlim(-15,10)+ylim(0,8)+
  #scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  scale_color_identity() +  
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10(pvalue)",
    color = "Significance"
  )+
  theme(    # 移除网格线
    panel.border = element_rect(     # 添加四周边框
      colour = "black",              # 边框颜色
      fill = NA,                     # 填充设为空
      size = 1  ),                     # 边框线宽
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    panel.grid = element_blank(),            # 去除背景网格
    panel.background = element_blank(),
    legend.text = element_text(size = 15),   # 修改图例文字字体大小
    legend.title = element_text(size = 15)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", size = 0.5) +
  geom_text_repel(
    data = top_points,
    aes(label = key, color = Color),  # 确保 Color 列中包含正确的颜色分类
    size = 4,max.overlaps = 50
  ) +
  geom_segment(
    data = top_points, 
    aes(x = log2FoldChange, y = negLogP, xend = log2FoldChange, yend = negLogP + 0.5, color = Color),  # 设置连线的起始点和结束点
    size = 0.5, linetype = "solid"
  )


 
