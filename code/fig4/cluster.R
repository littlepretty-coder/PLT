library(tidyr)
library(pheatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
group = read.delim('./data/metadata.txt',sep='\t')
group$Timepoint0 <- ifelse(group$Group == "W/ Infection", 
                           paste0("W_", group$Timepoint), 
                           ifelse(group$Group == "W/O Infection", 
                                  paste0("WO_", group$Timepoint), 
                                  group$Timepoint))



OTUs <- read.delim('./data/属水平绝对定量.txt',sep='\t',check.names = F)
bacteria = OTUs[OTUs$Kingdom == 'Bacteria',]
fungi = OTUs[OTUs$Kingdom == 'Eukaryota',]
#archaea = OTUs[OTUs$Kingdom == 'Archaea',]
bac  = bacteria
bac = bac[,-c(1,3)]
bac = aggregate(. ~ Genus, data = bac, FUN = sum)
rownames(bac) = bac$Genus
bac = bac[,-1]

fun  = fungi
fun = fun[,-c(1,3)]
fun = aggregate(. ~ Genus, data = fun, FUN = sum)
rownames(fun) = fun$Genus
fun = fun[,-1]

bac_names = c("Anaerobutyricum", "Schaalia", "Segatella", "Mediterraneibacter",
        "Phocaeicola", "Intestinibacter", "Alistipes", "Escherichia",
        "Enterobacter", "Morganella", "Dysgonomonas", "Bifidobacterium",
        "Raoultella", "Bacteroides", "Lacticaseibacillus", "Lactobacillus",
        "Sarcina", "Enterocloster", "Weissella", "Akkermansia", "Hungatella",
        "Pseudomonas", "Enterococcus", "Roseburia", "Faecalibacterium",
        "Ligilactobacillus", "Thomasclavelia", "Citrobacter", "Megamonas",
        "Rothia", "Blautia", "Abiotrophia", "Parabacteroides", "Staphylococcus",
        "Mycobacterium", "Streptococcus", "Kluyvera", "Klebsiella", "Yersinia",
        "Haemophilus", "Veillonella", "Stenotrophomonas")

fun_names = c("Fusarium", "Malassezia", "Geosmithia", "Uncinocarpus",
              "Meyerozyma", "Kwoniela", "Saccharomyces", "Puccinia", "Botrytis",
              "Pyricularia", "Ramularia", "Trichosporon", "Mycotypha", "Letharia",
              "Colletotrichum", "Neurospora", "Candida", "Rhizophagus",
              "Phycomyces", "Suillus", "Talaromyces", "Tetrapisispora",
              "Alternaria", "Xylona", "Pichia", "Penicillium",
              "unclassified_Saccharomycetales", "Lodderomyces", "Aspergillus",
              "Nakaseomyces", "Drepanopeziza", "Dissoconium", "Diutina",
              "Kluyveromyces", "Clavispora", "Meira")
names = c(bac_names,fun_names)
fun0 = fun[rownames(fun) %in% fun_names,]
bac0 = bac[rownames(bac) %in% bac_names,]

fun0 = t(fun0)
bac0 = t(bac0)

cor_matrix <- cor(fun0, bac0, method = "spearman")
brks <- c(seq(-0.6, -0.01, by = 0.001), seq(0, 0.6, by = 0.001))
# 绘制热图
pheatmap(cor_matrix,
         color = colorRampPalette(c("#2171B5", "white", "#D94801"))(100), # 配色方案
         main = "Fungi-Bacteria Interaction Heatmap", # 热图标题
         # color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
         #           colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
         # legend_breaks = seq(-0.6, 0.6, 0.1),
         # breaks = brks,
         cluster_rows = TRUE, # 行聚类
         cluster_cols = TRUE) # 列聚类
#-----------------------------------------------------------------------

pheatmap_result = pheatmap(cor_matrix,
                           color = colorRampPalette(c("#2171B5", "white", "#D94801"))(100), # 配色方案
                           main = "Fungi-Bacteria Interaction Heatmap", # 热图标题
                           # color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
                           #           colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
                           # legend_breaks = seq(-0.6, 0.6, 0.1),
                           # breaks = brks,
                           cluster_rows = TRUE, # 行聚类
                           cluster_cols = TRUE) # 列聚类
# 提取行的聚类分簇
row_clusters <- cutree(pheatmap_result$tree_row, k = 4)  # 假设要将行聚类成3个簇

# 提取列的聚类分簇
col_clusters <- cutree(pheatmap_result$tree_col, k = 4)  # 假设要将列聚类成3个簇

# 查看聚类分簇结果
print(row_clusters)
print(col_clusters)



fungi_1 = c('Ramularia','Trichosporon','Clavispora','Diutina',
            'Malassezia','Geosmithia',
            'Pichia','Neurospora','Pyricularia')
fungi_2 = c('Botrytis','Rhizophagus','Aspergillus','Fusarium',
            'Suillus','Colletotrichum','Penicillium','Letharia',
            'Talaromyces','Puccinia','Tetrapisispora')
fungi_3 = c('Dissoconium','Saccharomyces','Xylona','Meyerozyma',
            'Nakaseomyces','Kluyveromyces','Lodderomyces')
fungi_4 = c('Candida','Uncinocarpus','Drepanopeziza','Alternaria',
            'Mycotypha','Phycomyces','Meira','unclassified_Saccharomycetales')

bac_1 = c('Sarcina','Akkermansia','Segatella','Alistipes','Enterocloster',
          'Hungatella','Faecalibacterium','Roseburia','Blautia','Mediterraneibacter',
          'Anaerobutyricum','Intestinibacter','Dysgonomonas','Weissella',
          'Ligilactobacillus')
bac_2 = c('Thomasclavelia','Bacteroides','Phocaeicola','Parabacteroides')
bac_3 = c('Pseudomonas','Kluyvera','Stenotrophomonas','Schaalia',
          'Bifidobacterium','Veillonella','Morganella','Megamonas',
          'Staphylococcus','Lactobacillus','Abiotrophia','Haemophilus',
          'Lacticaseibacillus')
bac_4 = c('Mycobacterium','Klebsiella','Yersinia','Escherichia',
          'Raoultella','Citrobacter','Enterobacter','Streptococcus',
          'Enterococcus','Rothia')
fun_names = c(fungi_1,fungi_2,fungi_3,fungi_4)
bac_names = c(bac_1,bac_2,bac_3,bac_4)

all = rbind(fun0,bac0)

all = na.omit(all)
names = c(fun_names,bac_names)
mapping <- list(
  fungi_1 = fungi_1,
  fungi_2 = fungi_2,
  fungi_3 = fungi_3,
  fungi_4 = fungi_4,
  bacteria_1 = bac_1,
  bacteria_2 = bac_2,
  bacteria_3 = bac_3,
  bacteria_4 = bac_4
)

# ================================= 开始重抽样 ==========================================


# ================================= 真菌-细菌相关性聚类的一致性评估 ==========================================

library(mclust)
library(aricode)
library(clusteval)

set.seed(123)
n_iterations <- 50

# 1. 获取真菌和细菌的属名
fungi_genera <- colnames(fun0)  # 真菌属名
bacteria_genera <- colnames(bac0)  # 细菌属名

cat("真菌属数量:", length(fungi_genera), "\n")
cat("细菌属数量:", length(bacteria_genera), "\n")

# 2. 创建您的自定义分组（基于生物学知识）
# 真菌的分组（您之前定义的）
fungi_groups <- list(
  fungi_1 = fungi_1,
  fungi_2 = fungi_2,
  fungi_3 = fungi_3,
  fungi_4 = fungi_4
)

# 细菌的分组（您之前定义的）
bacteria_groups <- list(
  bacteria_1 = bac_1,
  bacteria_2 = bac_2,
  bacteria_3 = bac_3,
  bacteria_4 = bac_4
)

# 3. 为每个菌属分配原始分组标签
original_fungi_clusters <- c()
for(i in 1:length(fungi_groups)) {
  for(genus in fungi_groups[[i]]) {
    if(genus %in% fungi_genera) {
      original_fungi_clusters[genus] <- i
    }
  }
}

original_bacteria_clusters <- c()
for(i in 1:length(bacteria_groups)) {
  for(genus in bacteria_groups[[i]]) {
    if(genus %in% bacteria_genera) {
      original_bacteria_clusters[genus] <- i
    }
  }
}

original_row_clusters <- original_fungi_clusters
original_col_clusters <- original_bacteria_clusters

cat("\n原始真菌聚类结果:\n")
print(table(original_row_clusters))

cat("\n原始细菌聚类结果:\n")
print(table(original_col_clusters))

# ========== 开始Bootstrap重抽样 ==========
cat("\n开始Bootstrap重抽样...\n")
cat("每次抽取", n_samples, "个样本（有放回）\n")
pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

# 存储每次bootstrap的聚类结果
bootstrap_row_clusters <- list()
bootstrap_col_clusters <- list()
bootstrap_cor_matrices <- list()

# 存储一致性结果
consistency_results <- data.frame(
  Iteration = 1:n_iterations,
  Row_ARI = NA,      # 真菌聚类一致性
  Col_ARI = NA, 
  Row_Jaccard = NA  ,
  Row_RI = NA  ,# 细菌聚类Jaccard# 细菌聚类一致性
  Row_NMI = NA,
  Col_NMI = NA,
  Col_Jaccard = NA ,
  Col_RI = NA  # 细菌聚类Jaccard
)

for(i in 1:n_iterations) {
  # ========== 1. 对样本进行重抽样 ==========
  sampled_samples <- sample(sample_names, size = n_samples, replace = TRUE)
  
  # 处理重复样本名
  if(any(duplicated(sampled_samples))) {
    unique_samples <- make.unique(sampled_samples)
  } else {
    unique_samples <- sampled_samples
  }
  
  # ========== 2. 提取重抽样后的真菌和细菌数据 ==========
  # 注意：需要根据原始样本名找到对应的列
  bootstrap_fungi <- matrix(0, nrow = length(fungi_genera), ncol = n_samples)
  rownames(bootstrap_fungi) <- fungi_genera
  colnames(bootstrap_fungi) <- unique_samples
  
  bootstrap_bacteria <- matrix(0, nrow = length(bacteria_genera), ncol = n_samples)
  rownames(bootstrap_bacteria) <- bacteria_genera
  colnames(bootstrap_bacteria) <- unique_samples
  
  for(j in 1:n_samples) {
    original_sample <- sampled_samples[j]  # 原始样本名
    new_sample <- unique_samples[j]        # 新的列名
    
    # 从原始数据中提取该样本的数据
    if(original_sample %in% colnames(all)) {
      bootstrap_fungi[, j] <- as.numeric(all[fungi_genera, original_sample])
      bootstrap_bacteria[, j] <- as.numeric(all[bacteria_genera, original_sample])
    }
  }
  
  # ========== 3. 计算真菌-细菌相关性矩阵 ==========
  # 注意：这里计算的是真菌和细菌之间的相关性，与原始pheatmap一致
  bootstrap_cor <- cor(t(bootstrap_fungi), t(bootstrap_bacteria), 
                       method = "spearman", 
                       use = "pairwise.complete.obs")
  heatmap_result00 = pheatmap(bootstrap_cor)
  bootstrap_cor_matrices[[i]] <- bootstrap_cor
  
  # ========== 4. 对相关性矩阵进行聚类（与pheatmap完全一致） ==========
  # 对行（真菌）进行聚类
  # bootstrap_row_dist <- dist(bootstrap_cor, method = "euclidean")
  # bootstrap_row_hc <- hclust(bootstrap_row_dist, method = "complete")
  bootstrap_row_clusters[[i]] <- cutree(heatmap_result00$tree_row, k = 4) 
  
  # 对列（细菌）进行聚类
  # bootstrap_col_dist <- dist(t(bootstrap_cor), method = "euclidean")
  # bootstrap_col_hc <- hclust(bootstrap_col_dist, method = "complete")
  bootstrap_col_clusters[[i]] <- cutree(heatmap_result00$tree_col, k = 4)
  
  # ========== 5. 计算与原始聚类的一致性 ==========
  # 真菌聚类一致性
  common_rows <- intersect(names(original_row_clusters), 
                           names(bootstrap_row_clusters[[i]]))
  if(length(common_rows) > 1) {
    consistency_results$Row_ARI[i] <- adjustedRandIndex(
      original_row_clusters[common_rows],
      bootstrap_row_clusters[[i]][common_rows]
    )
    consistency_results$Row_NMI[i] <- NMI(
      original_row_clusters[common_rows],
      bootstrap_row_clusters[[i]][common_rows]
    )
    consistency_results$Row_RI[i] <- RI(
      original_row_clusters[common_rows],
      bootstrap_row_clusters[[i]][common_rows]
    )
    consistency_results$Row_Jaccard[i] <- external_validation( original_row_clusters[common_rows], 
                                                               bootstrap_row_clusters[[i]][common_rows], method = "jaccard_index") 
  }
  
  # 细菌聚类一致性
  common_cols <- intersect(names(original_col_clusters), 
                           names(bootstrap_col_clusters[[i]]))
  if(length(common_cols) > 1) {
    consistency_results$Col_ARI[i] <- adjustedRandIndex(
      original_col_clusters[common_cols],
      bootstrap_col_clusters[[i]][common_cols]
    )
    consistency_results$Col_NMI[i] <- NMI(
      original_col_clusters[common_cols],
      bootstrap_col_clusters[[i]][common_cols]
    )
    consistency_results$Col_RI[i] <- RI(
      original_col_clusters[common_cols],
      bootstrap_col_clusters[[i]][common_cols]
    )
    consistency_results$Col_Jaccard[i] <- external_validation( original_col_clusters[common_cols], 
                         bootstrap_col_clusters[[i]][common_cols], method = "jaccard_index") 
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

# ================================= 结果分析 ==========================================

cat("\n\n========== 真菌-细菌相关性聚类一致性评估 ==========\n")

cat("\n--- 真菌聚类一致性（行聚类）---\n")
cat("ARI (均值 ± SD):", 
    round(mean(consistency_results$Row_ARI, na.rm = TRUE), 3), "±",
    round(sd(consistency_results$Row_ARI, na.rm = TRUE), 3), "\n")
cat("NMI (均值 ± SD):", 
    round(mean(consistency_results$Row_NMI, na.rm = TRUE), 3), "±",
    round(sd(consistency_results$Row_NMI, na.rm = TRUE), 3), "\n")
cat("Jaccard (均值 ± SD):", 
    round(mean(consistency_results$Row_Jaccard, na.rm = TRUE), 3), "±",
    round(sd(consistency_results$Row_Jaccard, na.rm = TRUE), 3), "\n")

cat("\n--- 细菌聚类一致性（列聚类）---\n")
cat("ARI (均值 ± SD):", 
    round(mean(consistency_results$Col_ARI, na.rm = TRUE), 3), "±",
    round(sd(consistency_results$Col_ARI, na.rm = TRUE), 3), "\n")
cat("NMI (均值 ± SD):", 
    round(mean(consistency_results$Col_NMI, na.rm = TRUE), 3), "±",
    round(sd(consistency_results$Col_NMI, na.rm = TRUE), 3), "\n")
cat("Jaccard (均值 ± SD):", 
    round(mean(consistency_results$Col_Jaccard, na.rm = TRUE), 3), "±",
    round(sd(consistency_results$Col_Jaccard, na.rm = TRUE), 3), "\n")

# 3. 可视化三个指标的比较
par(mfrow = c(2, 4))

# 真菌的三个指标
hist(consistency_results$Row_ARI, main = "Fungi ARI", 
     xlab = "ARI", col = "steelblue", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Row_ARI, na.rm = TRUE), col = "red", lwd = 2)

hist(consistency_results$Row_NMI, main = "Fungi NMI", 
     xlab = "NMI", col = "lightblue", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Row_NMI, na.rm = TRUE), col = "red", lwd = 2)

hist(consistency_results$Row_Jaccard, main = "Fungi Jaccard", 
     xlab = "Jaccard", col = "cornflowerblue", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Row_Jaccard, na.rm = TRUE), col = "red", lwd = 2)

hist(consistency_results$Row_RI, main = "Fungi Jaccard", 
     xlab = "Jaccard", col = "cornflowerblue", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Row_RI, na.rm = TRUE), col = "red", lwd = 2)

# 细菌的三个指标
hist(consistency_results$Col_ARI, main = "Bacteria ARI", 
     xlab = "ARI", col = "coral", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Col_ARI, na.rm = TRUE), col = "red", lwd = 2)

hist(consistency_results$Col_NMI, main = "Bacteria NMI", 
     xlab = "NMI", col = "peachpuff", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Col_NMI, na.rm = TRUE), col = "red", lwd = 2)

hist(consistency_results$Col_Jaccard, main = "Bacteria Jaccard", 
     xlab = "Jaccard", col = "lightsalmon", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Col_Jaccard, na.rm = TRUE), col = "red", lwd = 2)

hist(consistency_results$Row_RI, main = "Fungi Jaccard", 
     xlab = "Jaccard", col = "cornflowerblue", breaks = 20, xlim = c(0, 1))
abline(v = mean(consistency_results$Row_RI, na.rm = TRUE), col = "red", lwd = 2)
# 4. 创建箱线图比较三个指标
library(ggplot2)
library(tidyr)

# 准备绘图数据
plot_data <- consistency_results %>%
  dplyr::select(Row_ARI, Row_NMI, Row_Jaccard,Row_RI, Col_ARI, Col_NMI, Col_Jaccard, Col_RI) %>%
  pivot_longer(cols = everything(), 
               names_to = "Metric", 
               values_to = "Value") %>%
  mutate(
    Type = ifelse(grepl("Row", Metric), "Fungi", "Bacteria"),
    Metric = gsub("Row_|Col_", "", Metric)
  )

# 箱线图
plot_data$Metric = factor(plot_data$Metric,levels = c('RI','NMI','ARI','Jaccard'))
p <- ggplot(plot_data, aes(x = Metric, y = Value, fill = Type)) +
  # 箱线图
  geom_boxplot(alpha = 0.8, 
               outlier.shape = 21, 
               outlier.size = 2,
               outlier.fill = "gray50",
               outlier.alpha = 0.5,
               outlier.stroke = 0.3,
               width = 0.6,
               position = position_dodge(width = 0.8),
               lwd = 0.8,
               fatten = 2) +
  
  # 添加半透明的抖动点显示数据分布
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.8,
                                              jitter.height = 0),
              alpha = 0.2, 
              size = 1.2,
              aes(color = Type),
              show.legend = FALSE) +
  
  # 添加中位数标签
  # geom_text(data = median_data,
  #           aes(x = Metric, y = median_value + 0.03, 
  #               label = sprintf("%.2f", median_value),
  #               group = Type),
  #           position = position_dodge(width = 0.8),
  #           size = 3.5,
  #           fontface = "bold") +
  
  # 颜色设置
scale_fill_manual(values = c("Fungi" = "#4477AA", "Bacteria" = "#EE6677")) +
  scale_color_manual(values = c("Fungi" = "#4477AA", "Bacteria" = "#EE6677")) +
  
  # 设置Y轴：0到1，每隔0.2一个刻度
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
    expand = c(0.02, 0.02)  # 稍微减少边缘空白
  ) +
  
  # # 添加参考线（可选）
  # geom_hline(yintercept = c(0.5, 0.7, 0.9), 
  #            linetype = "dashed", 
  #            color = "gray70", 
  #            alpha = 0.5,
  #            linewidth = 0.3) +
  
  # 标签
  labs(
    title = "Comparison of Clustering Consistency Metrics",
    subtitle = paste0("Based on ", n_iterations, " Bootstrap iterations"),
    x = "Metric",
    y = "Consistency Value",
    fill = "Group",
    caption = "Higher values indicate better consistency"
  ) +
  
  # 主题美化
  theme_minimal(base_size = 12) +
  theme(
    # 标题
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", 
                              margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40",
                                 margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0.5, size = 9, color = "gray60",
                                margin = margin(t = 10)),
    
    # 轴
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11),
    axis.ticks = element_line(color = "gray70"),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 网格线
    # panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    # panel.background = element_rect(fill = "white", color = NA),
    # plot.background = element_rect(fill = "white", color = NA),
    
    # 图例
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(1, "cm"),
    legend.background = element_rect(fill = "white", color = "gray90", 
                                     linewidth = 0.3),
    legend.margin = margin(t = 10, b = 5, l = 10, r = 10),
    
    # 边框
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),
    
    # 箱线图特定
    strip.background = element_rect(fill = "gray95", color = "gray70"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  
  # 添加显著性标记（如果需要）
  # stat_compare_means(aes(group = Type), 
  #                    label = "p.signif", 
  #                    method = "wilcox.test",
  #                    hide.ns = TRUE) +
  
  # 添加背景色条区分不同指标（可选）
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           alpha = 0.03, fill = "gray80") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           alpha = 0.03, fill = "gray80")

# 显示图形
print(p)

print(p)

# 5. 计算指标间的相关性
cat("\n\n=== 指标间相关性 ===\n")
metric_cor <- cor(consistency_results[, c("Row_ARI", "Row_NMI", "Row_Jaccard",'Row_RI' ,
                                          "Col_ARI", "Col_NMI", "Col_Jaccard",'Col_RI')], 
                  use = "pairwise.complete.obs")
pheatmap(metric_cor,cellwidth = 10,cellheight = 10)
print(round(metric_cor, 3))

# 6. 为每个菌属计算Jaccard一致性（基于共现）
cat("\n\n=== 基于Jaccard的菌属共现稳定性 ===\n")

# 创建真菌共现矩阵
fungi_cooccurrence <- matrix(0, nrow = length(fungi_genera), ncol = length(fungi_genera))
rownames(fungi_cooccurrence) <- fungi_genera
colnames(fungi_cooccurrence) <- fungi_genera

for(i in 1:n_iterations) {
  if(!is.null(bootstrap_row_clusters[[i]])) {
    current_clusters <- bootstrap_row_clusters[[i]]
    
    # 对每个聚类内的菌属对，增加共现计数
    for(cluster_id in 1:4) {
      genera_in_cluster <- names(current_clusters[current_clusters == cluster_id])
      genera_in_cluster <- intersect(genera_in_cluster, fungi_genera)
      
      if(length(genera_in_cluster) > 1) {
        for(g1 in genera_in_cluster) {
          for(g2 in genera_in_cluster) {
            if(g1 != g2) {
              fungi_cooccurrence[g1, g2] <- fungi_cooccurrence[g1, g2] + 1
            }
          }
        }
      }
    }
  }
}

# 转换为Jaccard相似度
fungi_jaccard <- fungi_cooccurrence / (successful_iter * 2 - fungi_cooccurrence)

# 查看与某个菌属最常共现的菌属
example_genus <- fungi_genera[1]
top_cooccur <- sort(fungi_jaccard[example_genus, ], decreasing = TRUE)[1:5]
cat("\n与", example_genus, "最常共现的真菌属（Jaccard相似度）:\n")
print(top_cooccur)

# 7. 保存结果
write.csv(consistency_results, "bootstrap_clustering_consistency_with_jaccard.csv", row.names = FALSE)

