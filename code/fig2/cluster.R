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

fun0 = fun[rownames(fun) %in% fun_names,]
bac0 = bac[rownames(bac) %in% bac_names,]

fun0 = t(fun0)
bac0 = t(bac0)

cor_matrix <- cor(fun0, bac0, method = "spearman")
brks <- c(seq(-0.6, -0.01, by = 0.001), seq(0, 0.6, by = 0.001))
# 绘制热图
pheatmap(cor_matrix,
         #color = colorRampPalette(c("#2171B5", "white", "#D94801"))(100), # 配色方案
         main = "Fungi-Bacteria Interaction Heatmap", # 热图标题
         color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
                   colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
         legend_breaks = seq(-0.6, 0.6, 0.1),
         breaks = brks,
         cluster_rows = TRUE, # 行聚类
         cluster_cols = TRUE) # 列聚类
#-----------------------------------------------------------------------


cor_test_results <- psych::corr.test(fun0, bac0, method = "spearman", adjust = "none") 

cor_matrix <- cor_test_results$r   # 相关性矩阵
p_matrix <- cor_test_results$p      # p 值矩阵

# 设置不显著的相关性为 NA（例如 p > 0.05）
cor_matrix[p_matrix > 0.05] <- 0

# 设定颜色映射，使得 0 为白色
brks <- c(seq(-0.6, -0.01, by = 0.001), seq(0, 0.6, by = 0.001))

# 绘制热图
pheatmap(cor_matrix,
         main = "Fungi-Bacteria Interaction Heatmap",  # 热图标题
         color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
                   colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
         legend_breaks = seq(-0.6, 0.6, 0.1),
         breaks = brks,
         cluster_rows = TRUE,  # 行聚类
         cluster_cols = TRUE)   # 使得不显著的值显示为白色










#----------------------------------------------------------------------
pheatmap(cor_matrix,
         clustering_distance_rows = "manhattan", # 行距离：曼哈顿距离
         clustering_distance_cols = "manhattan", # 列距离：曼哈顿距离
         clustering_method = "average",         # 平均聚类方法
         main = "Manhattan + Average Clustering")

# # 显著性计算
# 
# p_values <- matrix(NA, nrow = ncol(fun0), ncol = ncol(bac0))
# for (i in 1:ncol(fun0)) {
#   for (j in 1:ncol(bac0)) {
#     test <- cor.test(fun0[, i], bac0[, j], method = "spearman")
#     p_values[i, j] <- test$p.value
#   }
# }
# 
# # 筛选显著性（假设显著性水平为 0.05）
# significance <- ifelse(p_values < 0.05, "*", "")
# 
# # 绘制带标注的热图
# pheatmap(cor_matrix,
#          color = colorRampPalette(c("blue", "white", "red"))(100),
#          main = "Fungi-Bacteria Interaction Heatmap",
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          display_numbers = significance) # 添加显著性标注
fungi_1 = c('Ramularia','Trichosporon','Clavispora','Diutina',
            'Malassezia','Geosmithia',
            'Pichia','Neurpspora','Pyricularia')
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

names = c(fun_names,bac_names)

#--------------------------------------------------------------------------------
bac_da = read.delim('./fig2/volcano1/细菌genus热图数据.txt')
fun_da = read.delim('./fig2/volcano1/真菌genus热图数据.txt')

# bac_da0 = bac_da[,c(1,10:16)]
# fun_da0 = fun_da[,c(1,10:16)]

bac_da1 = bac_da[bac_da$key %in% bac_names,]
fun_da1 = fun_da[fun_da$key %in% fun_names,]


fun_wt01 = fun_da1[,c(1,10,11)]
fun_wt01$key <- factor(fun_wt01$key, levels = fungi_order)
fun_wt01 <- reshape2::melt(fun_wt01, id.vars = "key", 
                          variable.name = "Group", 
                          value.name = "Median")

fun_wt01 <- fun_wt01 %>%
  mutate(Median = ifelse(Group == "median_wt0", -Median, Median))

fun_wt01$ValueType <- ifelse(fun_wt01$Median >= 0, "W_T1", "W_T0")
mapping <- list(
  m1 = fungi_1,
  m2 = fungi_2,
  m3 = fungi_3,
  m4 = fungi_4
)
fun_wt01$Module <- 0
for (name in names(mapping)) {
  fun_wt01$Module[fun_wt01$key %in% mapping[[name]]] <- name
}
fun_wt01_top8 <- fun_wt01 %>%
  group_by(key) %>%
  summarise(total_abundance = sum(abs(Median))) %>%
  top_n(8, total_abundance) %>%
  inner_join(fun_wt01, by = "key")
# 假设数据框为 fun_wt01，并已正确设置
ggplot(fun_wt01_top8, aes(x = Median, y = key, fill = Module, alpha = ValueType)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_manual(values = c("m1" = "#DBB428", "m2" = "#D4562E",
                               "m3" = "#FFFFB3", "m4" = "#FCCDE5")) +  # 自定义颜色
  scale_alpha_manual(values = c("W_T1" = 1, "W_T0" = 0.7)) +  # 设置透明度
  labs(x = "Median", y = "Genus", fill = "Model", alpha = "Group") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank(),            # 去除背景网格
        panel.background = element_blank(),
        legend.text = element_text(size = 10),   # 修改图例文字字体大小
        legend.title = element_text(size = 10))      # 去除背景色
#------------------------------------------------------

#----------------------------------------------------------------------------------

bac_wt01 = bac_da1[,c(1,10,11)]
bac_wt01$key <- factor(bac_wt01$key, levels = bacteria_order)
bac_wt01 <- reshape2::melt(bac_wt01, id.vars = "key", 
                           variable.name = "Group", 
                           value.name = "Median")

#bac_wt01$Median = log2(bac_wt01$Median)
bac_wt01 <- bac_wt01 %>%
  mutate(Median = ifelse(Group == "median_wt0", -Median, Median))

bac_wt01$ValueType <- ifelse(bac_wt01$Median >= 0, "W_T1", "W_T0")
mapping <- list(
  m1 = bac_1,
  m2 = bac_2,
  m3 = bac_3,
  m4 = bac_4
)
bac_wt01$Module <- 0
for (name in names(mapping)) {
  bac_wt01$Module[bac_wt01$key %in% mapping[[name]]] <- name
}
bac_wt01_top8 <- bac_wt01 %>%
  group_by(key) %>%
  summarise(total_abundance = sum(abs(Median))) %>%
  top_n(8, total_abundance) %>%
  inner_join(bac_wt01, by = "key")
# 假设数据框为 fun_wt01，并已正确设置
ggplot(bac_wt01_top8, aes(x = Median, y = key, fill = Module, alpha = ValueType)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  # scale_x_break(c(-1e+07,-3e+06),#截断位置及范围
  #                space = 0.05,#间距大小
  #                scales = 1.5)+
  # scale_x_break(c(-1e+05,-1e+03),
  #               scales = 1.5,
  #               space = 0.3)+#上下显示比例，大于1上面比例大，小于1下面比例大 
  scale_fill_manual(values = c("m1" = "#A51C36", "m2" = "#7ABBDB",
                               "m3" = "#84BA42", "m4" = "#682487")) +  # 自定义颜色
  scale_alpha_manual(values = c("W_T1" = 1, "W_T0" = 0.7)) +  # 设置透明度
  labs(x = "Median", y = "Genus", fill = "Model", alpha = "Group") +
  theme(axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank(),            # 去除背景网格
        panel.background = element_blank(),
        legend.text = element_text(size = 10),   # 修改图例文字字体大小
        legend.title = element_text(size = 10))     # 去除背景色
#-----------------------------------------------------------------
bacteria_order <- c(
  "Streptococcus", "Veillonella", "Bifidobacterium", "Klebsiella", "Rothia",
  "Enterococcus", "Schaalia", "Blautia", "Escherichia", "Haemophilus",
  "Mediterraneibacter", "Enterobacter", "Bacteroides", "Enterocloster", 
  "Faecalibacterium", "Citrobacter", "Mycobacterium", "Phocaeicola",
  "Thomasclavelia", "Roseburia", "Staphylococcus", "Lactobacillus",
  "Anaerobutyricum", "Ligilactobacillus", "Pseudomonas", "Hungatella",
  "Parabacteroides", "Lacticaseibacillus", "Yersinia", "Alistipes", 
  "Segatella", "Weissella", "Akkermansia", "Megamonas", "Abiotrophia", 
  "Raoultella", "Intestinibacter", "Morganella", "Kluyvera",
  "Stenotrophomonas", "Sarcina", "Dysgonomonas"
)


#-----------------------------------------------------------------------------------

group = group[group$Timepoint %in% c('T0','T1','T2','T3'),]
group0 = group[,c(1,10)]
rownames(group0) = group$SampleID
group0 = group0[,-1,drop=F]


# 按分组计算真菌和细菌的相关性矩阵
groups <- unique(group0$Timepoint0)  # 获取所有组别
correlation_results <- list()  # 存储每个分组的相关性矩阵

for (group in groups) {
  #group = 'W_T0'
  group_samples <- rownames(group0)[group0$Timepoint0 == group]
  
  # 提取真菌和细菌丰度矩阵
  fungi_group <- fun0[group_samples, , drop = FALSE]
  bacteria_group <- bac0[group_samples, , drop = FALSE]
  
  # 计算相关性矩阵（真菌 vs 细菌）
  cor_matrix <- cor(bacteria_group, method = "spearman")
  correlation_results[[group]] <- cor_matrix  # 保存结果
}

# 绘制热图
for (group in groups) {
  #group = 'WO_T0'
  cor_matrix <- correlation_results[[group]]
  
  # 绘制热图
  pheatmap(cor_matrix,
           width = 12,
           height = 8,
           main = paste("Correlation Heatmap -", group),
           filename = paste("./fig2/1229/corr/bac_corr_", group,'.pdf'),
           cluster_rows = TRUE,  # 是否对行聚类
           cluster_cols = TRUE,  # 是否对列聚类
           #display_numbers = TRUE,  # 在热图上显示相关系数
           fontsize_number = 7)  # 数字字体大小
}
