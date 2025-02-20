rm(list=ls())
options(stringsAsFactors = F)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
groups = read.delim("./data/metadata.txt", sep='\t')
#groups = groups[(groups$Group == 'W/O Infection'),]
#groups = groups[(groups$is20 == 0)&(groups$Group == 'W/O Infection'),]

OTUs <- read.delim('./data/属水平绝对定量.txt',sep='\t',check.names = F)


bacteria = OTUs[OTUs$Kingdom == 'Bacteria',]
fungi = OTUs[OTUs$Kingdom == 'Eukaryota',]
archaea = OTUs[OTUs$Kingdom == 'Archaea',]

index = read.delim('./fig1/1C/fungi.diversity.index.txt',sep='\t')
index = index[,c(1,2)]


df = fungi
df = df[,-c(1,3)]
df = aggregate(df[,c(2:ncol(df))],by = list(df$Genus),sum)
colnames(df)[1] = 'Genus'
rownames(df) = df$Genus
df = df[,-1]
df = df[,colnames(df)%in%groups$SampleID]
# 定义函数，输入为丰度表
find_dominant_genus <- function(abundance_table) {
  sample_names <- colnames(abundance_table)
  result_df <- data.frame(Sample = character(), Dominant_Genus = character(), stringsAsFactors = FALSE)
  for (sample in sample_names) {
    sample_data <- abundance_table[, sample]
    dominant_index <- which.max(sample_data)
    dominant_genus <- rownames(abundance_table)[dominant_index]
    result_df <- rbind(result_df, data.frame(Sample = sample, Dominant_Genus = dominant_genus))
  }
  return(result_df)
}
dominant_genera_df <- find_dominant_genus(df)
groups = merge(groups,dominant_genera_df,by.y = 'Sample',by.x = 'SampleID')
groups = merge(groups,index,by.y = 'sample',by.x = 'SampleID')
#groups$Dominant_Genus <- ifelse(groups$Shannon > 4, "diversity", groups$Dominant_Genus)

mycol=c("#A6CEE3", "#CCEBC5", "#FDBF6F", "#FFFF99", "#66c2A5", "#FDDAEC", "#BF5B17", "#666666", "#1B9E77",
        "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#1F78B4", "#B2DF8A", "#33A02C",
        "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#ded54c",
        "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4",
        "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33",
        "#F781BF", "#999999", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
        "#BC80BD", "#FFED6F")



ggplot(groups, aes(x = days, y = ID, color = Dominant_Genus.y)) +
  geom_point(size = 5, aes(fill = Dominant_Genus.y), stroke = 1.2) +
  #geom_text(aes(label = SampleID), vjust = -1, size = 3) +
  scale_color_manual(values=mycol) +  # 设置颜色调色板
  #scale_shape_manual(values = 1:50) +  # 设置形状
  theme_minimal() +
  labs(x = "Days", y = "Patient", title = "Sample Overview by Patient, Days, and Main Taxa") +
  theme(panel.grid = element_blank())


get_dominant_or_diversity <- function(category_column) {
  freq_table <- table(category_column)  # 统计各类属的频率
  dominant <- names(which.max(freq_table))  # 找出出现次数最多的属
  if (max(freq_table) > length(category_column) / 2) {
    return(dominant)  # 如果该属超过一半，返回该属
  } else {
    return("diversity")  # 否则返回 "diversity"
  }
}

result <- groups %>%
  group_by(ID) %>%
  summarise(
    bac = get_dominant_or_diversity(Dominant_Genus.x),
    fungi = get_dominant_or_diversity(Dominant_Genus.y)
  ) %>%
  ungroup()
table(result$fungi)
# Aspergillus      Botrytis       Candida    Clavispora     diversity    Geosmithia    Malassezia  Nakaseomyces    Neurospora 
# 3             1            16             3            59             3             1             3            17 
# Pyricularia     Ramularia Saccharomyces  Trichosporon 
# 1             2            11             9 


dominant_genera_bac = dominant_genera_df
dominant_genera_fungi = dominant_genera_df
save(dominant_genera_bac,dominant_genera_fungi,file = './fig1/1D/domiant.RData')
load('./fig1/1D/domiant.RData')
write.table(dominant_genera_bac,'./figure/suppdata/fig1/sample_dominant_genera_bac.csv',sep=',',row.names = F)
write.table(dominant_genera_fungi,'./figure/suppdata/fig1/sample_dominant_genera_fungi.csv',sep=',',row.names = F)
