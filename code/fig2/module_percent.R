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
all = fun0
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
all$Module <- 0
for (name in names(mapping)) {
  all$Module[rownames(all) %in% mapping[[name]]] <- name
}
all = all[,c(ncol(all),1:(ncol(all)-1))]
tmpp = aggregate(all[,2:ncol(all)],by = list(all$Module),sum)

group = read.delim('./data/metadata.txt',sep='\t')
group$Timepoint0 <- ifelse(group$Group == "W/ Infection",
                           paste0("W_", group$Timepoint),
                           ifelse(group$Group == "W/O Infection",
                                  paste0("WO_", group$Timepoint),
                                  group$Timepoint))
group = group[,c(1,10)]


a_long <- melt(tmpp, id.vars = "Group.1", variable.name = "SampleID", value.name = "Abundance")
b = group
# 2. 将样本和分组信息合并
a_long <- a_long %>%
  left_join(b, by = "SampleID")
colnames(a_long)[4] = 'Group'
a_long = as.data.frame(a_long)
# 3. 按分组和物种汇总丰度
library(data.table)
grouped_data <- a_long %>%
  group_by(Group, Group.1) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop")
custom_colors <- c(
  "fungi_1" = "#DBB428", "fungi_2" = "#D4562E",
  "fungi_3" = "#FFFFB3", "fungi_4" = "#FCCDE5",
  "bacteria_1" = "#A51C36", "bacteria_2" = "#7ABBDB",
  "bacteria_3" = "#84BA42", "bacteria_4" = "#682487"
)
# 4. 计算每个分组内物种的百分比
grouped_data = grouped_data[!(grouped_data$Group %in% c('WO_','W_')),]
grouped_data <- grouped_data %>%
  group_by(Group) %>%
  mutate(Percentage = TotalAbundance / sum(TotalAbundance) * 100)

ggplot(grouped_data, aes(x = Group, y = Percentage, fill = Group.1)) +
  geom_bar(stat = "identity", position = "stack") + # 堆叠柱状图
  labs(x = "Group", y = "Percentage (%)", fill = "Species", 
       title = "Percentage of Species in Each Group") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # 显示为百分比格式
  theme_minimal() +
  scale_fill_manual(values = custom_colors) +          
  theme(
    axis.text = element_text(size = 10),       # 坐标轴文字大小
    axis.title = element_text(size = 12),      # 坐标轴标题大小
    legend.text = element_text(size = 10),     # 图例文字大小
    legend.title = element_text(size = 12)     # 图例标题大小
  )
write.table(grouped_data,'./figure/suppdata/fig2/module_fungal_abundance_percentage.csv',sep=',',row.names = F)
