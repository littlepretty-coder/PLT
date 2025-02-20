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

metab = read.delim('./metab/metab.txt',check.names = F)
fun0 = fun0[rownames(fun0)%in%metab$SampleID,]
bac0 = bac0[rownames(bac0)%in%metab$SampleID,]
rownames(metab) = metab$SampleID
metab = metab[,-c(1,ncol(metab))]

tmp = read.delim('./metab/wilcox.txt')
tmp$fdr = p.adjust(tmp$kw_pvalue)
tmp = tmp[tmp$fdr < 0.01,]
names = c(
  "Acetylcarnitine",
  "N-methyltryptamine",
  "Glycerophosphocholine",
  "Citric acid",
  "Hydroxyphenyllactic acid",
  "3,4-DiHydroxyhydrocinnamic acid",
  "7-Dehydrocholic acid",
  "5-Hydroxylysine",
  "γ-Aminobutyric acid",
  "N-methyl-aspartic acid",
  "Conjugated linoleic acids",
  "Dihomo-γ-linolenic acid",
  "5Z,8Z,11Z,14Z,17Z-Eicosapentaenoic acid",
  "Diethanolamine",
  "13-Hydroxy-4Z,7Z,10Z,14E,16Z,19Z-docosahexaenoic acid",
  "α-Linolenic acid",
  "Glycolithocholic acid-3-sulfate",
  "Acetylcholine",
  "γ-Linolenic acid",
  "Propionylcarnitine"
)
metab = metab[,colnames(metab) %in% names,]

group = merge(metab[,1:2],group[,c(1,10)],by.x = 'row.names',by.y = 'SampleID')
group = group[,-c(2:3)]
#wo0 = group[group$Timepoint0 == 'W_T3',]$Row.names

bac00 = bac0[rownames(bac0)%in%wo0,]
metab0 = metab[rownames(metab)%in%wo0,]
cor_matrix <- cor(fun0, metab, method = "spearman")

# 计算每对相关性的p值矩阵
p_value_matrix <- matrix(NA, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    p_value_matrix[i, j] <- cor.test(fun0[, i], metab[, j], method = "spearman")$p.value
  }
}

# 显著性标注（设定p值的阈值）
sig_matrix <- ifelse(p_value_matrix < 0.001, "***",
                     ifelse(p_value_matrix < 0.01, "**",
                            ifelse(p_value_matrix < 0.05, "*", '')))


breaks_values <- seq(-0.4, 0.4, length.out = 101)
# 绘制热图
pheatmap(cor_matrix,
         color = colorRampPalette(c("#2171B5", "white", "#D94801"))(100),
         main = "Fungi-Metabolites Interaction Heatmap",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #display_numbers = sig_matrix,
         fontsize = 10,breaks = breaks_values)
