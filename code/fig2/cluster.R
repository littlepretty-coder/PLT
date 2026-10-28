library(tidyr)
library(pheatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
bacteria = OTUs[OTUs$Kingdom == 'Bacteria',]
fungi = OTUs[OTUs$Kingdom == 'Eukaryota',]
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
pheatmap(cor_matrix,

         main = "Fungi-Bacteria Interaction Heatmap",
         color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
                   colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
         legend_breaks = seq(-0.6, 0.6, 0.1),
         breaks = brks,
         cluster_rows = TRUE, 
         cluster_cols = TRUE) 
#-----------------------------------------------------------------------


cor_test_results <- psych::corr.test(fun0, bac0, method = "spearman", adjust = "none") 

cor_matrix <- cor_test_results$r
p_matrix <- cor_test_results$p  

cor_matrix[p_matrix > 0.05] <- 0
brks <- c(seq(-0.6, -0.01, by = 0.001), seq(0, 0.6, by = 0.001))

pheatmap(cor_matrix,
         main = "Fungi-Bacteria Interaction Heatmap", 
         color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
                   colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
         legend_breaks = seq(-0.6, 0.6, 0.1),
         breaks = brks,
         cluster_rows = TRUE, 
         cluster_cols = TRUE) 




