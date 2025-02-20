library(microeco)
library(fdrtool)
library(meconetcomp)
library(magrittr)
library(igraph)
library(ggplot2)

group = read.delim('./data/metadata.txt',sep='\t')
group$Timepoint0 <- ifelse(group$Group == "W/ Infection", 
                           paste0("W_", group$Timepoint), 
                           ifelse(group$Group == "W/O Infection", 
                                  paste0("WO_", group$Timepoint), 
                                  group$Timepoint))
group = group[,c(1,10)]
colnames(group)[2] = 'Group'
rownames(group) = group$SampleID


df <- read.delim('./data/属水平绝对定量.txt',sep='\t',check.names = F)
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
df = df[df$Genus %in% names,]
rownames(df) = df$Species
df = df[,-c(1,2,3)]
rownames(df) = gsub('-','_',rownames(df))
patterns <- c("\\[", "\\]", "\\(", "\\)", "\\.", "\\'", "\\/", "\\:", "\\=", "\\+")
rownames(df) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(df) )


bac = read.delim('./fig2/new/细菌数据.txt')
fun = read.delim('./fig2/new/真菌数据.txt')
info = rbind(bac,fun)
info$fdr = p.adjust(info$kw_pvalue)
info = info[info$key %in% rownames(df),]
info = info[info$fdr < 0.01,]
df = df[rownames(df)%in%info$key,]


tax = read.delim('./data/tax.txt')
rownames(tax) = tax$Species
rownames(tax) = gsub('-','_',rownames(tax))
patterns <- c("\\[", "\\]", "\\(", "\\)", "\\.", "\\'", "\\/", "\\:", "\\=", "\\+")
rownames(tax) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(tax) )
tax = tax[rownames(tax) %in% rownames(df),]
tax = tidy_taxonomy(tax)



mt <- microtable$new(sample_table = group, otu_table = df, tax_table = tax)

mt


#--------------------------group network---------------------------------------------------------
networks <- list()

# W_T0


tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T0")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T0 <- tmp

# W_T1

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T1")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T1 <- tmp

# W_T2

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T2")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T2 <- tmp

# W_T3

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T3")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T3 <- tmp

# WO_T0

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "WO_T0")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$WO_T0 <- tmp

# WO_T1

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "WO_T1")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$WO_T1 <- tmp

# WO_T3

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "WO_T3")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$WO_T3 <- tmp

#--------------------------------------------------------------------------------
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
vectors <- list(
  fungi_1 = fungi_1,
  fungi_2 = fungi_2,
  fungi_3 = fungi_3,
  fungi_4 = fungi_4,
  bac_1 = bac_1,
  bac_2 = bac_2,
  bac_3 = bac_3,
  bac_4 = bac_4
)
class <- do.call(rbind, lapply(names(vectors), function(name) {
  data.frame(Element = vectors[[name]], Source = name, stringsAsFactors = FALSE)
}))
class$Element = paste0('g__',class$Element)
colnames(class) = c('Genus','Module')
#-----------------------------------------------------------------------------------------


networks %<>%cal_module(undirected_method ="cluster_fast_greedy")
tmp <- cal_network_attr(networks)

networks %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

color = data.frame(module = c('bac_1','bac_2','bac_3','bac_4',
                              'fungi_1','fungi_2','fungi_3','fungi_4'))
#extended_palette <- colorRampPalette(brewer.pal(8, "Paired"))
color$color = c("#A51C36", "#7ABBDB", "#84BA42", "#682487", "#DBB428", 
                "#D4562E", "#FFFFB3",'#FCCDE5')

module_colors <- setNames(color$color, color$module)



for(i in names(networks)){
  #i='W_T0'
  nodes = networks[[i]]$res_node_table
  edges = networks[[i]]$res_edge_table
  nodes = merge(nodes,class,by = 'Genus')
  nodes = nodes[,c(2,1,3:ncol(nodes))]
  abun = as.data.frame(t(networks[[i]]$data_abund))
  abun = abun[rownames(abun)%in%nodes$name,]
  nodes$abun = rowSums(abun)
  nodes <- nodes %>% select(-module)
  nodes = merge(nodes,color,by.x = 'Module',by.y = 'module')
  nodes = nodes[,c(2,1,3:ncol(nodes))]
  colnames(nodes)[1] = 'Id'
  colnames(edges)[1:2] = c('source','target')
  write.table(edges, paste0('./fig2/1229/net/',i,'_network.edge_list.csv'), sep = ',', row.names = FALSE, quote = FALSE)
  write.table(nodes, paste0('./fig2/1229/net/',i,'_network.node_list.csv'), sep = ',', row.names = FALSE, quote = FALSE)
  #networks[[i]]$save_network(filepath =paste0(tmp, "/", i, ".gexf"))
  
}






save(networks,file = './fig2/1229/net/networks.RData')
networks[["W_T0"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["W_T1"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["W_T2"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["W_T3"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["WO_T0"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["WO_T1"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["WO_T3"]]$plot_network(method ="networkD3", node_color ="Kingdom")
#-----------------------------compare node----------------------------------------

tmp <- node_comp(networks, property = "name")

# 使用trans_venn类来计算节点的交集

tmp1 <- trans_venn$new(tmp, ratio = "numratio")

g1 <- tmp1$plot_venn(fill_color = FALSE)
g1 <- tmp1$plot_bar()
g1
ggsave("./fig2/net/node_overlap.pdf", g1, width = 15, height = 6)

# 使用常用的 jaccard distance 来量化网络间的差别

tmp$cal_betadiv(method = "jaccard")

tmp$beta_diversity$jaccard
jac_node = as.vector(unlist(tmp$beta_diversity$jaccard))
#---------------------------------compare edge-------------------------------------------------

tmp <- edge_comp(networks)

# 通过trans_venn产生交集

tmp1 <- trans_venn$new(tmp, ratio = "numratio")

g1 <- tmp1$plot_venn(fill_color = FALSE)
g1 <- tmp1$plot_bar()
ggsave("./fig2/net/edge_overlap.pdf", g1, width = 7, height = 6)

# 使用jaccard distance量化网络间在边上的差别

tmp$cal_betadiv(method = "jaccard")

tmp$beta_diversity$jaccard
jac_edge = as.vector(unlist(tmp$beta_diversity$jaccard))
#--------------------------------------------------------------------------------
library(tidyr)
library(ggplot2)
library(ggpubr)
tmpp = as.data.frame(cbind(jac_edge,jac_node))
df_long <- pivot_longer(tmpp, cols = everything(), 
                        names_to = "Group", values_to = "Value")

ggboxplot(df_long, x = "Group", y = "Value", fill = "Group") +
  stat_compare_means(method = "t.test", label = "p.format") + # 添加t检验的显著性p值
  labs(title = "Boxplot with Significance",
       x = "Group", y = "Value") +
  theme_minimal()



#-----------------------------common edge to a new network---------------------------------------------

tmp <- edge_comp(networks)

tmp1 <- trans_venn$new(tmp)

# 转换集合信息到新的 microtable 对象

tmp2 <- tmp1$trans_comm()

# 使用subset_network提取结果，这里提取三个网络的交集

# 使用 colnames(tmp2$otu_table) 来查找名字

Intersec_all <-subset_network(networks, venn = tmp2, name ="W_T0&W_T1&W_T2&W_T3&WO_T0&WO_T1&WO_T3")
#Intersec_all$tidy_dataset()
Intersec_all$plot_network(method ="networkD3", node_color ="Phylum")


nodes = Intersec_all$res_node_table
edges = Intersec_all$res_edge_table
abun = as.data.frame(t(Intersec_all$data_abund))
abun = abun[rownames(abun)%in%nodes$name,]
nodes$abun = rowSums(abun)
colnames(nodes)[1] = 'Id'
colnames(edges)[1:2] = c('source','target')
write.table(edges, paste0('./fig2/net/',i,'_network.edge_list.csv'), sep = ',', row.names = FALSE, quote = FALSE)
write.table(nodes, paste0('./fig2/net/',i,'_network.node_list.csv'), sep = ',', row.names = FALSE, quote = FALSE)



# 
# #-------------------------------------phylo distance----------------------------------------------
# 
# 
# 
# node_names <- unique(unlist(lapply(networks, function(x){colnames(x$data_abund)})))
# 
# # 复制数据集，然后过滤
# filter_soil_amp <- clone(mt)
# filter_soil_amp$otu_table <- filter_soil_amp$otu_table[node_names, ]
# filter_soil_amp$tidy_dataset()
# # 获得系统发生距离矩阵
# phylogenetic_distance_soil <- as.matrix(cophenetic(filter_soil_amp$phylo_tree))
# # 调用edge_node_distance类
# tmp <- edge_node_distance$new(network_list = soil_amp_network, dis_matrix = phylogenetic_distance_soil, label = c("+", "-"))
# # 差异检验
# tmp$cal_diff(method = "anova")
# # Fig.2c
# g1 <- tmp$plot(boxplot_add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
# ggsave("soil_amp_phylo_distance.pdf", g1, width = 7, height = 6)
# 
# 
# #----------module----------
# 
# tmp <- edge_node_distance$new(network_list = soil_amp_network, dis_matrix = phylogenetic_distance_soil,
#                               label = "+", with_module = TRUE, module_thres = 10)
# tmp$cal_diff(method = "anova")
# g1 <- tmp$plot(boxplot_add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")

#-------------------------------节点的拓扑学属性------------------------------

tmp <- "./fig2/1229/net/node_roles"

dir.create(tmp)

for(i in names(networks)){
  #i=1
  networks[[i]]$res_node_table %<>% .[!is.na(.$taxa_roles), ]
  
  g1 <- networks[[i]]$plot_taxa_roles(add_label = TRUE, add_label_group = c("Module hubs", "Connectors"), label_text_size = 4, label_text_italic = TRUE)
  
  ggsave(paste0(tmp, "/", i, ".pdf"), g1, width =7, height =5)
  
}

#-----------------------------不同网络的节点分类学属性------------------------------


networks_edgetax <- edge_tax_comp(networks, taxrank = "Genus", label = "+", rel =T)

# 过滤比例较低的类别

networks_edgetax <- networks_edgetax[apply(networks_edgetax, 1, mean) > 0.01, ]

# Fig.2d

g1 <- pheatmap::pheatmap(networks_edgetax, display_numbers = TRUE,cluster_cols = F)

ggsave("./fig2/1229/net/edge_tax_comp.genus.pdf", g1, width = 10, height =15)

#-----------------------------与非生物因素的关联----------------------------------

data = read.delim('./fig2/net/env.txt')
load('./fig2/1229/net/networks.RData')
tmp <- "./fig2/1229/net/module_eigen"

dir.create(tmp)


tmp0 <- clone(mt)
tmp0$sample_table = merge(tmp0$sample_table,data[,-c(2:5)] ,by = 'SampleID')
rownames(tmp0$sample_table) = tmp0$sample_table$SampleID


for(i in names(networks)){
  #i='W_T0'
  networks[[i]]$res_node_table <-  merge(networks[[i]]$res_node_table,class,by = 'Genus')
  networks[[i]]$res_node_table = networks[[i]]$res_node_table[,c(2,1,3:ncol(networks[[i]]$res_node_table))]
  networks[[i]]$res_node_table <- networks[[i]]$res_node_table %>%
    mutate(module = Module)
  rownames(networks[[i]]$res_node_table) = networks[[i]]$res_node_table$name
  }

# 创建 trans_env 对象来进行相关分析，选择一些环境因子
# tmp00 <- trans_env$new(dataset = tmp0, env_cols =2:40)
# #tmp00$data_env[['Group']] <- as.numeric(as.factor(tmp00$data_env[['Group']]) )
# tmp00$cal_diff(group = 'Group',method = 'KW')
# tmp00$plot_diff(dataset = tmp00, env_cols =3:40)
# head(tmp00$res_diff)
# tmp <- list()
# for(i in colnames(tmp00$data_env)[2:39]){
#   tmp[[i]] <- tmp00$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
# }
# plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 5))
# #ggsave('1.pdf',width = 20,height = 20)
# tmp00$cal_ordination(method = "dbRDA", use_measure = "bray")
# # show the orginal results
# tmp00$trans_ordination()
# tmp00$plot_ordination(plot_color = "Group")
# # the main results of RDA are related with the projection and angles between arrows
# # adjust the length of the arrows to show them better
# tmp00$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# # t1$res_rda_trans is the transformed result for plotting
# tmp00$plot_ordination(plot_color = "Group")





for(i in names(networks)){
  #i='WO_T0'
  # module eigengene 分析
  networks[[i]]$cal_eigen()
  
  tmp1 <- clone(mt)
  tmp1$sample_table = merge(tmp1$sample_table,data[,-c(2:5)] ,by = 'SampleID')
  rownames(tmp1$sample_table) = tmp1$sample_table$SampleID
  tmp1$sample_table %<>% base::subset(Group == i)
  
  tmp1$tidy_dataset()
  
  # 创建 trans_env 对象来进行相关分析，选择一些环境因子

  tmp2 <- trans_env$new(dataset = tmp1, env_cols =2:40)
  # tmp2$data_env[['Group']] <- as.numeric(as.factor(tmp2$data_env[['Group']]) )
  # tmp2$cal_diff(group = 'Group',method = 'KW')
  #tmp2$plot_diff(dataset = tmp1, env_cols =3:40)
  tmp2$cal_cor(add_abund_table = networks[[i]]$res_eigen,cor_method = 'spearman')
  #dev.new()
  g1 <- tmp2$plot_cor(cluster_ggplot = "row")
  #dev.off()
  #plot(g1)
  g1
  ggsave(paste0(tmp, "/", i, ".pdf"), g1, width = 15, height =5)
  
}



#---------------不同模块的物种组成和功能组成--------------------------------------

tmp_dir <- "./fig2/1229/net/module_composition"

dir.create(tmp_dir)
color = data.frame(Genus = names)
extended_palette <- colorRampPalette(brewer.pal(12, "Set1"))
color$color = extended_palette(nrow(color))
color <- rbind(color, data.frame(Genus = 'Others', color = 'grey', stringsAsFactors = FALSE))
genus_colors <- setNames(color$color, color$Genus)
for(i in names(networks)){
  #i='W_T0'
  # 为了进行完整演示，重新获取下 res_node_table
  
  #networks[[i]]$get_node_table()
  
  # trans_comm 转换节点信息表
  
  tmp <- networks[[i]]$trans_comm(use_col ="module", abundance =FALSE)
  
  # 转换成二元信息
  
  tmp$otu_table[tmp$otu_table >0] <- 1
  
  # 过滤一些物种数目太少的模块
  
  tmp1 <- networks[[i]]$res_node_table$module %>% table %>% .[. >= 5] %>% names
  
  tmp$sample_table <- tmp$sample_table[rownames(tmp$sample_table) %in% tmp1, ]
  
  tmp$tidy_dataset()
  
  tmp$cal_abund()
  
  # 使用trans_abund来进行可视化
  
  
  tmp2 <- trans_abund$new(tmp, taxrank = "Genus", ntaxa = 15)
  
  tmp2$data_abund$Sample %<>%factor(., levels =rownames(tmp$sample_table))
  
  # 显示比例信息
  
    g1 <- tmp2$plot_bar(xtext_type_hor =TRUE) +ylab("OTUs ratio (%)")+
      scale_fill_manual(values = genus_colors)
  g1
  ggsave(paste0(tmp_dir, "/taxa_", i, ".pdf"), g1, width =7, height = 5)
  
  # # 功能预测
  # 
  # tmp3 <- trans_func$new(tmp)
  # 
  # tmp3$cal_spe_func(prok_database ="FAPROTAX")
  # 
  # # 获取百分比数据
  # 
  # tmp3$cal_spe_func_perc(abundance_weighted =FALSE)
  # 
  # # 作图
  # 
  # g1 <- tmp3$plot_spe_func_perc(order_x = tmp$sample_names())
  # g1
  # ggsave(paste0(tmp_dir, "/trait_", i, ".pdf"), g1, width =9, height = 6)
  
}


#------------------------------网络子集---------------------------------------------


