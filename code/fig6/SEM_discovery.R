library(lavaan)
library(dplyr)
library(tibble)
metab = read.delim('./metab/metab.txt',check.names = F)
metabolites <- c(
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
metab = metab[,c('SampleID',metabolites)]
colnames(metab)[1] = 'Sample'
#-------------------------------------------------------
micro = read.delim('./data/属水平绝对定量.txt',check.names = F)
fungi = c('Trichosporon','Aspergillus','Saccharomyces','Candida',
            'Malassezia')
bac = c('Akkermansia','Alistipes','Enterocloster',
         'Faecalibacterium','Roseburia','Blautia','Mediterraneibacter',
       'Intestinibacter',
          'Ligilactobacillus','Bacteroides','Bifidobacterium','Veillonella',
       'Staphylococcus','Lactobacillus','Mycobacterium','Klebsiella','Escherichia',
          'Enterobacter','Streptococcus',
          'Enterococcus')
names = c(fungi,bac)
micro = micro[micro$Genus %in% names,]
fungi_names = micro[micro$Kingdom == 'Eukaryota',]$Species
bacteria_names = micro[micro$Kingdom == 'Bacteria',]$Species
#micro = micro[micro$Species %in% micros,]
micro = micro[,-c(2,3)]
micro = as.data.frame(t(micro))
colnames(micro) = micro[1,]
micro = micro[-1,]
micro <- rownames_to_column(micro, var = "Sample")
micro <- micro %>%
  mutate(across(.cols = 2:ncol(micro), .fns = as.numeric))
micro = micro[micro$Sample %in% metab$Sample,]
#----------------------------------------------------------------
metadata = read.delim('./data/metadata.txt',check.names = F)
metadata = metadata[,c(1,4)]
metadata = metadata[metadata$SampleID %in% micro$Sample,]
colnames(metadata)[1] = 'Sample'
#----------------------------------------------------------------







# 微生物丰度表
microbiome_data <-micro
# 代谢物丰度表
metabolite_data <-metab 
# 感染状态数据
infection_status <- metadata

# 合并所有数据
data <- merge(microbiome_data, metabolite_data, by = "Sample")
data <- merge(data, infection_status, by = "Sample")
data0 = data
data0 <- data0 %>%
  mutate(Group = recode(Group, "W/ Infection" = 0, "W/O Infection" = 1))

data0[,2:(ncol(data0)-1)] = log1p(data0[,2:(ncol(data0)-1)])
#write.table(data,'test0.csv',sep=',',row.names = F)
# 结构方程模型公式
results <- data.frame(bacteria = character(0),
                      fungi = character(0),
                      metabolites = character(0),
                      moderate_effect_p = numeric(0),
                      indirect_effect_p = numeric(0))
bacteria_names0 = c(  'Veillonella_parvula',
                      'Bifidobacterium_longum',
                      'Faecalibacterium_prausnitzii',
                      'Klebsiella_pneumoniae',
                      'Escherichia_coli',
                      'Streptococcus_salivarius'
                      )
total_iterations <- length(bacteria_names) * length(fungi_names) * length(metabolites)
for (i in bacteria_names){
  #i = 'Lactobacillus_acidophilus'
  for (j in fungi_names){
    ##j = 'Candida_albicans'
    for (x in metabolites){
     # x = 'Acetylcholine'
      tmp0 = data0[,c(i,j,x,'Group')]
      colnames(tmp0) = c('bacteria','fungi','metabolites','Group')
      model <- '
  
    metabolites ~ c1*bacteria + c2*fungi + c3*bacteria:fungi
    Group ~ c4*metabolites + c5*bacteria + c6*fungi + c7*bacteria:fungi
  
    indirect_effect := c1 * c4  # 计算菌a通过代谢物b对感染状态的间接效应
  
    moderate_effect := c3 * c4  # 菌a对感染状态的直接效应，考虑菌b的调节作用
  '
      fit <- sem(model, data = tmp0)
      summary(fit)
      param_estimates <- parameterEstimates(fit)
      
      # 提取moderate_effect的p值
      moderate_effect_p_value <- param_estimates[param_estimates$label == "moderate_effect", "pvalue"]
      indirect_effect_p_value <- param_estimates[param_estimates$label == "indirect_effect", "pvalue"]
      # 将结果保存到data frame中
      results <- rbind(results, data.frame(bacteria = i, fungi = j, metabolites = x,
                                           moderate_effect_p = moderate_effect_p_value,
                                           indirect_effect_p = indirect_effect_p_value))
    
      
      
    }
  }
}
#-------------------------------------------------------------------
model <- '

  Citric acid ~ c1*Faecalibacterium_prausnitzii + c2*Malassezia_restricta + 
  c3*Faecalibacterium_prausnitzii:Malassezia_restricta
  
  Group ~ c4*Citric acid + c5*Faecalibacterium_prausnitzii
  + c6*Malassezia_restricta + c7*Faecalibacterium_prausnitzii:Malassezia_restricta

  indirect_effect := c1 * c4  # 计算菌a通过代谢物b对感染状态的间接效应
  

  moderate_effect := c3 * c4  # 菌a对感染状态的直接效应，考虑菌b的调节作用
  
  # 总效应
  total_effect := indirect_effect + moderate_effect + c5  # 总效应包括直接效应和间接效应
'
# 拟合模型
library(blavaan)
fit <- sem(model, data = data0)

# 输出模型结果
summary(fit)
# library(semPlot)
semPaths(fit,whatLabels = 'est',layout = 'tree',style = 'lisrel')
save(fit,file = './metab/SEM.RData')
