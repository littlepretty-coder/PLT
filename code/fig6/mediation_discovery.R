# 加载必要的 R 包
if (!require("mediation")) install.packages("mediation")
if (!require("dplyr")) install.packages("dplyr")
library(mediation)
library(dplyr)
library(readxl)
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
# 将分组信息转换为因子
merged_data = read.delim('./data/种水平绝对定量.txt',check.names = F)
merged_data = merged_data[merged_data$Genus %in% names,]
group = read.delim('./data/metadata.txt')
group$Timepoint0 <- ifelse(group$Group == "W/ Infection", 
                           paste0("W_", group$Timepoint), 
                           ifelse(group$Group == "W/O Infection", 
                                  paste0("WO_", group$Timepoint), 
                                  group$Timepoint))
merged_data = merged_data[,-c(2,3)]
merged_data = as.data.frame(t(merged_data))
colnames(merged_data) = merged_data[1,]
merged_data = merged_data[-1,]

metab = read_excel('./metab/all.xlsx')
merged_data = merged_data[rownames(merged_data) %in% colnames(metab)[6:206],]

metab = metab[,-c(2:5)]
metab = as.data.frame(t(metab))
colnames(metab) = metab[1,]
metab = metab[-1,]
tmp0 = merge(metab,merged_data,by= 'row.names')
tmp0 = merge(tmp0,group[,c(1,4)],by.x= 'Row.names',by.y = 'SampleID')

tmp0$Group <- as.factor(ifelse(tmp0$Group == "W/ Infection", 1, 0))
colnames(tmp0)[1] = 'SampleID'
merged_data = tmp0
#bac1'#A51C36'good,bac2'#7ABBDB'bad,bac3'#84BA42'good,bac1'#682487'bad
#-------------------------- -------------------------------------------------------
colnames(merged_data) <- gsub("[^a-zA-Z0-9_]", "", colnames(merged_data))
merged_data[,c(2:(ncol(merged_data)-1))] <- data.frame(lapply(merged_data[,c(2:(ncol(merged_data)-1))], function(x) {
  if (is.character(x)) { # 如果是字符串类型
    as.numeric(x) # 转换为数值
  } else {
    x # 保留原来的值
  }
}))
colnames(merged_data) <- sapply(colnames(merged_data), function(name) {
  if (grepl("^[0-9]", name)) { # 检查是否以数字开头
    paste0("X", name) # 添加前缀 "X"
  } else {
    name # 保持原列名
  }
})
#-------------------------------------------------------------------------------------

significant_results <- list()
#-----------------------------------------------------------------------------------------
bacteria_names = colnames(merged_data)[404:1666]
#metabolite_names = colnames(merged_data)[2:403]
metabolites <- c(
  "γ-Linolenic acid2", "α-Linolenic acid1", 
  "Linoleyl carnitine",
  "Phenol",
  "Propionylcarnitine",
  "Guanine",
  "Choline",
  "N-acetylglutamic acid",
  "Riboflavin-5'-monophosphate",
  "Methionine sulfoxide",
  "Glycerophosphocholine",
  'Glycolithocholic acid-3-sulfate','Propionylcarnitine','Acetylcholine',
  'Glycerophosphocholine','Lactose',
  'N-methyl-aspartic acid',
  'Choline'
)
metabolite_names <- gsub("[^a-zA-Z0-9_]", "", metabolites)
metabolite_names<- sapply(metabolite_names, function(name) {
  if (grepl("^[0-9]", name)) { # 检查是否以数字开头
    paste0("X", name) # 添加前缀 "X"
  } else {
    name # 保持原列名
  }
})
bacteria_names <- bacteria_names[bacteria_names != 'Group']
for (bacteria in bacteria_names) {
  #bacteria = 'Streptococcus_sp_FDAARGOS_521'
  for (metabolite in metabolite_names) {
    #metabolite = 'Propionylcarnitine'
    tmp = merged_data[,c('SampleID',bacteria,metabolite,'Group')]
    #colnames_before_scale <- colnames(tmp)
    
    tmp[,2:3]<- log1p(tmp[,2:3])
    #tmp[[metabolite]] <- scale(tmp[[metabolite]])
    
    # 恢复列名
    #colnames(tmp) <- colnames_before_scale

    if ( (sum(tmp[[bacteria]]) == 0)){
      next
    }
    mediator_model <- lm(as.formula(paste(metabolite, "~", bacteria)), data = tmp)
    mediator_summary = summary(mediator_model)

    # 拟合结果模型
    outcome_model <- glm(as.formula(paste("Group ~", bacteria, "+", metabolite)),
                         data = tmp, family = binomial(link = "logit"))
    outcome_summary = summary(outcome_model)

    # 进行中介效应分析
    mediation_result <- mediate(mediator_model, outcome_model,
                                treat = bacteria, mediator = metabolite,
                                boot = TRUE, sims = 1000)

    # 检查 p 值是否显著（中介效应的 p 值）
    if (summary(mediation_result)$d.avg.p<0.05) {
      significant_results[[paste(bacteria, metabolite, sep = "_")]] <- list(
        Bacteria = bacteria,
        Metabolite = metabolite,
        ACME = mediation_result$d.avg,  # 平均中介效应
        ACME_p = summary(mediation_result)$d.avg.p,  # 中介效应的 p 值
        ADE = mediation_result$z.avg,  # 平均直接效应
        ADE_p = summary(mediation_result)$z.avg.p,  # 直接效应的 p 值
        Mediated = mediation_result$n0,  # 总效应
        Mediated_p = mediation_result$n0.p  ,# 总效应的 p 值
        mediator = summary(mediator_model),
        a = mediator_summary$coefficients[bacteria, "Estimate"],
        a_p = mediator_summary$coefficients[bacteria, "Pr(>|t|)"],
        b = outcome_summary$coefficients[metabolite, "Estimate"],
        b_p = outcome_summary$coefficients[metabolite, "Pr(>|z|)"],
        c = outcome_summary$coefficients[bacteria, "Estimate"],
        c_p = outcome_summary$coefficients[bacteria, "Pr(>|z|)"]
      )
    }
  }
}

# 将结果转换为数据框
if (T) {
  df <- do.call(rbind, lapply(names(significant_results), function(name) {
    cbind(name, as.data.frame(t(significant_results[[name]])))
  }))
  df = as.data.frame(df)
  df <- data.frame(lapply(df, function(col) {
    if (is.list(col)) {
      # 将 list 转换为字符（用逗号连接）
      sapply(col, function(x) paste(x, collapse = ", "))
    } else {
      col
    }
  }), stringsAsFactors = FALSE)

  # 保存结果到文件
  write.csv(df, "./metab/mediation_results0219core.csv", row.names = FALSE)
} else {
  print("没有显著的中介效应结果。")
}
df0 = df[grepl('Faeca',df$Bacteria),]
df1 = df[grepl('Bifido',df$Bacteria),]
