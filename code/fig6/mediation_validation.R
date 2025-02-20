# 加载必要的 R 包
if (!require("mediation")) install.packages("mediation")
if (!require("dplyr")) install.packages("dplyr")
library(mediation)
library(dplyr)
library(readxl)
names = merged_data$Species
# 将分组信息转换为因子
merged_data = read.delim('./validation/种水平绝对定量.txt',check.names = F)
merged_data = merged_data[,-1]
merged_data = merged_data[merged_data$Species %in% names,]
group = read.delim('./validation/metadata.txt')
merged_data = as.data.frame(t(merged_data))
colnames(merged_data) = merged_data[1,]
merged_data = merged_data[-1,]

metab = read_excel('./validation/all.xlsx')
colnames(metab) =  gsub("-(\\d+)$", "-D\\1", colnames(metab))
merged_data = merged_data[rownames(merged_data) %in% colnames(metab)[11:84],]

metab = metab[,-c(2:6)]
metab = as.data.frame(t(metab))
colnames(metab) = metab[1,]
metab = metab[-1,]
tmp0 = merge(metab,merged_data,by= 'row.names')
tmp0 = merge(tmp0,group,by.x= 'Row.names',by.y = 'SampleID')

#tmp0$Group <- as.factor(ifelse(tmp0$Group == "W/ Infection", 1, 0))
colnames(tmp0)[1] = 'SampleID'
colnames(tmp0)[ncol(tmp0)] = 'Group'
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
bacteria_names = colnames(merged_data)[380:1641]
#metabolite_names = colnames(merged_data)[2:403]
metabolites <- c("Decanoylcarnitine", "γ-Linolenic acid2", "α-Linolenic acid1", 
                 "Octanoylcarnitine", "Taurodeoxycholic acid", "Taurochenodeoxycholic acid", 
                 "Propionylcarnitine", "Taurohyodeoxycholic acid", "Isovalerylcarnitine", 
                 "Valerylcarnitine", "Guanine", "5-Aminopentanoic acid", "Glucosamine", 
                 "β-Muricholic acid", "P-cresyl sulfate", "3-(4-Hydroxyphenyl)propionic acid", 
                 "Vanillic acid", "Ketoleucine", "7-Dehydrocholic acid", "trans-Ferulic acid")
#metab = metab[-c(66:74),]
metab = metab[metabolites]
metabolites_name <- gsub("[^a-zA-Z0-9_]", "", metabolites)
metabolites_name<- sapply(metabolites_name, function(name) {
  if (grepl("^[0-9]", name)) { # 检查是否以数字开头
    paste0("X", name) # 添加前缀 "X"
  } else {
    name # 保持原列名
  }
})
for (bacteria in bacteria_names) {
  #bacteria = 'Streptococcus_sp_FDAARGOS_521'
  for (metabolite in metabolites_name) {
    #metabolite = "Decanoylcarnitine"
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
    a = summary(mediation_result)

    # 检查 p 值是否显著（中介效应的 p 值）
    if (T) {
      significant_results[[paste(bacteria, metabolite, sep = "_")]] <- list(
        Bacteria = bacteria,
        Metabolite = metabolite,
        ACME = mediation_result$d.avg,  # 平均中介效应
        ACME_p = mediation_result$d.avg.p,  # 中介效应的 p 值
        ADE = mediation_result$z.avg,  # 平均直接效应
        ADE_p = mediation_result$z.avg.p  # 直接效应的 p 值
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
  write.csv(df, "./validation/mediation_results_all.csv", row.names = FALSE)
} else {
  print("没有显著的中介效应结果。")
}
df0 = df[grepl('Faeca',df$Bacteria),]
df1 = df[grepl('Bifido',df$Bacteria),]
