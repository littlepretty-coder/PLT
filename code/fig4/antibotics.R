library(readxl)
data = read_excel('./review/antibiotic/processed_data(抗生素_用法分割_提取数字_顺延重复)-TJ189.xlsx')
data0 = data[,c(3,8,13)]
data0$使用剂量 = gsub('B10万单位','10mg',data0$使用剂量)
data0$使用剂量 = gsub('B15万单位','15mg',data0$使用剂量)
data0$使用剂量 = gsub('B75万单位','75mg',data0$使用剂量)
data0$使用剂量 = gsub('50万单位','50mg',data0$使用剂量)
data0$num <- as.numeric(gsub("[^0-9.]", "", data0$使用剂量))  # 移除非数字字符 [6,8](@ref)

# 步骤2: 提取单位部分
data0$unit <- gsub("[0-9.]", "", data0$使用剂量)  # 移除数字，保留单位字符 [9](@ref)

# 步骤3: 按单位转换数值
data0$converted_value <- ifelse(
  data0$unit == "mg", 
  data0$num / 1000,   # mg → g: 除以1000
  data0$num            # g: 直接保留
)

tmp = aggregate(data0$converted_value,by = list(data0$ID,data0$顺延重复数字),sum)
library(dplyr)
result <- tmp %>%
  group_by(Group.1) %>%
  summarise(mean = mean(x, na.rm = TRUE)) %>%
  ungroup()
group = read.delim('./data/metadata1.txt')
group = group[,c(2,4)]
group = unique(group)
info = merge(result,group,by.x = 'Group.1',by.y = 'ID')
info$factor = as.factor(info$Group)
#wilcox.test(info$x,info$factor)
# info$factor <- as.factor(info$factor)
# # 检查正态性和方差齐性（ANOVA前提假设）
# shapiro.test(info$x)              # 正态性检验（若p>0.05满足）
# bartlett.test(x ~ factor, data = info) 
# anova(x ~ factor, data = info)
wilcox_result = wilcox.test(mean ~ factor, data = info)
library(broom) 
result_df <- tidy(wilcox_result)
print(result_df)
library(writexl)
write_xlsx(result_df, "./Groups_wilcoxon_test_results.xlsx")
library(vegan)
adonis2_result = adonis2(info$x ~ info$factor, data = info, method = "euclidean")
result_df <- tidy(adonis2_result)
print(result_df)
write_xlsx(result_df, "./Groups_adonis2_test_results.xlsx")
library(ggplot2)
library(ggpubr)
info$factor = factor(info$factor,levels = c('W/O Infection','W/ Infection'))
ggplot(info, aes(x = factor, y = mean)) +
  geom_boxplot(
    width = 0.5, 
    outlier.shape = NA,  # 隐藏异常值
    outlier.colour = NA, # 双重保险
    color = "black",     # 箱线边框颜色
    fill = c("grey80", "#B7B7EB")  # 手动指定填充色
  ) +
  geom_jitter(
    width = 0.1, 
    size = 2,
    alpha = 0.6          # 散点透明度（避免重叠遮挡）
  ) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",
    label.x = 1.5,
    size = 5
  ) +
  labs(
    title = "Compairsions of antibiotic dosage between Groups", 
    x = "Group", 
    y = "Overall antibiotic dosage(g)"
  ) +
  theme_classic() +      # 经典主题（自带边框 + 无网格）
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),  # 标题居中
    axis.line = element_line(color = "black", linewidth = 0.5),  # 加粗坐标轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # 添加外边框
  )
#------------------------------------------------------------------------------
library(dplyr)
data = read_excel('./review/antibiotic/processed_data(抗生素_用法分割_提取数字_顺延重复)-TJ189.xlsx')
data$使用剂量 = gsub('B10万单位','10mg',data$使用剂量)
data$使用剂量 = gsub('B15万单位','15mg',data$使用剂量)
data$使用剂量 = gsub('B75万单位','75mg',data$使用剂量)
data$使用剂量 = gsub('50万单位','50mg',data$使用剂量)
dat = read_excel('./review/antibiotic/小儿肝患者感染时间点信息-129+20例（天津+仁济）.xlsx')
time = dat[,c(2,6)]

a = time
colnames(a) = c('ID','infection_day')
a = a[a$infection_day != '无',]
a$infection_day = as.numeric(a$infection_day)
b = data[,c(3,13,8)]
colnames(b) = c('ID','drug_day','Dose')
b$num <- as.numeric(gsub("[^0-9.]", "", b$Dose))  # 移除非数字字符 [6,8](@ref)

# 步骤2: 提取单位部分
b$unit <- gsub("[0-9.]", "", b$Dose)  # 移除数字，保留单位字符 [9](@ref)

# 步骤3: 按单位转换数值
b$dose_g <- ifelse(
  b$unit == "mg", 
  b$num / 1000,   # mg → g: 除以1000
  b$num            # g: 直接保留
)
b = b[,c(1,2,6)]

# 合并数据框a和b，关联样本ID
b_with_status <- b %>%
  left_join(a, by = "ID") %>%        # 按样本ID合并感染时间
  mutate(
    infection_status = ifelse(drug_day < infection_day, 0, 1)  # 判断用药时间
  ) %>%
  select(-infection_day)  

b_with_status = na.omit(b_with_status)
b_with_status$infection_status = as.factor(b_with_status$infection_status)
wilcox_result = wilcox.test(dose_g ~ infection_status, data = b_with_status)

ggplot(b_with_status, aes(x = infection_status, y = dose_g)) +
  geom_boxplot(
    width = 0.5, 
    outlier.shape = NA,  # 隐藏异常值
    outlier.colour = NA, # 双重保险
    color = "black",     # 箱线边框颜色
    fill = c("grey80", "#B7B7EB")  # 手动指定填充色
  ) +
  geom_jitter(
    width = 0.1, 
    size = 2,
    alpha = 0.6          # 散点透明度（避免重叠遮挡）
  ) +
  stat_compare_means(
    method = "anova", 
    label = "p.format",
    label.x = 1.5,
    size = 5
  ) +
  labs(
    title = "Compairsions of antibiotic dosage between Groups", 
    x = "Group", 
    y = "Overall antibiotic dosage(g)"
  ) +
  theme_classic() +      # 经典主题（自带边框 + 无网格）
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),  # 标题居中
    axis.line = element_line(color = "black", linewidth = 0.5),  # 加粗坐标轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # 添加外边框
  )



library(tidyr)
b_with_status0 = aggregate(b_with_status$dose_g,by = list(b_with_status$ID,b_with_status$drug_day,
                                                          b_with_status$infection_status),sum)
colnames(b_with_status0) = c("ID" ,"drug_day" ,"infection_status","dose_g")

mean_dose <- b_with_status0 %>%
  group_by(ID, infection_status) %>%
  summarise(mean_dose = mean(dose_g, na.rm = TRUE)) %>%
  ungroup()
wide_dose <- mean_dose %>%
  pivot_wider(
    names_from = infection_status,
    values_from = mean_dose,
    names_prefix = "mean_dose_"
  )
wilcox_result = wilcox.test(wide_dose$mean_dose_0, wide_dose$mean_dose_1, paired = T)
library(broom) 
result_df <- tidy(wilcox_result)
print(result_df)
library(writexl)
write_xlsx(result_df, "./B&A_wilcoxon_test_results.xlsx")
# model <- glm(mean_dose ~ infection_status, data = mean_dose)
# anova_result = anova(model)
# result_df <- tidy(anova_result)
# write_xlsx(result_df, "./B&A_anova_test_results.xlsx")
library(ggplot2)
library(ggpubr)

ggplot(wide_dose, aes(x = factor(1), y = mean_dose_0, fill = "Before Infection(W/ Infection)")) +
  geom_boxplot(width = 0.5) +
  geom_boxplot(aes(x = factor(2), y = mean_dose_1, fill = "After Infection(W/ Infection)"), width = 0.5) +
  geom_segment(
    aes(x = 1, xend = 2, y = mean_dose_0, yend = mean_dose_1),
    color = "gray", alpha = 0.6
  ) +
  geom_point(aes(x = 1, y = mean_dose_0), color = "black") +
  geom_point(aes(x = 2, y = mean_dose_1), color = "black") +
  stat_compare_means(
    paired = TRUE,
    label = "p.format",
    label.x = 1.5,
    size = 5
  ) +
  scale_fill_manual(values = c("#9BBBE1", "grey50")) +
  labs(x = "", y = "Average antibiotic(g)", fill = "Infection Status") +
  theme_classic() +      # 经典主题（自带边框 + 无网格）
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),  # 标题居中
    axis.line = element_line(color = "black", linewidth = 0.5),  # 加粗坐标轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # 添加外边框
  )
#-------------------------------------------------------------------------------------

library(dplyr)
data = read_excel('./review/antibiotic/processed_data(抗生素_用法分割_提取数字_顺延重复)-TJ189.xlsx')
data$使用剂量 = gsub('B10万单位','10mg',data$使用剂量)
data$使用剂量 = gsub('B15万单位','15mg',data$使用剂量)
data$使用剂量 = gsub('B75万单位','75mg',data$使用剂量)
data$使用剂量 = gsub('50万单位','50mg',data$使用剂量)
dat = read_excel('./review/antibiotic/小儿肝患者感染时间点信息-129+20例（天津+仁济）.xlsx')
time = dat[,c(2,6)]

a = time
colnames(a) = c('ID','infection_day')
a = a[a$infection_day != '无',]
a$infection_day = as.numeric(a$infection_day)
b = data[,c(3,13,8,7)]
colnames(b) = c('ID','drug_day','Dose','type')
b$num <- as.numeric(gsub("[^0-9.]", "", b$Dose))  # 移除非数字字符 [6,8](@ref)

# 步骤2: 提取单位部分
b$unit <- gsub("[0-9.]", "", b$Dose)  # 移除数字，保留单位字符 [9](@ref)

# 步骤3: 按单位转换数值
b$dose_g <- ifelse(
  b$unit == "mg", 
  b$num / 1000,   # mg → g: 除以1000
  b$num            # g: 直接保留
)
b = b[,c(1,2,4,7)]

# 合并数据框a和b，关联样本ID
b_with_status <- b %>%
  left_join(a, by = "ID") %>%        # 按样本ID合并感染时间
  mutate(
    infection_status = ifelse(drug_day < infection_day, 0, 1)  # 判断用药时间
  ) %>%
  select(-infection_day)  

b_with_status = na.omit(b_with_status)
b_with_status0 = aggregate(b_with_status$dose_g,by = list(b_with_status$ID,b_with_status$drug_day,
                                                          b_with_status$infection_status,b_with_status$type),sum)
colnames(b_with_status0) = c("ID" ,"drug_day" ,"infection_status",'type',"dose_g")

mean_dose <- b_with_status0 %>%
  group_by(ID,type, infection_status) %>%
  summarise(mean_dose = mean(dose_g, na.rm = TRUE),
            .groups = "drop") %>%
  ungroup()
wide_dose <- mean_dose %>%
  pivot_wider(
    names_from = c(type, infection_status),
    values_from = mean_dose,
    names_prefix = "mean_"
  )
# 提取抗细菌药在感染前后的剂量
antibacterial_data <- mean_dose %>%
  filter(type == "抗细菌") %>%
  select(ID, infection_status, mean_dose) %>%
  pivot_wider(
    names_from = infection_status,
    values_from = mean_dose,
    names_prefix = "mean_"
  )
antibacterial_data[is.na(antibacterial_data)] = 0
# Wilcoxon配对检验（非参数）
wilcox_antibacterial <- wilcox.test(
  antibacterial_data$mean_0,
  antibacterial_data$mean_1,
  paired = TRUE
)

# 提取抗真菌药在感染前后的剂量
antifungal_data <- mean_dose %>%
  filter(type == "抗真菌") %>%
  select(ID, infection_status, mean_dose) %>%
  pivot_wider(
    names_from = infection_status,
    values_from = mean_dose,
    names_prefix = "mean_"
  )
antifungal_data[is.na(antifungal_data)] = 0
# Wilcoxon配对检验（非参数）
wilcox_antifungal <- wilcox.test(
  antifungal_data$mean_0,
  antifungal_data$mean_1,
  paired = TRUE
)

#-----------------------------------------------------------------------------------




data = read_excel('./processed_data(抗生素_用法分割_提取数字_顺延重复)-TJ189.xlsx')
data$使用剂量 = gsub('B10万单位','10mg',data$使用剂量)
data$使用剂量 = gsub('B15万单位','15mg',data$使用剂量)
data$使用剂量 = gsub('B75万单位','75mg',data$使用剂量)
data$使用剂量 = gsub('50万单位','50mg',data$使用剂量)

b = data[,c(3,13,8,7)]
colnames(b) = c('ID','drug_day','Dose','type')
b$num <- as.numeric(gsub("[^0-9.]", "", b$Dose))  # 移除非数字字符 [6,8](@ref)

# 步骤2: 提取单位部分
b$unit <- gsub("[0-9.]", "", b$Dose)  # 移除数字，保留单位字符 [9](@ref)

# 步骤3: 按单位转换数值
b$dose_g <- ifelse(
  b$unit == "mg", 
  b$num / 1000,   # mg → g: 除以1000
  b$num            # g: 直接保留
)
b = b[,c(1,2,4,7)]
b = b[b$ID %in% group$ID,]
group = read.delim('../data/metadata1.txt')
group = group[,c(2,4)]
group = unique(group)
info = merge(b,group,by = 'ID')
info0 = aggregate(info$dose_g,by = list(info$ID,info$drug_day,info$type,info$Group),sum)
colnames(info0) = c('ID','drug_day','type','Group','Dose')
aaa = info0[info0$type == '抗真菌',]
aaa$factor = as.factor(aaa$Group)
wilcox.test(Dose ~ Group, data = aaa)
