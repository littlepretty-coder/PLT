data = read.delim('./data/属水平绝对定量.txt')
data = data[,-c(1,3)]
tmp = aggregate(data[,2:ncol(data)],by = list(data$Genus),sum)
tmp = tmp[tmp$Group.1 %in% names,]
tmp$sum = rowSums(tmp[,2:ncol(tmp)])
tmp0 = tmp[,c(1,ncol(tmp))]
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
tmp0$Module <- 0
for (name in names(mapping)) {
  tmp0$Module[tmp0$Group.1 %in% mapping[[name]]] <- name
}
df_with_percentage <- tmp0 %>%
  group_by(Module) %>%
  mutate(
    Total = sum(sum),
    Percentage = (sum / Total) * 100
  ) %>%
  arrange(Module, desc(sum)) %>%  # 按 Module 分组后按丰度降序排列
  mutate(Rank = row_number())    # 添加排名

# 2. 合并排名以外的物种为 "其他"
df_top3 <- df_with_percentage %>%
  mutate(Group.1 = ifelse(Rank > 3, "Other", Group.1)) %>% # 合并“其他”物种
  group_by(Module, Group.1) %>% # 重新分组
  summarise(
    sum = sum(sum),
    Percentage = sum(Percentage),
    .groups = "drop"
  )
colorb = read.delim('./fig1/1D/bacteria_genus_color.txt')
colorf = read.delim('./fig1/1E/pca_fungi_genus_color.txt')
colnames(colorb) = colnames(colorf)
color = rbind(colorb,colorf)
color

df_top3 <- df_top3 %>%
  left_join(color, by = c("Group.1" = "common_genera")) # 使用物种名称匹配颜色
df_top3[df_top3[, 2] == "Other", 5] <- "gray94"
# 4. 绘制饼图并使用对应颜色
ggplot(df_top3, aes(x = "", y = Percentage, fill = Group.1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Species Percentage (Top 3 + Other) in Each Module",
       x = NULL, y = NULL, fill = "Species") +
  scale_fill_manual(values = setNames(df_top3$color, df_top3$Group.1)) + # 设置自定义颜色
  theme_void() +
  theme(legend.position = "right") +
  facet_wrap(~ Module)
dev.off()
dev.new()
write.table(df_top3,'./figure/suppdata/fig2/bacterial_module_abundance_piedata.csv',sep=',',row.names = F)
