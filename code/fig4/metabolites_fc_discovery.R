library(data.table)
metab = fread('./metab/wilcox.txt')
metab$fc = log2(metab$median_w/metab$median_wo)
metab = metab[order(metab$fc),]
metab0 = metab[c(1:10,(nrow(metab)-9):nrow(metab)),]
metab00 = metab0[,c(1,28,33)]
metab00$fdr = p.adjust(metab00$p_w_wo)
metab00$stars <- ifelse(metab00$fdr < 0.001, "***",
                    ifelse(metab00$fdr < 0.01, "**",
                           ifelse(metab00$fdr < 0.05, "*", "")))

library(ggplot2)
ggplot(metab00, aes(x = reorder(key,fc))) +
  geom_bar(aes(y=fc,fill=key),stat = "identity", fill = "steelblue") +
  geom_text(aes(y = fc, label = stars), size = 10) +  # 添加星号
  theme_minimal() +
  coord_flip()+
  labs(title = "Fold Change with Significance",
       x = "Group", y = "Fold Change (FC)")
info = read_excel('./metab/all.xlsx')
info = info[info$`Metabolite Name` %in% metab00$key,]
names = metab00$key
