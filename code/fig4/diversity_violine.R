#test
#rm(list=ls())
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)

mygroup$group<-factor(mygroup$group,levels = unique(mygroup$group))
df2<-merge(df1,mygroup,by='sample')
myindex<-'Richness'
#myindex<-'Richness'
df2$new<-df2[,myindex]

stat.test <- df2 %>%
  wilcox_test(new ~ group) 

# Box plots
stat.test <- stat.test %>% 
  add_xy_position(x = "group", dodge = 0.8,step.increase = 0.5)

stat.test$y.position
#mycol<-c("#EFCE87" ,"#ED8D5A")
mycol<-c("#CFDD97" ,"#EFCE87",'#FFF2AE','#BEBADA','#95c4cc')
stat.test$p.adj.signif<-cut(stat.test$p,breaks =c(0,0.001,0.01,0.05,1),
                            labels = c('***','**','*','ns')  )

ggviolin(df2, x = "group", y = "new", fill = "group",
         add=c('boxplot'),shape=21,
         add.params = list(width=0.1))+
  geom_jitter(position = position_jitter(width = 0.2), size = 0.2, alpha = 1) +
  stat_pvalue_manual(stat.test,   label = "p", 
                     size = 5,
                     hide.ns = F)+ 
  scale_fill_manual(values = mycol)+
  labs(y=myindex,x=NULL)+
  #ylim(-3,300)+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())


#----------------------------------------------------------------------------------------------------------------

