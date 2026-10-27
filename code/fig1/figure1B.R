
library(readxl)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


#display.brewer.all()
blues <- brewer.pal(4, "Blues")
tmp = aggregate(OTUs[,3:ncol(OTUs)],by = list(OTUs$Kingdom),sum)
rownames(tmp) = tmp[,1]
tmp = tmp[,-1]
tmp = as.data.frame(t(tmp))
tmp = merge(tmp,groups,by.x = 'row.names',by.y = 'SampleID')
colnames(tmp)[1] = 'SampleID'



bacteria = tmp[,c(3,5:10)]
bacteria$log10 = log10(bacteria$Bacteria)

fungi = tmp[,c(4,5:10)]
fungi$log10 = log10(fungi$Eukaryota)

archaea = tmp[,c(2,5:10)]
archaea$log10 = log10(archaea$Archaea)

bacteria = bacteria[,-5]
fungi = fungi[,-5]
archaea = archaea[,-5]



p <- ggplot(bacteria, aes(x = days, y=log10,color = Group))
p + geom_point()+xlim(-1,40)+ylim(-1,12)+theme_minimal() +
  #geom_smooth(aes(group = Group),method = 'lm', formula = y ~ x, se = T)+
  #stat_cor(aes(group = Group), method = "pearson")+
  scale_color_manual(values = c('black','black'))  +
  theme(panel.border = element_blank(),                
        #axis.line.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) + 
  theme(panel.grid = element_blank()) +
  xlab("Day") +             
  ylab("log10(total abundance)") +
ggtitle("Bacteria")
# 
p <- ggplot(fungi, aes(x = days, y=log10,color = Group))
p + geom_point()+xlim(-1,40)+ylim(-1,12)+theme_minimal() +
  geom_smooth(aes(group = Group),method = 'lm',  se = T)+
   #stat_compare_means(aes(group = Group), method = "t.test")
  stat_cor(aes(group = Group), method = "pearson")+
  scale_color_manual(values = c('red','blue'))  +
  theme(panel.border = element_blank(),                
        #axis.line.x = element_blank(),
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black")) +  
  theme(panel.grid = element_blank()) +
  xlab("Day") +             
  ylab("log10(total abundance)") +
  ggtitle("Fungi")
# 
# 
p <- ggplot(archaea, aes(x = days, y=log10,color = Group))
p + geom_point()+xlim(-1,40)+ylim(-1,12)+theme_minimal() +
  geom_smooth(aes(group = Group),method = 'lm', formula = y ~ x, se = T)+
  stat_cor(aes(group = Group), method = "pearson")+
  scale_color_manual(values = c('red','blue'))  +
  theme(panel.border = element_blank(),    
        #axis.line.x = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +  
  theme(panel.grid = element_blank()) +
  xlab("Day") +    
  ylab("log10(total abundance)") +
  ggtitle("Archaea")

#-----------------------c("#D7DF23","grey" ,"#557F3A","grey"  ))----------------------------------------------------------------------------------
theme_set(theme_gray(base_family = "SimSun"))
p <- ggplot(bacteria, aes(x = days, y=log10,color = group2))
p + geom_point(size=3)+xlim(-1,40)+ylim(0,12)+theme_minimal() +
  #geom_smooth(data =bacteria ,method = 'lm', formula = y ~ x, se = T)+
  #stat_cor(data =bacteria, method = "pearson")+
  scale_color_manual(values=c("black","black" ,"black","black"  )) +
  #scale_color_manual(values = c("#c72228","#F5867F" ,"#0C4E9B","#6B98C4"  )) +
  #scale_color_manual(values = rainbow(length(unique(bacteria$group1)))) +
  theme(panel.border = element_blank(),              
        #axis.line.x = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +  
  theme(panel.grid = element_blank()) +
  xlab("Day") +    
  ylab("log10(total abundance)") +  
  ggtitle("Bacteria")
ggsave('./fig1/1B/bacteria1.pdf',width = 7,height = 5)
p <- ggplot(fungi, aes(x = days, y=log10,color = group2))
p + geom_point(size=3)+xlim(-1,40)+ylim(0,12)+theme_minimal() +  
  #geom_smooth(data =fungi,method = 'lm',  se = T)+
  #stat_compare_means(aes(group = Group), method = "t.test")
  #stat_cor(data =fungi, method = "pearson")+
  scale_color_manual(values =c("black","black" ,"black","black"  )) +
  #scale_color_manual(values = c('red','blue')) +
  theme(panel.border = element_blank(),              
        #axis.line.x = element_blank(), 
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black")) + 
  theme(panel.grid = element_blank()) +
  xlab("Day") +            
  ylab("log10(total abundance)") +  
  ggtitle("Fungi")
#ggsave('./1B/Fungi1.pdf',width = 7,height = 5)

p <- ggplot(archaea, aes(x = days, y=log10,color = group2))
p + geom_point(size=3)+xlim(-1,40)+ylim(-1,12)+theme_minimal() +  
  #geom_smooth(data = archaea,method = 'lm', formula = y ~ x, se = T)+
  #stat_cor(data = archaea, method = "pearson")+
  scale_color_manual(values =c("black","black" ,"black","black"  )) +
  #scale_color_manual(values = c('red','blue')) +
  theme(panel.border = element_blank(),               
        #axis.line.x = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +  
  theme(panel.grid = element_blank()) +
  xlab("Day") +             
  ylab("log10(total abundance)") +  
  ggtitle("Archaea")

