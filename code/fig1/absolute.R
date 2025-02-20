options (warn = -1)#test
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)
groups = read.delim("./data/metadata.txt", sep='\t')
groups = groups[(groups$Group == 'W/ Infection'),]


OTUs <- read.delim('./data/属水平绝对定量.txt',sep='\t',check.names = F)

bacteria = OTUs[OTUs$Kingdom == 'Bacteria',]
fungi = OTUs[OTUs$Kingdom == 'Eukaryota',]
archaea = OTUs[OTUs$Kingdom == 'Archaea',]


# 第一个向量
fungi_absolute <- c("Trichosporon", "Candida", "Saccharomyces", "Meyerozyma", 
                    "Nakaseomyces", "Ramularia", "Neurospora", "Aspergillus", 
                    "Clavispora", "Pyricularia", "Penicillium", "Geosmithia", 
                    "Lodderomyces", "Colletotrichum", "Fusarium")

fungi_relative <- c("Trichosporon", "Candida", "Neurospora", "Saccharomyces", 
                    "Ramularia", "Aspergillus", "Pyricularia", "Penicillium", 
                    "Nakaseomyces", "Botrytis", "Geosmithia", "Colletotrichum", 
                    "Fusarium", "Rhizophagus", "Malassezia")
fungi_union <- union(fungi_absolute, fungi_relative)

bacteria_absolute <- c("Klebsiella", "Veillonella", "Bacteroides", "Bifidobacterium", 
                       "Pseudomonas", "Mediterraneibacter", "Streptococcus", "Blautia", 
                       "Phocaeicola", "Enterocloster", "Enterococcus", "Rothia", 
                       "Yersinia", "Thomasclavelia", "Escherichia")

bacteria_relative <- c("Klebsiella", "Bifidobacterium", "Veillonella", "Streptococcus", 
                       "Bacteroides", "Enterococcus", "Yersinia", "Rothia", 
                       "Escherichia", "Phocaeicola", "Pseudomonas", "Enterobacter", 
                       "Schaalia", "Faecalibacterium", "Mediterraneibacter")
bacteria_union <- union(bacteria_absolute, bacteria_relative)



df = fungi
df = df[,-c(1,3)]
df = aggregate(df[,c(2:ncol(df))],by = list(df$Genus),sum)
colnames(df)[1] = 'Genus'
rownames(df) = df$Genus
df = df[,-1]
df = df[,colnames(df)%in%groups$SampleID]
data = df
df1 = data[rownames(data)%in%fungi_union,]
#df1<-sweep(df1,MARGIN = 2,colSums(df1),'/')


df1$name = rownames(df1)
df1 = df1[,c(ncol(df1),1:(ncol(df1)-1))]
df1 = df1[,-ncol(df1)]
df = df1
#df[,-1]<-sweep(df[,-1],MARGIN = 2,colSums(df[,-1]),'/')
colnames(df)[1]<-'Tax'
df1<-melt(df,id.vars = 'Tax')
colnames(df1)<-c('variable','samples','value')
df1$variable<-factor(df1$variable,levels = rev(df$Tax))
groups<-groups[,c(1,5)]
colnames(groups)<-c('samples','group')
groups <- groups[trimws(groups$group) != "", ]
df2<-merge(df1,groups,by='samples')
df3<-df2

df4<-as.data.frame(df3 %>% group_by(group,variable) %>% dplyr::summarise(mymean=mean(value),mysd=sd(value)))
df5<-as.data.frame(df4 %>%mutate(ybegin=mymean-mysd,yend=mymean+mysd))
df5$group<-factor(df5$group,levels=groups$group[!duplicated(groups$group)])
df6<-as.data.frame(df5 %>% dplyr::arrange(group,desc(variable)))
df7<-as.data.frame(df6 %>% group_by(group) %>%
                     mutate(ybegin2=cumsum(mymean)-mysd,yend2=cumsum(mymean)+mysd))


library(ggplot2)
library(ggalluvial)
color = read.delim('./fig1/1E/pca_fungi_genus_color.txt',sep='\t')
colnames(color) = c('common_genera','color')
#groups = merge(groups,color,by.x = 'Dominant_Genus',by.y = 'common_genera')
color_mapping <- setNames(color$color, color$common_genera)
# colset<-df[1:18,1,drop=F]
# colset$color = mycol22[1:18]
# colset$Tax = rownames(colset)
# colnames(colset)[1:2]<-c('Tax','Color')
# 
# colset$Tax<-factor(colset$Tax,levels = unique(colset$Tax))
# colset<-colset%>%arrange(Tax)
# mycol<- colset$Color
# names(mycol)<- unique(colset$Tax)
# df7<-merge(df7,colset,by.x='variable',by.y='Tax')
#df7$variable<-factor(df7$variable,colset$Tax)

df7$group <- factor(df7$group,levels = c('T0','T1','T2','T3'))
df7$variable<-factor(df7$variable,levels = rev(rownames(df)))
ggplot(df7, aes(x=group,stratum=variable,
                alluvium=variable,fill=variable,label=variable,y=mymean)) + 
  scale_fill_manual(values = color_mapping)+
  geom_stratum(aes(fill=variable),color='black',width = 0.5,size=0.2)+####width
  geom_flow(aes(fill=variable)) +
  #ylim(0,8e+9)+
  #geom_errorbar(aes(ymax=ybegin2 , ymin= yend2 ),width=.2) +
  #geom_text(stat='stratum',size=2)+
  #facet_grid(.~group,scale="free")+
  guides(fill= guide_legend(reverse = TRUE))+
  #scale_y_continuous(limits = c(0,1.2),breaks = c(0,0.25,0.5,0.75,1)) +
  theme(panel.background = element_rect(fill='white', colour='white'), 
        panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        axis.text.x  = element_text(size=24, angle = 0,colour="black", face = "bold",hjust = 0.5),  
        axis.title.x = element_text(vjust=0.1, face = "bold",size=24),
        axis.text.y = element_text(size=24, colour="black"),
        axis.title.y = element_text(vjust=0.2, size = 24, face = "bold"))+
  labs(  y = "Abusolute abundance")

#ggsave(paste0('./1D/绝对定量relative/bacteria.barplot.未感染组.ab_re.pdf'),width = 8,height = 8)


