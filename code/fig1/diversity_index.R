#test
library(optparse)
library(vegan)
library(ape)
library(picante)



OTUs <- read.delim('./data/种水平绝对定量.txt',sep='\t',check.names = F)

groups = read.delim("./data/metadata.txt", sep='\t')

bacteria = OTUs[OTUs$Kingdom == 'Bacteria',]
fungi = OTUs[OTUs$Kingdom == 'Eukaryota',]
archaea = OTUs[OTUs$Kingdom == 'Archaea',]



genus=bacteria
rownames(genus) = genus$Species
genus = genus[,-c(1,2)]

#genus <- genus[, match(groups$SampleID, names(genus))]
df=genus
#t = ncol(df)
#df = df[,-t]
df<-sweep(df,MARGIN = 2,colSums(df),'/')
Shannon<-diversity(df, index = "shannon", MARGIN = 2)
Simpson<-diversity(df, index = "simpson", MARGIN = 2)
Richness <- specnumber(df,MARGIN = 2) #spe.rich == sobs
index<-as.data.frame(cbind(Shannon,Simpson,Richness))
tdf<-t(df)
tdf<-ceiling(as.data.frame(t(df)))
obs_chao_ace<-t(estimateR(tdf))
obs_chao_ace<-obs_chao_ace[rownames(index),]
index$Chao1<-obs_chao_ace[,2]
index$Ace<-obs_chao_ace[,4]
index$Sobs<-obs_chao_ace[,1]
write.table(cbind(sample=c(rownames(index)),index),paste0('./fig1/1C/bacteria.diversity.index00.txt'),  row.names = F,sep = '\t',quote = F)
