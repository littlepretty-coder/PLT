clean = read.delim('./data/data_clean.txt',sep='\t',check.names = F)
clean = clean[rowSums(clean)>0,]

data = read.delim('./data/decontam/四界物种信息-20cases-去人类.txt',sep='\t',check.names = F)
filtered_df <- data[grepl("\\|s_", data$Taxonomy), ]
df_extracted <- filtered_df %>%
  mutate(
    Kingdom = str_extract(Taxonomy, "(?<=k__)[^|]*"),
    Phylum = str_extract(Taxonomy, "(?<=p__)[^|]*"),
    Class = str_extract(Taxonomy, "(?<=c__)[^|]*"),
    Order = str_extract(Taxonomy, "(?<=o__)[^|]*"),
    Family = str_extract(Taxonomy, "(?<=f__)[^|]*"),
    Genus = str_extract(Taxonomy, "(?<=g__)[^|]*"),
    Species = str_extract(Taxonomy, "(?<=s__)[^|]*")
  )
tax = df_extracted[,c(384:390)]
write.table(tax,'./data/tax.txt',sep='\t',row.names = F)
count = df_extracted[(df_extracted$Genus %in% c('Salinibacter','Haloarcula','Trichoderma')),] 
count = count[,c(389,2:383)]
count0 = aggregate(count[,c(2:383)],by = list(count$Genus),sum)
colnames(count0)[1] = 'G'

count1 = read_excel('./data/标准株计算比值及种属的绝对丰度(10000000,50000,20000)/标准株计算比值及种属的绝对丰度(10000000,50000,20000)/标准株计算比值及种属的绝对丰度(10000000,50000,20000).xlsx')
count1 = count1[1:3,]
count1$G = gsub('g__','',count1$G)

count = merge(count0,count1,by= 'G')


sort_by_patient_and_day <- function(df) {
  split_data <- do.call(rbind, strsplit(colnames(df), "-"))
  split_data <- data.frame(patient = split_data[, 1], day = split_data[, 2], stringsAsFactors = FALSE)
  split_data$day <- ifelse(split_data$day == "SQ", 0, as.numeric(sub("D", "", split_data$day)))
  sorted_order <- order(split_data$patient, split_data$day)
  df_sorted <- df[, sorted_order]
  return(df_sorted)
}



count <- sort_by_patient_and_day(count)
count = count[c(2,1,3),]
rownames(count) = count$G
count = count[,-1]

coefficients <- c(10000000, 50000, 20000)
coff <- sweep(count, 1, coefficients, "/")
rownames(coff) = c('Bacteria','Archaea','Eukaryota')


clean_tax = merge(tax[,c(1,7)],clean,by.x = 'Species',by.y = 'row.names')
clean_tax = clean_tax[clean_tax$Kingdom != 'Viruses',]


merged_df <- merge(clean_tax, coff, by.x = "Kingdom",by.y = 'row.names')
for (i in 3:(ncol(clean_tax))) {
  #i=3
  # 样本列名（对应样本丰度和标准化系数）
  sample_name <- colnames(clean_tax)[i]
  
  # 对于每个样本列，除以标准化系数
  merged_df[[sample_name]] <- merged_df[[paste0(sample_name, ".x")]] / merged_df[[paste0(sample_name, ".y")]]
}
selected_columns <- colnames(merged_df)[!grepl("\\.x$|\\.y$", colnames(merged_df))]

# 创建新的数据框，只保留列名中不含 ".x" 和 ".y" 的列
filtered_df <- merged_df[, selected_columns]
write.table(filtered_df,'./data/种水平绝对定量.txt',sep='\t',row.names = F)
write.table(clean_tax,'./data/种水平相对定量.txt',sep='\t',row.names = F)



#--------------------------------------------------------------------------------------------------

species_abusolute = read.delim('./data/种水平绝对定量.txt',sep='\t',check.names = F)
species_relative = read.delim('./data/种水平相对定量.txt',sep='\t',check.names = F)
tax = read.delim('./data/tax.txt',sep='\t')

genus_ab = merge(tax[,c(6,7)],species_abusolute,by = 'Species')
genus_re = merge(tax[,c(6,7)],species_relative,by = 'Species')
write.table(genus_ab,'./data/属水平绝对定量.txt',sep='\t',row.names = F)
write.table(genus_re,'./data/属水平相对定量.txt',sep='\t',row.names = F)

#---------------------------------分组整理------------------------------------------------------------------
metadata = as.data.frame(colnames(species_abusolute)[3:723])
colnames(metadata) = 'SampleID'
split_data <- do.call(rbind, strsplit(metadata$SampleID, "-"))
metadata$ID <- split_data[, 1]
metadata$days <- ifelse(split_data[, 2] == "SQ", -1, as.numeric(sub("D", "", split_data[, 2])))

metadata1 = read_excel("./109/109data/四界种信息-绝对定量.xlsx", sheet = "Sample")
metadata = merge(metadata,metadata1,by = 'SampleID',all=T)
metadata <- sort_by_patient_and_day(metadata, "SampleID")
write.table(metadata,'./data/metadata.txt',sep='\t',row.names = F)


metadata = read.delim('./data/metadata.txt',sep='\t')
metadata$group1 <- ifelse(metadata$days == -1, "移植前",
                          ifelse(metadata$days > 0 , "移植后",NA))


metadata$group2 <- ifelse(metadata$days == -1 & metadata$Group == 'W/ Infection', "感染移植前",
                          ifelse(metadata$days > 0 & metadata$Group == 'W/ Infection', "感染移植后",
                                 ifelse(metadata$days > 0 & metadata$Group == 'W/O Infection', "未感染移植后",
                                        ifelse(metadata$days == -1 & metadata$Group == 'W/O Infection', "未感染移植前",
                                               NA))))
write.table(metadata,'./data/metadata.txt',sep='\t',row.names = F)
