
clean = clean[rowSums(clean)>0,]
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

  sample_name <- colnames(clean_tax)[i]

  merged_df[[sample_name]] <- merged_df[[paste0(sample_name, ".x")]] / merged_df[[paste0(sample_name, ".y")]]
}
selected_columns <- colnames(merged_df)[!grepl("\\.x$|\\.y$", colnames(merged_df))]

filtered_df <- merged_df[, selected_columns]



