
library(microeco)
library(meconetcomp)
library(magrittr)
library(dplyr)
library(tidyr)
library(igraph)

set.seed(43)

b_with_status[is.na(b_with_status)]=0
b_with_status0 = aggregate(b_with_status$dose_g,by = list(b_with_status$ID,b_with_status$drug_day,
                                                          b_with_status$infection_status,b_with_status$type),sum)
colnames(b_with_status0) = c("ID" ,"drug_day" ,"infection_status",'type',"dose_g")
b_with_status0$SampleID = paste0(b_with_status0$ID,'-D',b_with_status0$drug_day)
b_with_status0 <- b_with_status0 %>%
  group_by(SampleID) %>%
  summarise(
    dose_g = sum(dose_g, na.rm = TRUE),

    infection_status =first(infection_status), 
    .groups = 'drop'
  )
meta0 = merge(meta,b_with_status0,by = 'SampleID',all.x = T)
meta0[is.na(meta0)]=0
meta = meta0[,c(1,2,3,4,5,11,12)]

PATIENT_COL <- "ID"  
TIME_COL    <- "days" 


ABX_COL <- "dose_g"    
INF_COL <- "infection_status"      

meta2 <- meta %>%
  select(SampleID, Group, all_of(PATIENT_COL), all_of(TIME_COL),
         all_of(c(ABX_COL, INF_COL))) %>%
  left_join(baseline_df %>% select(SampleID, Assigned_CST_uncertain), by = "SampleID")

KEEP_CST <- c("fungi_2", "uncertain", "bacteria_2")
meta2 <- meta2 %>% filter(Assigned_CST_uncertain %in% KEEP_CST)

table(meta2$Assigned_CST_uncertain)

fungi_1 = c('Ramularia','Trichosporon','Clavispora','Diutina','Malassezia','Geosmithia','Pichia','Neurpspora','Pyricularia')
fungi_2 = c('Botrytis','Rhizophagus','Aspergillus','Fusarium','Suillus','Colletotrichum','Penicillium','Letharia','Talaromyces','Puccinia','Tetrapisispora')
fungi_3 = c('Dissoconium','Saccharomyces','Xylona','Meyerozyma','Nakaseomyces','Kluyveromyces','Lodderomyces')
fungi_4 = c('Candida','Uncinocarpus','Drepanopeziza','Alternaria','Mycotypha','Phycomyces','Meira','unclassified_Saccharomycetales')

bac_1 = c('Sarcina','Akkermansia','Segatella','Alistipes','Enterocloster','Hungatella','Faecalibacterium','Roseburia','Blautia','Mediterraneibacter',
          'Anaerobutyricum','Intestinibacter','Dysgonomonas','Weissella','Ligilactobacillus')
bac_2 = c('Thomasclavelia','Bacteroides','Phocaeicola','Parabacteroides')
bac_3 = c('Pseudomonas','Kluyvera','Stenotrophomonas','Schaalia','Bifidobacterium','Veillonella','Morganella','Megamonas','Staphylococcus','Lactobacillus',
          'Abiotrophia','Haemophilus','Lacticaseibacillus')
bac_4 = c('Mycobacterium','Klebsiella','Yersinia','Escherichia','Raoultella','Citrobacter','Enterobacter','Streptococcus','Enterococcus','Rothia')

names_all <- c(fungi_1,fungi_2,fungi_3,fungi_4,bac_1,bac_2,bac_3,bac_4)

df2 <- df %>%
  filter(Genus %in% names_all) %>%
  select(Species, Genus, Kingdom, all_of(meta2$SampleID))

rownames(df2) <- df2$Species
df2 <- df2[, -(1:3)]
rownames(df2) <- gsub('-', '_', rownames(df2))
patterns <- c("\\[", "\\]", "\\(", "\\)", "\\.", "\\'", "\\/", "\\:", "\\=", "\\+")
rownames(df2) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(df2))

tax <- read.delim('./data/tax.txt')
rownames(tax) <- tax$Species
rownames(tax) <- gsub('-', '_', rownames(tax))
rownames(tax) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(tax))
tax <- tax[rownames(tax) %in% rownames(df2), ]
tax <- tidy_taxonomy(tax)

sample_table <- meta2 %>%
  select(SampleID, Group, all_of(PATIENT_COL), all_of(TIME_COL),
         any_of(c(ABX_COL, INF_COL)),
         Assigned_CST_uncertain)

rownames(sample_table) <- sample_table$SampleID

mt <- microtable$new(sample_table = sample_table,
                     otu_table = df2,
                     tax_table = tax)

build_cst_network <- function(mt_obj, cst_name,
                              cor_p = 0.001, cor_cut = 0.6,
                              min_samples = 30) {
  tmp <- clone(mt_obj)
  tmp$sample_table %<>% subset(Assigned_CST_uncertain == cst_name)
  tmp$tidy_dataset()
  if (nrow(tmp$sample_table) < min_samples) stop(paste("Too few samples for", cst_name))
  
  net <- trans_network$new(
    dataset = tmp,
    cor_method = "spearman",
    use_WGCNA_pearson_spearman = TRUE,
    nThreads = 10
  )
  net$cal_network(COR_p_thres = cor_p, COR_cut = cor_cut)
  net
}

networks <- list(
  fungi_2    = build_cst_network(mt, "fungi_2",    cor_p = 0.001, cor_cut = 0.6, min_samples = 30),
  uncertain  = build_cst_network(mt, "uncertain",  cor_p = 0.001, cor_cut = 0.6, min_samples = 30),
  bacteria_2 = build_cst_network(mt, "bacteria_2", cor_p = 0.001, cor_cut = 0.6, min_samples = 30)
)

topo_stats <- meconetcomp::cal_network_attr(networks)
topo_df <- topo_stats %>% as.data.frame() %>% tibble::rownames_to_column("CST")


cat("\nCST-level network topology (no bootstrap):\n")
print(topo_df[, c("CST","Vertex","Edge","Average_degree","Density","Modularity","Centralization")])
topo_df = t(topo_df)
colnames(topo_df) = topo_df[1,]
topo_df = topo_df[-1,]
topo_df = rownames_to_column(as.data.frame(topo_df),var = 'CST')

cst_metric_map <- topo_df %>%
  select(CST, Density, Modularity, Centralization)

make_time_order <- function(x) {
  if (all(grepl("^T\\d+$", x))) return(as.integer(sub("^T", "", x)))
  nx <- suppressWarnings(as.numeric(x))
  if (sum(!is.na(nx)) > 0) return(nx)
  as.integer(factor(x, levels = unique(x)))
}

interval_df <- sample_table %>%
  mutate(time_order = make_time_order(.data[[TIME_COL]])) %>%
  arrange(.data[[PATIENT_COL]], time_order) %>%
  group_by(.data[[PATIENT_COL]]) %>%
  mutate(
    CST_t  = Assigned_CST_uncertain,
    CST_t1 = lead(Assigned_CST_uncertain),
  
    abx_interval = lead(.data[[ABX_COL]]),
    inf_interval = lead(.data[[INF_COL]])
  ) %>%
  filter(!is.na(CST_t1)) %>%
  mutate(

    transition_any = as.integer(CST_t != CST_t1),
    trans_f2_u = as.integer(CST_t == "fungi_2"   & CST_t1 == "uncertain"),
    trans_u_b2 = as.integer(CST_t == "uncertain" & CST_t1 == "bacteria_2"),
    
    from_f2 = as.integer(CST_t == "fungi_2"),
    from_u  = as.integer(CST_t == "uncertain")
  ) %>%
  ungroup()
interval_df <- interval_df %>%
  left_join(cst_metric_map, by = c("CST_t" = "CST")) %>%
  rename(
    density_CSTt = Density,
    modularity_CSTt = Modularity,
    centralization_CSTt = Centralization
  )


p_f2u <- with(interval_df, sum(trans_f2_u, na.rm = TRUE) / sum(from_f2, na.rm = TRUE))
p_ub2 <- with(interval_df, sum(trans_u_b2, na.rm = TRUE) / sum(from_u,  na.rm = TRUE))

cat("\nTransition probability (no bootstrap):\n")
cat("P(fungi_2 -> uncertain) =", p_f2u, "\n")
cat("P(uncertain -> bacteria_2) =", p_ub2, "\n")

df_f2 <- interval_df %>% filter(CST_t == "fungi_2")
interval_df$density_CSTt <- as.numeric(interval_df$density_CSTt)

if (requireNamespace("lme4", quietly = TRUE)) {
  formula_str <- paste("trans_f2_u ~ density_CSTt + abx_interval + inf_interval + (1 |", PATIENT_COL, ")")
  
  fit_f2 <- lme4::glmer(
    as.formula(formula_str),
    data = interval_df,
    family = binomial()
  )
  cat("\nGLMM: fungi_2 -> uncertain\n")
  print(summary(fit_f2))
} else {
  fit_f2 <- glm(
    trans_f2_u ~ density_CSTt + abx_interval + inf_interval,
    data = df_f2,
    family = binomial()
  )
  cat("\nGLM: fungi_2 -> uncertain\n")
  print(summary(fit_f2))
}

df_u <- interval_df %>% filter(CST_t == "uncertain")

if (requireNamespace("lme4", quietly = TRUE)) {
  fit_u <- lme4::glmer(
    trans_u_b2 ~ density_CSTt + abx_interval + inf_interval + (1 | .data[[PATIENT_COL]]),
    data = df_u,
    family = binomial()
  )
  cat("\nGLMM: uncertain -> bacteria_2\n")
  print(summary(fit_u))
} else {
  fit_u <- glm(
    trans_u_b2 ~ density_CSTt + abx_interval + inf_interval,
    data = df_u,
    family = binomial()
  )
  cat("\nGLM: uncertain -> bacteria_2\n")
  print(summary(fit_u))
}

if (requireNamespace("lme4", quietly = TRUE)) {
  fit_f2_int <- lme4::glmer(
    trans_f2_u ~ density_CSTt * abx_interval + inf_interval + (1 | .data[[PATIENT_COL]]),
    data = df_f2,
    family = binomial()
  )
  cat("\nGLMM with interaction: fungi_2 -> uncertain\n")
  print(summary(fit_f2_int))
}


