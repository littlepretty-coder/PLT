# =========================
# 0) Libraries
# =========================
library(microeco)
library(meconetcomp)
library(magrittr)
library(dplyr)
library(tidyr)
library(igraph)

set.seed(43)

# =========================
# 1) Load data (你已有的部分，基本不改)
# =========================
baseline_df <- read.delim('./nc返修/figure3/baseline.txt', sep = '\t')

meta <- read.delim('./data/metadata.txt', sep = '\t')

# 你之前构造的 Group（这里保留）
meta$Timepoint0 <- ifelse(meta$Group == "W/ Infection",
                          paste0("W_", meta$Timepoint),
                          ifelse(meta$Group == "W/O Infection",
                                 paste0("WO_", meta$Timepoint),
                                 meta$Timepoint))
b_with_status[is.na(b_with_status)]=0
b_with_status0 = aggregate(b_with_status$dose_g,by = list(b_with_status$ID,b_with_status$drug_day,
                                                          b_with_status$infection_status,b_with_status$type),sum)
colnames(b_with_status0) = c("ID" ,"drug_day" ,"infection_status",'type',"dose_g")
b_with_status0$SampleID = paste0(b_with_status0$ID,'-D',b_with_status0$drug_day)
#b_with_status0 = aggregate(b_with_status0$dose_g,by = list(b_with_status0$SampleID),sum)
b_with_status0 <- b_with_status0 %>%
  group_by(SampleID) %>%
  summarise(
    dose_g = sum(dose_g, na.rm = TRUE),
    # 明确指定要保留的其他列
    #drug_day =first(infection_status),  # 假设CST在同一SampleID下相同
    infection_status =first(infection_status), 
    .groups = 'drop'
  )
meta0 = merge(meta,b_with_status0,by = 'SampleID',all.x = T)
meta0[is.na(meta0)]=0
meta = meta0[,c(1,2,3,4,5,11,12)]

# 你原来只取了两列：SampleID + Group
# ⚠️ 注意：后面 bootstrap 需要 PatientID + 时间点顺序
# 所以这里建议保留更多列（至少 patient + time）
# 下面用占位列名，你需要按实际列名替换
PATIENT_COL <- "ID"   # <<< 改成你真实的病人ID列名
TIME_COL    <- "days"   # <<< 改成你真实的时间点列名（T0/T1/T2/T3或数字）

# 如果你有抗生素/感染变量（可选）
ABX_COL <- "dose_g"     # <<< 没有就留着，后面会自动跳过
INF_COL <- "infection_status"       # <<< 没有就留着，后面会自动跳过

# 合并 CST label
meta2 <- meta %>%
  select(SampleID, Group, all_of(PATIENT_COL), all_of(TIME_COL),
         all_of(c(ABX_COL, INF_COL))) %>%
  left_join(baseline_df %>% select(SampleID, Assigned_CST_uncertain), by = "SampleID")

# 只保留三个可用 CST
KEEP_CST <- c("fungi_2", "uncertain", "bacteria_2")
meta2 <- meta2 %>% filter(Assigned_CST_uncertain %in% KEEP_CST)

table(meta2$Assigned_CST_uncertain)

# 读绝对定量丰度（你已有）
df <- read.delim('./data/属水平绝对定量.txt', sep = '\t', check.names = FALSE)

# 你自己的 genus 过滤列表（你已有）
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

# tax table（你已有）
tax <- read.delim('./data/tax.txt')
rownames(tax) <- tax$Species
rownames(tax) <- gsub('-', '_', rownames(tax))
rownames(tax) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(tax))
tax <- tax[rownames(tax) %in% rownames(df2), ]
tax <- tidy_taxonomy(tax)

# microtable：sample_table 必须 rownames = SampleID
sample_table <- meta2 %>%
  select(SampleID, Group, all_of(PATIENT_COL), all_of(TIME_COL),
         any_of(c(ABX_COL, INF_COL)),
         Assigned_CST_uncertain)

rownames(sample_table) <- sample_table$SampleID

mt <- microtable$new(sample_table = sample_table,
                     otu_table = df2,
                     tax_table = tax)

# =========================
# 4) Build CST networks (full data, no bootstrap)
# =========================
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
write.csv(topo_df, "./nc返修/figure3/all/network_topology_stats_no_bootstrap.csv", row.names = FALSE)

cat("\nCST-level network topology (no bootstrap):\n")
print(topo_df[, c("CST","Vertex","Edge","Average_degree","Density","Modularity","Centralization")])
topo_df = t(topo_df)
colnames(topo_df) = topo_df[1,]
topo_df = topo_df[-1,]
topo_df = rownames_to_column(as.data.frame(topo_df),var = 'CST')
# 为 interval 模型准备：把你要用的网络指标做成映射表
# （建议主分析用 Density；你也可换 Modularity/Centralization）
cst_metric_map <- topo_df %>%
  select(CST, Density, Modularity, Centralization)

# =========================
# 5) Interval table (patient trajectories, full data)
# =========================
make_time_order <- function(x) {
  # 支持 T0/T1/T2/T3
  if (all(grepl("^T\\d+$", x))) return(as.integer(sub("^T", "", x)))
  
  # 支持数字
  nx <- suppressWarnings(as.numeric(x))
  if (sum(!is.na(nx)) > 0) return(nx)
  
  # fallback：按出现顺序
  as.integer(factor(x, levels = unique(x)))
}

interval_df <- sample_table %>%
  mutate(time_order = make_time_order(.data[[TIME_COL]])) %>%
  arrange(.data[[PATIENT_COL]], time_order) %>%
  group_by(.data[[PATIENT_COL]]) %>%
  mutate(
    CST_t  = Assigned_CST_uncertain,
    CST_t1 = lead(Assigned_CST_uncertain),
    
    # interval 暴露：用下一次采样点的 ABX/INF 代表 t->t+1 区间（可按你的暴露定义替换）
    abx_interval = lead(.data[[ABX_COL]]),
    inf_interval = lead(.data[[INF_COL]])
  ) %>%
  filter(!is.na(CST_t1)) %>%
  mutate(
    # outcomes
    transition_any = as.integer(CST_t != CST_t1),
    trans_f2_u = as.integer(CST_t == "fungi_2"   & CST_t1 == "uncertain"),
    trans_u_b2 = as.integer(CST_t == "uncertain" & CST_t1 == "bacteria_2"),
    
    # denom flags
    from_f2 = as.integer(CST_t == "fungi_2"),
    from_u  = as.integer(CST_t == "uncertain")
  ) %>%
  ungroup()

# 把 CST-level 网络指标贴到 interval（按起始状态 CST_t）
interval_df <- interval_df %>%
  left_join(cst_metric_map, by = c("CST_t" = "CST")) %>%
  rename(
    density_CSTt = Density,
    modularity_CSTt = Modularity,
    centralization_CSTt = Centralization
  )

write.csv(interval_df, "./nc返修/figure3/all/interval_table_no_bootstrap.csv", row.names = FALSE)

# =========================
# 6) Basic transition probabilities (descriptive)
# =========================
p_f2u <- with(interval_df, sum(trans_f2_u, na.rm = TRUE) / sum(from_f2, na.rm = TRUE))
p_ub2 <- with(interval_df, sum(trans_u_b2, na.rm = TRUE) / sum(from_u,  na.rm = TRUE))

cat("\nTransition probability (no bootstrap):\n")
cat("P(fungi_2 -> uncertain) =", p_f2u, "\n")
cat("P(uncertain -> bacteria_2) =", p_ub2, "\n")

# =========================
# 7) Models: does CST network density predict transitions?
# =========================

# ---- Model A: fungi_2 -> uncertain (subset where CST_t == fungi_2)
df_f2 <- interval_df %>% filter(CST_t == "fungi_2")
interval_df$density_CSTt <- as.numeric(interval_df$density_CSTt)

# 如果 lme4 可用：混合效应 logistic（病人随机截距）
if (requireNamespace("lme4", quietly = TRUE)) {
  # 动态构建公式
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

# ---- Model B: uncertain -> bacteria_2 (subset where CST_t == uncertain)
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

# ---- Optional: interaction with antibiotics (addresses "after antibiotics")
if (requireNamespace("lme4", quietly = TRUE)) {
  fit_f2_int <- lme4::glmer(
    trans_f2_u ~ density_CSTt * abx_interval + inf_interval + (1 | .data[[PATIENT_COL]]),
    data = df_f2,
    family = binomial()
  )
  cat("\nGLMM with interaction: fungi_2 -> uncertain\n")
  print(summary(fit_f2_int))
}

# =========================
# 8) Simple plots (descriptive)
# =========================
# 箱线图：不同起始 CST 的 interval 暴露差异（可选）
p_abx <- ggplot(interval_df, aes(x = CST_t, y = abx_interval)) +
  geom_boxplot() + theme_bw() +
  labs(title = "Interval antibiotics exposure by starting CST", x = "Starting CST", y = "ABX (interval)")
ggsave("./nc返修/figure3/all/abx_by_startCST_no_bootstrap.pdf", p_abx, width = 6, height = 5)

p_inf <- ggplot(interval_df, aes(x = CST_t, y = inf_interval)) +
  geom_boxplot() + theme_bw() +
  labs(title = "Interval infection exposure by starting CST", x = "Starting CST", y = "Infection (interval)")
ggsave("./nc返修/figure3/all/infection_by_startCST_no_bootstrap.pdf", p_inf, width = 6, height = 5)
