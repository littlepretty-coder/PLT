
day_col <- NULL
if ("days" %in% colnames(baseline_df)) {
  day_col <- "days"
} else if ("Day" %in% colnames(baseline_df)) {
  day_col <- "Day"
} else if ("time_days" %in% colnames(baseline_df)) {
  day_col <- "time_days"
} else if ("Timepoint_num" %in% colnames(baseline_df)) {
  day_col <- "Timepoint_num"
} else {
  stop("No day column found. Please add a numeric day column in metadata, e.g. 'days'.")
}

stopifnot(all(c("PatientID", day_col) %in% colnames(baseline_df)))

## 2) Prepare data for scatterpie
##    prob columns are exactly the module names
prob_cols <- names(mapping_fixed)  # from your MVN code; length = 8
missing_probs <- setdiff(prob_cols, colnames(baseline_df))
if (length(missing_probs) > 0) {
  stop("Missing probability columns in baseline_df: ", paste(missing_probs, collapse = ", "),
       "\nMake sure you did: baseline_df <- assign_df %>% bind_cols(as.data.frame(prob0))")
}

plot_df <- baseline_df %>%
  select(SampleID, PatientID, all_of(day_col),Timepoint, all_of(prob_cols)) %>%
  filter(!is.na(PatientID), !is.na(.data[[day_col]])) %>%
  mutate(
    days = as.numeric(.data[[day_col]]),
    PatientID = as.factor(PatientID),
    y = as.numeric(fct_inorder(PatientID))  # numeric y for spacing rows
  ) %>%
  arrange(PatientID, days)
if (!requireNamespace("scatterpie", quietly = TRUE)) {
  install.packages("scatterpie")
}
library(scatterpie)

r <-0.5
mycol<-c("#CFDD97" ,"#EFCE87",'#FFF2AE','#BEBADA','#95c4cc','grey','lightpink','orange','grey90')
p_pies <- ggplot(plot_df, aes(x = days, y = y)) +
  scatterpie::geom_scatterpie(
    aes(x = days, y = y),
    cols = prob_cols,
    pie_scale = r
  ) +
  scale_y_continuous(
    breaks = sort(unique(plot_df$y)),
    labels = levels(fct_inorder(plot_df$PatientID)),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "days",
    y = "Patient",
    title = "Patient CST probability trajectories (each sample as an 8-CST pie)"
  ) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = mycol) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1)  # keep pies circular

ggsave(file.path(outdir, "Patient_CST_probability_pies_by_day.png"),
       p_pies, width = 12, height = max(4, 0.35 * nlevels(plot_df$PatientID)), dpi = 300)

print(p_pies)


legend_df <- tibble(
  CST = factor(prob_cols, levels = prob_cols),
  v = 1
)

p_legend <- ggplot(legend_df, aes(x = CST, y = v, fill = CST)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = "CST color legend", x = NULL, y = NULL)

ggsave(file.path(outdir, "CST_color_legend.png"), p_legend, width = 9, height = 2.2, dpi = 300)
print(p_legend)

stopifnot("Assigned_CST_uncertain" %in% colnames(baseline_df))

prob_cols <- names(mapping_fixed) 
r <- 0.5                   

plot_df <- baseline_df %>%
  select(SampleID, PatientID, all_of(day_col), Assigned_CST_uncertain, all_of(prob_cols)) %>%
  filter(!is.na(PatientID), !is.na(.data[[day_col]])) %>%
  mutate(
    Days = as.numeric(.data[[day_col]]),
    PatientID = as.factor(PatientID),
    y = as.numeric(forcats::fct_inorder(PatientID))
  ) %>%
  arrange(PatientID, Days)

plot_df <- plot_df %>%
  rowwise() %>%
  mutate(
    .p = list(c_across(all_of(prob_cols))),
    .p2 = list({
      p <- as.numeric(unlist(.p))
      if (all(!is.finite(p))) p <- rep(0, length(p))
      
      ord <- order(p, decreasing = TRUE)
      
      if (Assigned_CST_uncertain == "uncertain") {

        idx <- ord[1:2]
        out <- rep(0, length(p))
        s <- sum(p[idx])
        if (s <= 0) {

          out[idx] <- 0.5
        } else {
          out[idx] <- p[idx] / s
        }
        out
      } else {

        idx <- ord[1]
        out <- rep(0, length(p))
        out[idx] <- 1
        out
      }
    })
  ) %>%
  ungroup()


plot_df[prob_cols] <- do.call(rbind, plot_df$.p2)
plot_df <- plot_df %>% select(-.p, -.p2)
 

stopifnot("Timepoint" %in% colnames(plot_df))
stopifnot("Days" %in% colnames(plot_df)) 

infect_df <- plot_df %>%
  mutate(Timepoint = as.character(Timepoint)) %>%
  filter(Timepoint == "T2") %>%
  group_by(PatientID) %>%
  summarise(infect_day = min(Days, na.rm = TRUE), .groups = "drop") %>%
  left_join(plot_df %>% distinct(PatientID, y), by = "PatientID") 

library(scatterpie)

p_pies <- ggplot(plot_df, aes(x = Days, y = y)) +
  scatterpie::geom_scatterpie(
    aes(x = Days, y = y),
    cols = prob_cols,
    pie_scale = r,color = NA  
  ) +
  scale_y_continuous(
    breaks = sort(unique(plot_df$y)),
    labels = levels(forcats::fct_inorder(plot_df$PatientID)),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "Days",
    y = "Patient",
    title = "Patient CST trajectories (certain = top1 only; uncertain = top2 proportions)"
  ) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = mycol) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1)

ggsave(file.path(outdir, "Patient_CST_pies_certain_top1_uncertain_top2.png"),
       p_pies, width = 12, height = max(4, 0.35 * nlevels(plot_df$PatientID)), dpi = 300)

print(p_pies)
