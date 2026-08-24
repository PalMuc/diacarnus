library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(jsonlite)

## -----------------------------
## 1. Read the raw inStrain compare table
##    This is the scaffold-level comparison table
## -----------------------------
comp <- read_csv(
  "~/Downloads/raw_data/comparisonsTable.csv.gz",
  show_col_types = FALSE
)

## -----------------------------
## 2. Clean sample names
## -----------------------------
comp <- comp %>%
  mutate(
    name1_clean = str_replace(name1, "\\.MAGs\\.dedup\\.q30\\.bam$", ""),
    name2_clean = str_replace(name2, "\\.MAGs\\.dedup\\.q30\\.bam$", "")
  )

## -----------------------------
## 3. Read scaffold-to-bin mapping
##    This maps scaffold -> MAG (bin / genome)
## -----------------------------
scaffold2bin <- fromJSON(path.expand("~/Downloads/raw_data/scaffold2bin.json"))

scaffold2bin_df <- tibble(
  scaffold = names(scaffold2bin),
  genome = unname(scaffold2bin)
)

## -----------------------------
## 4. Attach MAG / genome ID to each scaffold row
## -----------------------------
comp <- comp %>%
  left_join(scaffold2bin_df, by = "scaffold")

cat("Rows with no genome assignment:", sum(is.na(comp$genome)), "\n")

## -----------------------------
## 5. Read sample metadata
##    pop_info must contain at least: IID, pop
## -----------------------------
pop_info <- read.table(
  "~/Downloads/pop_info.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

pop_info <- pop_info[, c("IID", "pop")]

## -----------------------------
## 6. Restrict metadata to samples present in compare file
## -----------------------------
samples_in_comp <- unique(c(comp$name1_clean, comp$name2_clean))

meta <- pop_info %>%
  filter(IID %in% samples_in_comp)

## -----------------------------
## 7. Join metadata to both sides of each pair
## -----------------------------
comp2 <- comp %>%
  left_join(meta, by = c("name1_clean" = "IID")) %>%
  rename(depth1 = pop) %>%
  left_join(meta, by = c("name2_clean" = "IID")) %>%
  rename(depth2 = pop) %>%
  mutate(
    pair_type = case_when(
      depth1 == "Shallow"    & depth2 == "Mesophotic" ~ "Shallow–Mesophotic",
      depth1 == "Mesophotic" & depth2 == "Shallow"    ~ "Shallow–Mesophotic",
      depth1 == "Shallow"    & depth2 == "Shallow"    ~ "Shallow–Shallow",
      depth1 == "Mesophotic" & depth2 == "Mesophotic" ~ "Mesophotic–Mesophotic",
      TRUE ~ "Other"
    )
  )

write_tsv(comp2, "comparisonsTable_with_MAG_and_depth.tsv")

## -----------------------------
## 8. Keep only Shallow–Mesophotic pairs
## -----------------------------
comp_SM <- comp2 %>%
  filter(pair_type == "Shallow–Mesophotic") %>%
  filter(
    percent_genome_compared >= 0,
    compared_bases_count >= 10000
  )

write_tsv(comp_SM, "comparisonsTable_Shallow_Mesophotic.tsv")

## -----------------------------
## 9. Summarise per MAG (= genome)
## -----------------------------
summ <- comp_SM %>%
  filter(!is.na(genome)) %>%
  group_by(genome) %>%
  summarise(
    n_rows = n(),
    n_scaffolds = n_distinct(scaffold),
    n_sample_pairs = n_distinct(paste(name1_clean, name2_clean, sep = "__")),
    mean_popANI = mean(popANI, na.rm = TRUE),
    median_popANI = median(popANI, na.rm = TRUE),
    min_popANI = min(popANI, na.rm = TRUE),
    max_popANI = max(popANI, na.rm = TRUE),
    mean_percent_compared = mean(percent_genome_compared, na.rm = TRUE),
    median_percent_compared = median(percent_genome_compared, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    genome = str_remove(genome, "^odDiaEryt1\\.")
  )

write_tsv(summ, "Shallow_Mesophotic_MAG_summary_from_comparisonsTable.tsv")

## -----------------------------
## 10. Long format for plotting
## -----------------------------
summ_long <- summ %>%
  pivot_longer(
    cols = c(mean_popANI, median_popANI, min_popANI),
    names_to = "Statistic",
    values_to = "popANI"
  ) %>%
  mutate(
    Statistic = recode(
      Statistic,
      mean_popANI = "Mean",
      median_popANI = "Median",
      min_popANI = "Min"
    )
  )

write_tsv(summ_long, "Shallow_Mesophotic_MAG_summary_long.tsv")

## -----------------------------
## 11. Define threshold and MAGs to highlight
##     Keep original logic:
##     mean_popANI < 0.9995 and at least 6 sample pairs
## -----------------------------
threshold <- 0.9995

mags_below_any <- summ %>%
  filter(
    mean_popANI < threshold,
    n_sample_pairs >= 6
  ) %>%
  arrange(min_popANI)

write_tsv(
  mags_below_any,
  "MAGs_with_any_Shallow_Mesophotic_popANI_below_99.95_min6pairs.tsv"
)

low_ani_genomes <- unique(mags_below_any$genome)

# Use the factor levels already in the plot
gen_levels <- levels(summ_long$genome)

shade_df <- mags_below_any %>%
  mutate(
    genome = factor(genome, levels = gen_levels),
    y = as.numeric(genome),
    y_min = y - 0.5,
    y_max = y + 0.5
  )




## -----------------------------
## 12. Prepare plotting data
##     Keep ordering by Mean to match earlier plot
## -----------------------------
genome_order <- summ_long %>%
  filter(Statistic == "Mean") %>%
  arrange(popANI) %>%
  pull(genome)

summ_long_plot <- summ_long %>%
  filter(!is.na(popANI), !is.na(genome)) %>%
  mutate(
    genome = factor(genome, levels = genome_order),
    below_thresh = genome %in% low_ani_genomes
  )

## -----------------------------
## 13. Publication-style plot
## -----------------------------
p <- ggplot() +
  # soft shading for MAGs in mags_below_any
  geom_rect(
    data = shade_df,
    aes(
      xmin = -Inf, xmax = Inf,   # full row; or use xmin = 0.92, xmax = threshold
      ymin = y_min, ymax = y_max
    ),
    fill = "grey90",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = threshold,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  geom_point(
    data = summ_long_plot,
    aes(
      x = popANI,
      y = genome,
      color = Statistic,
      alpha = below_thresh
    ),
    size = 2.2
  ) +
  scale_color_manual(
    values = c(
      Mean = "#1f78b4",
      Median = "black",
      Min = "#ff7f00"
    ),
    name = NULL
  ) +
  scale_alpha_manual(
    values = c(`TRUE` = 1, `FALSE` = 0.28),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0.92, 1.00)) +
  labs(
    x = "popANI (Shallow vs Mesophotic)",
    y = "MAG",
    title = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

print(p)



threshold <- 0.9995  # 99.95%

mags_below_any <- summ %>%
  filter(
    mean_popANI < threshold,     # clearly distinct strains
    n_sample_pairs >= 6        # at least 10 Shallow–Mesophotic comparisons
  ) %>%
  arrange(min_popANI)

write_tsv(
  mags_below_any,
  "MAGs_with_any_Shallow_Mesophotic_popANI_below_99.95_min6pairs.tsv"
)
