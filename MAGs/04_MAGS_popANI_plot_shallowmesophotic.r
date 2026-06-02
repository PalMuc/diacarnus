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
scaffold2bin <- fromJSON(
  "~/Downloads/raw_data/scaffold2bin.json"
)

scaffold2bin_df <- tibble(
  scaffold = names(scaffold2bin),
  genome = unname(scaffold2bin)
)

## -----------------------------
## 4. Attach MAG / genome ID to each scaffold row
## -----------------------------
comp <- comp %>%
  left_join(scaffold2bin_df, by = "scaffold")

## Optional check
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

## Save annotated scaffold-level table
write_tsv(comp2, "comparisonsTable_with_MAG_and_depth.tsv")

## -----------------------------
## 8. Keep only Shallow–Mesophotic pairs
## -----------------------------

comp_SM <- comp2 %>%
  filter(pair_type == "Shallow–Mesophotic") %>%
  filter(percent_genome_compared >= 0,       # ommit non-comparisons
         compared_bases_count >= 1000)   # and at least 1000 bases

write_tsv(comp_SM, "comparisonsTable_Shallow_Mesophotic.tsv")

## -----------------------------
## 9. Summarise per MAG (= genome)
##    Since raw compare is scaffold-level, this aggregates all scaffolds
##    belonging to the same MAG
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

## 11. Plot (robust version)

# 11a. Decide order of MAGs (e.g. by mean popANI)
genome_order <- summ_long %>%
  filter(Statistic == "Mean") %>%
  arrange(popANI) %>%
  pull(genome)

summ_long_plot <- summ_long %>%
  filter(!is.na(popANI), !is.na(genome)) %>%
  mutate(genome = factor(genome, levels = genome_order))

# 11b. Plot
p <- ggplot(
  summ_long_plot,
  aes(x = popANI, y = genome, color = Statistic)
) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = c(Mean = "red", Median = "black", Min = "blue")) +
  coord_cartesian(xlim = c(0.92, 1.00)) +
  labs(
    x = "popANI",
    y = "MAG",
    title = "Shallow vs Mesophotic popANI per MAG",
    subtitle = "Derived from inStrain comparisonsTable.csv.gz using scaffold2bin mapping"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

print(p)

ggsave(
  filename = "Shallow_Mesophotic_popANI_perMAG_from_comparisonsTable.png",
  plot = p,
  width = 8,
  height = 10,
  dpi = 300
)


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
