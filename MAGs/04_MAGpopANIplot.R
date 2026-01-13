# Load packages
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

## 1. Read inStrain compare file
comp <- read.table("~/Downloads/MAGs_compare_clustered_new_genomeWide_compare.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Clean sample names (strip .MAGs.dedup.q30.bam)
comp <- comp %>%
  mutate(
    name1_clean = str_replace(name1, "\\.MAGs\\.dedup\\.q30\\.bam$", ""),
    name2_clean = str_replace(name2, "\\.MAGs\\.dedup\\.q30\\.bam$", "")
  )

## 2. Read pop_info (IID = sample ID, pop = Shallow / Mesophotic)
# Adapt sep = "\t" or sep = "" depending on the actual file; this assumes whitespace / tab.
pop_info <- read.table("~/Downloads/pop_info.txt",
                       header = TRUE, stringsAsFactors = FALSE)

# Keep only the columns we need (IID, pop)
pop_info <- pop_info[, c("IID", "pop")]

## 3. Restrict metadata to samples that appear in compare file
samples_in_comp <- unique(c(comp$name1_clean, comp$name2_clean))

meta <- pop_info %>%
  filter(IID %in% samples_in_comp)

## 4. Join depths to each side of the pair and define Shallow–Mesophotic pairs
comp2 <- comp %>%
  left_join(meta, by = c("name1_clean" = "IID")) %>%
  rename(depth1 = pop) %>%
  left_join(meta, by = c("name2_clean" = "IID")) %>%
  rename(depth2 = pop) %>%
  mutate(pair_type = case_when(
    depth1 == "Shallow"   & depth2 == "Mesophotic" ~ "Shallow–Mesophotic",
    depth1 == "Mesophotic" & depth2 == "Shallow"   ~ "Shallow–Mesophotic",
    TRUE ~ "Other"
  )) %>%
  filter(pair_type == "Shallow–Mesophotic")

## 5. Summarise popANI per genome (mean, median, min)
summ <- comp2 %>%
  group_by(genome) %>%
  summarise(
    mean_popANI   = mean(popANI, na.rm = TRUE),
    median_popANI = median(popANI, na.rm = TRUE),
    min_popANI    = min(popANI, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = ends_with("popANI"),
               names_to = "Statistic",
               values_to = "popANI") %>%
  mutate(Statistic = recode(Statistic,
                            mean_popANI   = "Mean",
                            median_popANI = "Median",
                            min_popANI    = "Min"))

## 6. Plot Mesophotic vs Shallow popANI (mean, median, min)
ggplot(summ,
       aes(x = popANI, y = genome, color = Statistic, group = genome)) +
  geom_point() +
  scale_color_manual(values = c(Mean = "red", Median = "black", Min = "blue")) +
  coord_cartesian(xlim = c(0.92, 1.00)) +   # adjust limits if needed
  labs(x = "popANI", y = "Genome",
       title = "Mesophotic vs Shallow popANI (mean, median, min)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )
