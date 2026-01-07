library(tidyverse)

# Read the tab-delimited file exported from inStrain
# Adjust file name as needed
df <- read.delim("~/Downloads/MAG_population_popANI_conANI_summary.tsv")

# Make sure pair_class is ordered as you like
df <- df %>%
  mutate(
    pair_class = factor(pair_class,
                        levels = c("Mesophotic", "Shallow", "Mesophotic_vs_Shallow"))
  )


df_long <- df %>%
  select(genome, pair_class,
         mean_popANI, median_popANI, min_popANI) %>%
  pivot_longer(cols = c(mean_popANI, median_popANI, min_popANI),
               names_to = "stat", values_to = "popANI")

all_comparissons <-ggplot(df_long, aes(x = pair_class, y = popANI, color = stat, group = stat)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_line(position = position_dodge(width = 0.4)) +
  facet_wrap(~ genome, ncol = 3) +
  scale_color_manual(values = c(mean_popANI = "red",
                                median_popANI = "black",
                                min_popANI = "blue"),
                     labels = c(mean_popANI = "Mean",
                                median_popANI = "Median",
                                min_popANI = "Min")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "popANI", x = "Pair class",
       color = "Statistic",
       title = "popANI summary statistics by genome and pair class")

###Only comparsson between pops###

df_mvs <- df %>%
  filter(pair_class == "Mesophotic_vs_Shallow")

# Long format for popANI
pop_long <- df_mvs %>%
  select(genome, mean_popANI, median_popANI, min_popANI) %>%
  pivot_longer(cols = c(mean_popANI, median_popANI, min_popANI),
               names_to = "stat", values_to = "popANI")

ggplot(pop_long,
       aes(x = genome, y = popANI, color = stat, group = stat)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
  coord_flip() +
  scale_color_manual(values = c(mean_popANI = "red",
                                median_popANI = "black",
                                min_popANI = "blue"),
                     labels = c(mean_popANI = "Mean",
                                median_popANI = "Median",
                                min_popANI = "Min")) +
  theme_bw() +
  labs(x = "Genome", y = "popANI",
       color = "Statistic",
       title = "Mesophotic vs Shallow popANI (mean, median, min)")
