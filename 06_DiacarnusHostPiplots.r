library(dplyr)
library(ggplot2)

pi_df <- read.delim("~/Downloads/Diacarnus_shallow_mesophotic_pi.txt")

pi_df_clean <- pi_df %>%
  filter(!is.na(avg_pi))


p_genome <- ggplot(pi_df_clean, aes(pop, avg_pi, fill = pop)) +
  geom_violin(trim = FALSE, linewidth = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  labs(
    title = "Genome-wide nucleotide diversity (π) per population\n(all callable windows, including π = 0)",
    x = "Population",
    y = "Nucleotide diversity (π)"
  ) +
  theme_minimal()

p_genome

pi_poly <- pi_df_clean %>%
  filter(avg_pi > 0)

p_poly <- ggplot(pi_poly, aes(pop, avg_pi, fill = pop)) +
  geom_violin(trim = FALSE, linewidth = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_y_log10() +
  labs(
    title = "Nucleotide diversity (π) per population\n(only polymorphic windows, π > 0, log10 scale)",
    x = "Population",
    y = "Nucleotide diversity (π, log10)"
  ) +
  theme_minimal()

p_poly
