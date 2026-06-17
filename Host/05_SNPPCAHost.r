library(ggplot2)
library(dplyr)

## -----------------------------
## 1. Read input files
## -----------------------------
eigenvectors <- read.table(
  "~/Downloads/Diacarnus_pca_mind.eigenvec",
  quote = "\"",
  comment.char = "",
  stringsAsFactors = FALSE
)

eigenvalues <- read.table(
  "~/Downloads/Diacarnus_pca_mind.eigenval",
  quote = "\"",
  comment.char = "",
  stringsAsFactors = FALSE
)

pop_info <- read.csv(
  "~/Downloads/pop_info.txt",
  sep = "",
  stringsAsFactors = FALSE
)

## -----------------------------
## 2. Clean column names
## -----------------------------
colnames(eigenvectors)[1:2] <- c("FID", "IID")

n_pcs <- ncol(eigenvectors) - 2
pc_names <- paste0("PC", seq_len(n_pcs))
colnames(eigenvectors)[3:ncol(eigenvectors)] <- pc_names

## -----------------------------
## 3. Join metadata and define depth
## -----------------------------
eigenvectors <- eigenvectors %>%
  left_join(pop_info, by = "IID")

# pop is your shallow/mesophotic indicator
eigenvectors$pop   <- factor(eigenvectors$pop)
eigenvectors$depth <- eigenvectors$pop  # alias for clarity

## -----------------------------
## 4. Axis labels with explained variance
## -----------------------------
pc1_lab <- paste0(
  "PC1 (",
  round(100 * eigenvalues[1, 1] / sum(eigenvalues[, 1]), 1),
  "%)"
)

pc2_lab <- paste0(
  "PC2 (",
  round(100 * eigenvalues[2, 1] / sum(eigenvalues[, 1]), 1),
  "%)"
)

## -----------------------------
## 5. PCA plot by depth (shallow/mesophotic)
## -----------------------------
p_pca <- ggplot(eigenvectors, aes(x = PC1, y = PC2, color = depth)) +
  geom_point(alpha = 0.7, size = 3) +
  stat_ellipse(level = 0.95, linewidth = 0.7) +
  xlab(pc1_lab) +
  ylab(pc2_lab) +
  theme_minimal()

print(p_pca)

## -----------------------------
## 6. Nonparametric tests (pop = depth groups)
## -----------------------------
# Kruskal–Wallis (works for 2+ groups; here you have 2)
kw_pc1 <- kruskal.test(PC1 ~ depth, data = eigenvectors)
kw_pc2 <- kruskal.test(PC2 ~ depth, data = eigenvectors)

print(kw_pc1)
print(kw_pc2)

# For two groups, also do Wilcoxon rank-sum
wilcox_pc1 <- wilcox.test(PC1 ~ depth, data = eigenvectors, exact = FALSE)
wilcox_pc2 <- wilcox.test(PC2 ~ depth, data = eigenvectors, exact = FALSE)

print(wilcox_pc1)
print(wilcox_pc2)

## -----------------------------
## 7. Boxplots for interpretation
## -----------------------------
p_pc1 <- ggplot(eigenvectors, aes(x = depth, y = PC1, fill = depth)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.7) +
  xlab("Depth (shallow vs mesophotic)") +
  ylab(pc1_lab) +
  theme_minimal() +
  theme(legend.position = "none")

p_pc2 <- ggplot(eigenvectors, aes(x = depth, y = PC2, fill = depth)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.7) +
  xlab("Depth (shallow vs mesophotic)") +
  ylab(pc2_lab) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_pc1)
print(p_pc2)

#Aditional testing using PERMANOVA
# host_dist_full: matrix of Euclidean distances in PCA space
host_dist <- as.dist(host_dist_full)

depth <- eigenvectors$pop  # shallow/mesophotic matched by IID order

permanova_host <- adonis2(host_dist ~ depth, permutations = 9999)
permanova_host

