library(ggplot2)
library(dplyr)

## -----------------------------
## 1. Read input files
## -----------------------------
###PCA plot is run both on the direct pca output and the LD prunned output ###

eigenvectors <- read.table(
  "~/Downloads/Diacarnus0526_pca_pruned_OZ_mind.eigenvec",
  quote = "\"",
  comment.char = "",
  stringsAsFactors = FALSE
)

eigenvalues <- read.table(
  "~/Downloads/Diacarnus0526_pca_pruned_OZ_mind.eigenval",
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

library(ggplot2)
library(dplyr)

# Define depth colors (adjust to match your other figures)
depth_cols <- c(
  Shallow   = "lightblue",  # blue
  Mesophotic = "purple3"  # green
)
p_pca <- ggplot(eigenvectors, aes(x = PC1, y = PC2, color = depth)) +
  geom_point(alpha = 0.8, size = 2.8) +
  stat_ellipse(level = 0.95, linewidth = 0.7, linetype = "solid") +
  scale_color_manual(values = depth_cols, name = "Depth") +
  xlab(pc1_lab) +
  ylab(pc2_lab) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )
p_pca

## -----------------------------
## 6. Overall test for a depth signal across all PCs (PERMANOVA)
## -----------------------------
library(vegan)

# Use all retained PCs as the multivariate space
pc_matrix <- eigenvectors[, pc_names]
euclid    <- dist(pc_matrix, method = "euclidean")

set.seed(1)  # reproducible permutations
permanova <- adonis2(euclid ~ depth, data = eigenvectors, permutations = 9999)
print(permanova)

#Test that PCA is not driven by coverage differences

cov <- data.frame(
  sample = c("GW34916","GW34917","GW34918","GW34919","GW34920","GW34921",
             "GW34922","GW34923","GW34924","GW34925","GW34926","GW34927",
             "GW34928","GW34929","GW34930","GW34931","GW34932","GW34933","GW34934"),
  mean_cov = c(6.84855,7.29631,3.90866,4.93968,6.98663,6.28008,
               6.79648,3.69138,2.53356,3.80445,3.21863,2.71965,
               2.93978,3.8784,3.49440,4.84000,4.37840,3.41079,4.75628)
)
m <- merge(df, cov, by.x="id", by.y="sample")
nrow(m)   # should be 19 — if fewer, an ID didn't match

wilcox.test(mean_cov ~ depth, data=m)
summary(lm(PC1 ~ mean_cov, data=m))
cor.test(m$PC1, m$mean_cov)
