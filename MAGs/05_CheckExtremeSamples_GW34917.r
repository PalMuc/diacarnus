library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggrepel)

## =========================
## PCA / metadata
## =========================

# PLINK 1.x eigenvec: no header, columns = FID IID PC1 PC2 PC3 ...
eigenvectors <- read.table(
  "~/Downloads/Diacarnus_pca_pixybase_mind05_geno01_maf005_pruned.eigenvec",
  header = FALSE, stringsAsFactors = FALSE
)

# name columns: adjust number of PCs (here assume 10; change if needed)
n_pc <- ncol(eigenvectors) - 2
pc_names <- paste0("PC", seq_len(n_pc))
colnames(eigenvectors) <- c("FID", "IID", pc_names)

eigenvalues <- read.table(
  "~/Downloads/Diacarnus_pca_pixybase_mind05_geno01_maf005_pruned.eigenval",
  header = FALSE, stringsAsFactors = FALSE
)

pop_info <- read.table(
  "~/Downloads/pop_info.txt",
  header = TRUE, stringsAsFactors = FALSE
)[, c("IID", "pop")]

## Host PCA scores + pop from pop_info
host_scores <- eigenvectors %>%
  select(IID, Host_PC1 = PC1, Host_PC2 = PC2) %>%
  inner_join(pop_info, by = "IID")

comp <- read.table(
  "~/Downloads/MAGs_compare_clustered_new_genomeWide_compare.tsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

## =========================
## 1) ALL MAGs together
## =========================

# (A) use ALL MAGs
comp_sel <- comp

# OR: restrict to certain genomes:
# mag_ids <- c("odDiaEryt1.Alphaproteobacteria_5.1",
#              "odDiaEryt1.Poribacteria_2.1")
# comp_sel <- comp %>% filter(genome %in% mag_ids)

# Clean sample IDs
comp_sel <- comp_sel %>%
  mutate(
    name1_clean = str_replace(name1, "\\.MAGs\\.dedup\\.q30\\.bam$", ""),
    name2_clean = str_replace(name2, "\\.MAGs\\.dedup\\.q30\\.bam$", "")
  )

# 1 - popANI as distance, aggregated across all MAGs per host pair
dist_tab <- comp_sel %>%
  transmute(
    sample1 = name1_clean,
    sample2 = name2_clean,
    dist = 1 - popANI
  ) %>%
  group_by(sample1, sample2) %>%
  summarise(dist = mean(dist, na.rm = TRUE), .groups = "drop")

# Wide distance matrix
mat <- dist_tab %>%
  pivot_wider(names_from = sample2, values_from = dist) %>%
  column_to_rownames("sample1") %>%
  as.matrix()

# Ensure symmetry and zeros on diagonal
mat[is.na(mat)] <- t(mat)[is.na(mat)]
mat[is.na(mat)] <- 0
mat <- (mat + t(mat)) / 2
diag(mat) <- 0

# PCoA
d    <- as.dist(mat)
pcoa <- cmdscale(d, k = 2, eig = TRUE)

bac_scores <- data.frame(
  IID     = rownames(pcoa$points),
  Bac_PC1 = pcoa$points[, 1],
  Bac_PC2 = pcoa$points[, 2]
)

combined <- host_scores %>%
  inner_join(bac_scores, by = "IID")

# Check GW34917
combined %>% filter(IID == "GW34917")

# (A) Bacterial ordination
ggplot(combined, aes(Bac_PC1, Bac_PC2, color = pop)) +
  geom_point(alpha = 0.7) +
  geom_point(
    data = subset(combined, IID == "GW34917"),
    color = "black", size = 3, shape = 8
  ) +
  geom_text_repel(
    data = subset(combined, IID == "GW34917"),
    aes(label = IID),
    color = "black", size = 3
  ) +
  theme_minimal()

# (B) Host–bacterium alignment
ggplot(combined, aes(Host_PC1, Bac_PC1, color = pop)) +
  geom_point(alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  geom_point(
    data = subset(combined, IID == "GW34917"),
    color = "black", size = 3, shape = 8
  ) +
  geom_text_repel(
    data = subset(combined, IID == "GW34917"),
    aes(label = IID),
    color = "black", size = 3
  ) +
  xlab("Host PC1") + ylab("MAGS PCoA1") +
  theme_minimal()

## =========================
## 2) Poribacteria only
## =========================

# Define Poribacteria genomes
# Check first: unique(comp$genome)
pori_ids <- c("odDiaEryt1.Poribacteria_2.1")

comp_pori <- comp %>%
  filter(genome %in% pori_ids) %>%
  mutate(
    name1_clean = str_replace(name1, "\\.MAGs\\.dedup\\.q30\\.bam$", ""),
    name2_clean = str_replace(name2, "\\.MAGs\\.dedup\\.q30\\.bam$", "")
  )

dist_tab_pori <- comp_pori %>%
  transmute(
    sample1 = name1_clean,
    sample2 = name2_clean,
    dist = 1 - popANI
  ) %>%
  group_by(sample1, sample2) %>%
  summarise(dist = mean(dist, na.rm = TRUE), .groups = "drop")

mat_pori <- dist_tab_pori %>%
  pivot_wider(names_from = sample2, values_from = dist) %>%
  column_to_rownames("sample1") %>%
  as.matrix()

mat_pori[is.na(mat_pori)] <- t(mat_pori)[is.na(mat_pori)]
mat_pori[is.na(mat_pori)] <- 0
mat_pori <- (mat_pori + t(mat_pori)) / 2
diag(mat_pori) <- 0

d_pori    <- as.dist(mat_pori)
pcoa_pori <- cmdscale(d_pori, k = 2, eig = TRUE)

pori_scores <- data.frame(
  IID      = rownames(pcoa_pori$points),
  Pori_PC1 = pcoa_pori$points[, 1],
  Pori_PC2 = pcoa_pori$points[, 2]
)

comb_pori <- host_scores %>%
  inner_join(pori_scores, by = "IID")

# Poribacteria ordination
ggplot(comb_pori, aes(Pori_PC1, Pori_PC2, color = pop)) +
  geom_point(alpha = 0.7) +
  geom_point(
    data = subset(comb_pori, IID == "GW34917"),
    color = "black", size = 3, shape = 8
  ) +
  geom_text_repel(
    data = subset(comb_pori, IID == "GW34917"),
    aes(label = IID),
    color = "black", size = 3
  ) +
  theme_minimal()

ggplot(comb_pori, aes(Pori_PC1, Pori_PC2, color = pop)) +
  geom_point(alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  geom_point(
    data = subset(comb_pori, IID == "GW34917"),
    color = "black", size = 3, shape = 8
  ) +
  geom_text_repel(
    data = subset(comb_pori, IID == "GW34917"),
    aes(label = IID),
    color = "black", size = 3
  ) +
  theme_minimal()
