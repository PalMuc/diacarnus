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

##EXTRA ANALYSIS: corroborate outliers are related to depht by hihg Fst regions##

#in bash
cd /home/lvilleg/Diacarnus_0526/03_mapping/dedup_q30

VCF=host_freebayes.variant.filtered.vcf.gz

# pop files (regenerate so the job is self-contained)
printf "GW%d\n" $(seq 34916 34927) > shallow.txt      # 12 shallow
printf "GW%d\n" $(seq 34928 34934) > mesophotic.txt   # 7 mesophotic

# guard: every pop-file sample must be present in the VCF
comm -23 <(sort shallow.txt mesophotic.txt) \
         <(bcftools query -l "$VCF" | sort) > missing_samples.txt
if [ -s missing_samples.txt ]; then
  echo "ERROR: these pop-file samples are not in the VCF:"; cat missing_samples.txt; exit 1
fi

vcftools --gzvcf "$VCF" \
  --weir-fst-pop shallow.txt \
  --weir-fst-pop mesophotic.txt \
  --out diacarnus_persnp_fst

#in R

library(dplyr)
library(readr)
library(ggplot2)


## 1. Per-SNP Fst from vcftools -------------------------------------
fst <- read_tsv("~/Downloads/diacarnus_persnp_fst.weir.fst",
                col_types = cols(CHROM = col_character(),
                                 POS   = col_double(),
                                 WEIR_AND_COCKERHAM_FST = col_double())) %>%
  rename(chrom = CHROM, pos = POS, fst = WEIR_AND_COCKERHAM_FST) %>%
  filter(!is.na(fst))

## 2. Merge Fst onto the SNPs pcadapt tested ------------------------
dat <- res_main %>%
  mutate(chrom = as.character(chrom), pos = as.numeric(pos)) %>%
  left_join(fst, by = c("chrom", "pos")) %>%            # left_join keeps all tested SNPs
  mutate(group = ifelse(q < 0.05, "outlier", "background"))

## 3. Diagnostic: how many outliers matched an Fst value? -----------
dat %>%
  mutate(has_fst = !is.na(fst)) %>%
  count(group, has_fst)
# tells you how many outlier / background SNPs got an Fst; the
# outlier + has_fst=FALSE count is your ~836 dropped SNPs

## 4. Restrict to SNPs with a valid Fst for testing -----------------
d2 <- dat %>% filter(!is.na(fst), group %in% c("outlier", "background"))
table(d2$group)

## 5. Descriptives --------------------------------------------------
d2 %>% group_by(group) %>%
  summarise(n = n(),
            median_fst = median(fst),
            mean_fst   = mean(fst),
            .groups = "drop")

## 6. Permutation test on the median (primary test) ----------------
obs   <- median(d2$fst[d2$group == "outlier"])
bg    <- median(d2$fst[d2$group == "background"])
n_out <- sum(d2$group == "outlier")
pool  <- d2$fst

set.seed(1)
null   <- replicate(10000, median(sample(pool, n_out)))
p_perm <- (sum(null >= obs) + 1) / (length(null) + 1)

cat("outlier median Fst:", round(obs, 4),
    " | background median:", round(bg, 4),
    " | permutation p:", signif(p_perm, 3), "\n")

## 7. Wilcoxon, correct direction (secondary) ----------------------
# 'less' tests whether the first alphabetical level (background) < outlier
wilcox.test(fst ~ group, data = d2, alternative = "less")

## 8. Plot ----------------------------------------------------------
ggplot(d2, aes(fst, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = obs, linetype = "dashed") +
  labs(x = "Per-SNP Weir & Cockerham Fst (shallow vs mesophotic)",
       y = "Density", fill = NULL) +
  theme_bw()
