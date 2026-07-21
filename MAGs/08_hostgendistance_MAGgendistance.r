library(tidyverse)

## -----------------------------
## 1. Input files and parameters
## -----------------------------
pca_file  <- "~/Downloads/Diacarnus0526_pca_pruned_OZ_mind.eigenvec"
mag_dir   <- "~/Downloads/output"
instrain_long_file <- file.path(mag_dir, "MAGs_compare_clustered_new_genomeWide_compare.tsv")

bp_cutoff <- 10000   # minimum compared bases per host pair

## -----------------------------
## 2. Host PCA -> host distance
## -----------------------------
pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

# Number of PCs to use for host genetic distance
k <- 5

pc_cols <- paste0("PC", 1:k)
pc_mat <- as.matrix(pca[, pc_cols])
rownames(pc_mat) <- pca$IID

host_dist_full <- as.matrix(dist(pc_mat, method = "euclidean"))

## -----------------------------
## 3. Read inStrain long file
## -----------------------------
instrain <- read.table(instrain_long_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

# Extract host IDs from name1/name2 (strip BAM suffix etc.)
instrain <- instrain %>%
  mutate(
    host1 = sub("\\.MAGs.*", "", name1),
    host2 = sub("\\.MAGs.*", "", name2)
  ) %>%
  select(genome, host1, host2, compared_bases_count)

## -----------------------------
## 4. Find all MAG distance files
## -----------------------------
mag_dist_files <- list.files(
  mag_dir,
  pattern = "^MAG_dist_.*\\.tsv$",
  full.names = TRUE
)

## -----------------------------
## 5. Loop over MAGs to get correlations
## -----------------------------
results <- list()

for (f in mag_dist_files) {
  mag_id <- sub("^MAG_dist_(.*)\\.tsv$", "\\1", basename(f))
  message("Processing: ", mag_id)
  
  # 5a. Read ANI distance matrix (already 1 - popANI)
  ani_mat <- tryCatch({
    as.matrix(
      read.table(
        f,
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        row.names = 1,
        stringsAsFactors = FALSE
      )
    )
  }, error = function(e) {
    warning("Failed to read ", basename(f), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(ani_mat)) next
  
  suppressWarnings(storage.mode(ani_mat) <- "numeric")
  
  # 5b. Keep only hosts present in both ANI and PCA
  common_hosts <- intersect(rownames(ani_mat), rownames(host_dist_full))
  
  if (length(common_hosts) < 4) {
    warning("Skipping ", mag_id, ": fewer than 4 overlapping hosts")
    next
  }
  
  ani_sub  <- ani_mat[common_hosts, common_hosts, drop = FALSE]
  host_sub <- host_dist_full[common_hosts, common_hosts, drop = FALSE]
  
  # 5c. For this MAG, get compared_bases_count per host pair
  instrain_mag <- instrain %>%
    filter(genome == mag_id) %>%
    mutate(
      h_low  = pmin(host1, host2),
      h_high = pmax(host1, host2)
    )
  
  # lookup table for compared_bases_count
  key_vec <- paste(instrain_mag$h_low, instrain_mag$h_high, sep = "||")
  bp_vec  <- instrain_mag$compared_bases_count
  names(bp_vec) <- key_vec
  
  get_bp <- function(h1, h2) {
    hl <- min(h1, h2)
    hh <- max(h1, h2)
    bp_vec[paste(hl, hh, sep = "||")]
  }
  
  # 5d. Mask ANI entries where compared bases < cutoff or missing
  for (i in seq_along(common_hosts)) {
    for (j in seq_along(common_hosts)) {
      if (i >= j) next  # upper triangle only
      h1 <- common_hosts[i]
      h2 <- common_hosts[j]
      bp <- get_bp(h1, h2)
      if (is.na(bp) || bp < bp_cutoff) {
        ani_sub[h1, h2] <- NA_real_
        ani_sub[h2, h1] <- NA_real_
      }
    }
  }
  
  # 5e. Use only unique host pairs (upper triangle, excluding diagonal)
  idx <- upper.tri(ani_sub, diag = FALSE)
  
  ani_vec  <- ani_sub[idx]
  host_vec <- host_sub[idx]
  
  # Remove NA pairs
  keep <- !(is.na(ani_vec) | is.na(host_vec))
  
  n_pairs_total   <- length(ani_vec)
  n_pairs_used    <- sum(keep)
  n_pairs_dropped <- n_pairs_total - n_pairs_used
  
  if (n_pairs_used < 3) {
    warning("Skipping ", mag_id, ": too few host pairs after bp cutoff")
    next
  }
  
  # 5f. Correlation test
  ct <- cor.test(host_vec[keep], ani_vec[keep], method = "pearson")
  
  results[[mag_id]] <- tibble(
    MAG               = mag_id,
    n_hosts           = length(common_hosts),
    n_pairs_total     = n_pairs_total,
    n_pairs_used      = n_pairs_used,
    n_pairs_dropped   = n_pairs_dropped,
    bp_cutoff         = bp_cutoff,
    cor_r             = unname(ct$estimate),
    p_value           = ct$p.value
  )
}

## -----------------------------
## 6. Combine, FDR, save results
## -----------------------------
results_df <- bind_rows(results) %>%
  arrange(p_value)

print(results_df, n = nrow(results_df))

results_df_fdr <- results_df %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj_fdr)

print(results_df_fdr, n = nrow(results_df_fdr))

sig_mags <- results_df_fdr %>%
  filter(p_adj_fdr < 0.05)

sig_mags

write.table(
  results_df_fdr,
  file = file.path(mag_dir, paste0("host_distance_vs_MAG_ANI_correlation_bp", bp_cutoff, ".tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## =============================
## 7. Helper: rebuild pairs for ONE MAG
## =============================

get_mag_pairs <- function(mag_id, mag_dir, host_dist_full, instrain, bp_cutoff = 10000) {
  f <- file.path(mag_dir, paste0("MAG_dist_", mag_id, ".tsv"))
  if (!file.exists(f)) {
    warning("File not found for MAG: ", mag_id)
    return(NULL)
  }
  
  ani_mat <- as.matrix(
    read.table(
      f,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      row.names = 1,
      stringsAsFactors = FALSE
    )
  )
  storage.mode(ani_mat) <- "numeric"
  
  common_hosts <- intersect(rownames(ani_mat), rownames(host_dist_full))
  if (length(common_hosts) < 4) {
    warning("Fewer than 4 overlapping hosts for MAG: ", mag_id)
    return(NULL)
  }
  
  ani_sub  <- ani_mat[common_hosts, common_hosts, drop = FALSE]
  host_sub <- host_dist_full[common_hosts, common_hosts, drop = FALSE]
  
  instrain_mag <- instrain %>%
    filter(genome == mag_id) %>%
    mutate(
      h_low  = pmin(host1, host2),
      h_high = pmax(host1, host2)
    )
  
  key_vec <- paste(instrain_mag$h_low, instrain_mag$h_high, sep = "||")
  bp_vec  <- instrain_mag$compared_bases_count
  names(bp_vec) <- key_vec
  
  get_bp <- function(h1, h2) {
    hl <- min(h1, h2)
    hh <- max(h1, h2)
    bp_vec[paste(hl, hh, sep = "||")]
  }
  
  for (i in seq_along(common_hosts)) {
    for (j in seq_along(common_hosts)) {
      if (i >= j) next
      h1 <- common_hosts[i]
      h2 <- common_hosts[j]
      bp <- get_bp(h1, h2)
      if (is.na(bp) || bp < bp_cutoff) {
        ani_sub[h1, h2] <- NA_real_
        ani_sub[h2, h1] <- NA_real_
      }
    }
  }
  
  idx <- upper.tri(ani_sub, diag = FALSE)
  ani_vec  <- ani_sub[idx]
  host_vec <- host_sub[idx]
  keep <- !(is.na(ani_vec) | is.na(host_vec))
  
  if (sum(keep) < 3) {
    warning("Too few pairs after filtering for MAG: ", mag_id)
    return(NULL)
  }
  
  tibble(
    host_dist = host_vec[keep],
    ani_dist  = ani_vec[keep]
  )
}

## =============================
## 8. Inspect specific MAGs with lm()
## =============================

#  Dehalococcoidia_6.1
mag_to_plot <- "odDiaEryt1.Dehalococcoidia_6.1"
df_pairs <- get_mag_pairs(mag_to_plot, mag_dir, host_dist_full, instrain, bp_cutoff)

if (!is.null(df_pairs)) {
  ggplot(df_pairs, aes(x = host_dist, y = ani_dist)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal()
  
  
  lm_fit <- lm(ani_dist ~ host_dist, data = df_pairs2)
  print(summary(lm_fit))
}
# Vicinamibacteria_1.1
mag_to_plot2 <- "odDiaEryt1.Vicinamibacteria_1.1"
df_pairs2 <- get_mag_pairs(mag_to_plot2, mag_dir, host_dist_full, instrain, bp_cutoff)

if (!is.null(df_pairs2)) {
  ggplot(df_pairs2, aes(x = host_dist, y = ani_dist)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal()
  
  lm_fit2 <- lm(ani_dist ~ host_dist, data = df_pairs2)
  print(summary(lm_fit2))
}

#  Thermodesulfobacteriota_1.11
mag_to_plot3 <- "odDiaEryt1.Thermodesulfobacteriota_1.1"
df_pairs <- get_mag_pairs(mag_to_plot, mag_dir, host_dist_full, instrain, bp_cutoff)

if (!is.null(df_pairs)) {
  ggplot(df_pairs, aes(x = host_dist, y = ani_dist)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal()
  
  lm_fit3 <- lm(ani_dist ~ host_dist, data = df_pairs)
  print(summary(lm_fit3))
}
