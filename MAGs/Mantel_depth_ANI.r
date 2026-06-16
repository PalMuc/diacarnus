# Load packages -----------------------------------------------------------
library(vegan)   # mantel() [web:740]
library(dplyr)

# 1. File paths -----------------------------------------------------------

# Directory containing all MAG_dist_*.tsv files (46 MAGs from inStrain)
mag_dir <- "Downloads/output"

# Output files
out_table <- "/home/lvilleg/Diacarnus_0526/mantel_depthgroup_7MAGs_results.tsv"
out_audit <- "/home/lvilleg/Diacarnus_0526/mantel_depthgroup_7MAGs_sample_audit.tsv"


# 2. Define candidate MAGs (from popANI screen) --------------------------

candidates <- c(
  "odDiaEryt1.Poribacteria_3.1",
  "odDiaEryt1.Arenicellales_2.1",
  "odDiaEryt1.Poribacteria_1.1",
  "odDiaEryt1.Porisulfidus_sp_1.1",
  "odDiaEryt1.Latescibacterota_1.1",
  "odDiaEryt1.Pseudomonadales_3.1",
  "odDiaEryt1.Acidimicrobiales_1.1"
)


# 3. Depth metadata (shallow vs mesophotic) ------------------------------

depth_meta <- data.frame(
  sample_id = c("GW34916","GW34917","GW34918","GW34919","GW34920",
                "GW34921","GW34922","GW34923","GW34924","GW34925",
                "GW34926","GW34927","GW34928","GW34929","GW34930",
                "GW34931","GW34932","GW34933","GW34934"),
  depth_raw = c("10","10","10","10","10",
                "10","10","10","10","10",
                "10","12","40","40","40",
                "38-39","38-39","34","34"),
  stringsAsFactors = FALSE
)

# Convert to numeric depth safely
depth_meta$depth <- as.numeric(depth_meta$depth_raw)
depth_meta$depth[depth_meta$depth_raw == "38-39"] <- 38.5

# Grouping: shallow (<=12 m) vs mesophotic (>32 m)
depth_meta$depth_group <- ifelse(depth_meta$depth <= 12,
                                 "shallow",
                                 ifelse(depth_meta$depth > 32, "mesophotic", NA))

depth_meta$depth_group <- factor(depth_meta$depth_group,
                                 levels = c("shallow", "mesophotic"))

message("Depth groups:")
print(table(depth_meta$depth_group, useNA = "ifany"))

message("Depth metadata sample IDs:")
message(paste(depth_meta$sample_id, collapse = ", "))


# 4. Find and filter MAG distance files ----------------------------------

mag_files_all <- list.files(mag_dir,
                            pattern = "^MAG_dist_.*\\.tsv$",
                            full.names = TRUE)

if (length(mag_files_all) == 0) {
  stop("No MAG_dist_*.tsv files found in: ", mag_dir)
}

mag_names_from_files <- sub("\\.tsv$", "", basename(mag_files_all))
mag_names_from_files <- sub("^MAG_dist_", "", mag_names_from_files)

keep <- mag_names_from_files %in% candidates
mag_files <- mag_files_all[keep]

message("Found ", length(mag_files), " MAG distance files matching candidates.")
if (length(mag_files) == 0) {
  stop("None of the candidate MAGs were found in ", mag_dir)
}


# 5. Helper to safely build audit rows -----------------------------------

make_audit_df <- function(mag_name, stage_name, ids) {
  if (length(ids) == 0) return(NULL)
  data.frame(
    MAG = mag_name,
    stage = stage_name,
    sample_id = ids,
    stringsAsFactors = FALSE
  )
}


# 6. Function: Mantel vs shallow/mesophotic for one MAG ------------------

run_one_group_mantel <- function(mag_path, depth_meta) {
  mag_name <- basename(mag_path)
  mag_name <- sub("^MAG_dist_", "", mag_name)
  mag_name <- sub("\\.tsv$", "", mag_name)
  
  message("\nProcessing: ", mag_name, " ...")
  
  mag_dist <- as.matrix(
    read.table(mag_path,
               header = TRUE,
               row.names = 1,
               sep = "\t",
               check.names = FALSE)
  )
  
  mag_ids <- rownames(mag_dist)
  
  message("  MAG matrix dimensions: ", paste(dim(mag_dist), collapse = " x "))
  message("  GW34934 present in MAG matrix? ", "GW34934" %in% mag_ids)
  
  # Align MAG IDs with metadata IDs
  common_ids <- intersect(depth_meta$sample_id, mag_ids)
  
  message("  Shared IDs before NA filtering (n=", length(common_ids), "): ",
          ifelse(length(common_ids) == 0, "none", paste(common_ids, collapse = ", ")))
  
  audit_before <- make_audit_df(mag_name, "shared_before_NA", common_ids)
  
  if (length(common_ids) < 4) {
    return(list(
      result = data.frame(
        MAG = mag_name,
        n_group_ids = nrow(depth_meta),
        n_mag_ids = length(mag_ids),
        n_shared = length(common_ids),
        n_samples = length(common_ids),
        n_removed_NA = NA_integer_,
        gw34934_in_mag = "GW34934" %in% mag_ids,
        gw34934_in_shared = "GW34934" %in% common_ids,
        gw34934_kept = FALSE,
        mantel_group_r = NA_real_,
        mantel_group_p = NA_real_,
        status = "too_few_shared_samples",
        stringsAsFactors = FALSE
      ),
      audit = audit_before
    ))
  }
  
  mag_common <- mag_dist[common_ids, common_ids, drop = FALSE]
  
  # Remove samples with any NA in MAG matrix
  na_rows <- apply(is.na(mag_common), 1, any)
  na_cols <- apply(is.na(mag_common), 2, any)
  keep_ids <- rownames(mag_common)[!(na_rows | na_cols)]
  removed_ids <- setdiff(common_ids, keep_ids)
  n_removed <- length(removed_ids)
  
  message("  Shared IDs after NA filtering (n=", length(keep_ids), "): ",
          ifelse(length(keep_ids) == 0, "none", paste(keep_ids, collapse = ", ")))
  message("  Removed due to NA (n=", n_removed, "): ",
          ifelse(length(removed_ids) == 0, "none", paste(removed_ids, collapse = ", ")))
  
  audit_removed <- make_audit_df(mag_name, "removed_due_to_NA", removed_ids)
  audit_kept    <- make_audit_df(mag_name, "kept_after_NA", keep_ids)
  audit_df <- bind_rows(audit_before, audit_removed, audit_kept)
  
  if (length(keep_ids) < 4) {
    return(list(
      result = data.frame(
        MAG = mag_name,
        n_group_ids = nrow(depth_meta),
        n_mag_ids = length(mag_ids),
        n_shared = length(common_ids),
        n_samples = length(keep_ids),
        n_removed_NA = n_removed,
        gw34934_in_mag = "GW34934" %in% mag_ids,
        gw34934_in_shared = "GW34934" %in% common_ids,
        gw34934_kept = "GW34934" %in% keep_ids,
        mantel_group_r = NA_real_,
        mantel_group_p = NA_real_,
        status = "too_few_after_NA_removal",
        stringsAsFactors = FALSE
      ),
      audit = audit_df
    ))
  }
  
  mag_common <- mag_common[keep_ids, keep_ids, drop = FALSE]
  
  if (any(is.na(mag_common))) {
    return(list(
      result = data.frame(
        MAG = mag_name,
        n_group_ids = nrow(depth_meta),
        n_mag_ids = length(mag_ids),
        n_shared = length(common_ids),
        n_samples = length(keep_ids),
        n_removed_NA = n_removed,
        gw34934_in_mag = "GW34934" %in% mag_ids,
        gw34934_in_shared = "GW34934" %in% common_ids,
        gw34934_kept = "GW34934" %in% keep_ids,
        mantel_group_r = NA_real_,
        mantel_group_p = NA_real_,
        status = "NA_remaining",
        stringsAsFactors = FALSE
      ),
      audit = audit_df
    ))
  }
  
  # Build depth group vector
  depth_group_vec <- depth_meta$depth_group[match(keep_ids, depth_meta$sample_id)]
  names(depth_group_vec) <- keep_ids
  
  # Check for missing group information
  if (any(is.na(depth_group_vec))) {
    return(list(
      result = data.frame(
        MAG = mag_name,
        n_group_ids = nrow(depth_meta),
        n_mag_ids = length(mag_ids),
        n_shared = length(common_ids),
        n_samples = length(keep_ids),
        n_removed_NA = n_removed,
        gw34934_in_mag = "GW34934" %in% mag_ids,
        gw34934_in_shared = "GW34934" %in% common_ids,
        gw34934_kept = "GW34934" %in% keep_ids,
        mantel_group_r = NA_real_,
        mantel_group_p = NA_real_,
        status = "missing_group_metadata",
        stringsAsFactors = FALSE
      ),
      audit = audit_df
    ))
  }
  
  # Binary group distance: 0 same group, 1 different group
  group_dist <- outer(depth_group_vec, depth_group_vec,
                      FUN = function(a, b) as.numeric(a != b))
  rownames(group_dist) <- colnames(group_dist) <- keep_ids
  
  # Mantel: MAG vs shallow/mesophotic group distance
  set.seed(123)
  m_group <- mantel(mag_common, group_dist,
                    method = "spearman",
                    permutations = 9999)   # standard Mantel setup [web:740][web:642]
  
  message("  Mantel group r: ", round(as.numeric(m_group$statistic), 6),
          "   p: ", as.numeric(m_group$signif))
  message("  Final n: ", length(keep_ids))
  
  list(
    result = data.frame(
      MAG = mag_name,
      n_group_ids = nrow(depth_meta),
      n_mag_ids = length(mag_ids),
      n_shared = length(common_ids),
      n_samples = length(keep_ids),
      n_removed_NA = n_removed,
      gw34934_in_mag = "GW34934" %in% mag_ids,
      gw34934_in_shared = "GW34934" %in% common_ids,
      gw34934_kept = "GW34934" %in% keep_ids,
      mantel_group_r = as.numeric(m_group$statistic),
      mantel_group_p = as.numeric(m_group$signif),
      status = "ok",
      stringsAsFactors = FALSE
    ),
    audit = audit_df
  )
}


# 7. Run over all candidate MAGs -----------------------------------------

all_outputs <- lapply(mag_files, run_one_group_mantel, depth_meta = depth_meta)

results <- bind_rows(lapply(all_outputs, `[[`, "result")) %>%
  arrange(mantel_group_p)

audit <- bind_rows(lapply(all_outputs, `[[`, "audit"))


# 8. Add FDR correction ---------------------------------------------------

results <- results %>%
  mutate(
    mantel_group_p_adj = p.adjust(mantel_group_p, method = "fdr")  # BH/FDR [web:717][web:722]
  )


# 9. Write outputs --------------------------------------------------------

write.table(results, out_table, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(audit,   out_audit, sep = "\t", quote = FALSE, row.names = FALSE)


# 10. Print summary -------------------------------------------------------

message("\nSummary of shallow-vs-mesophotic Mantel tests (7 MAGs):")
print(results)

results %>%
  arrange(mantel_group_p_adj) %>%
  select(MAG,
         n_samples,
         mantel_group_r,
         mantel_group_p,
         mantel_group_p_adj,
         status) %>%
  head(7)
