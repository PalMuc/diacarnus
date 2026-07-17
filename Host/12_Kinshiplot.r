library(tidyverse)
library(pheatmap)

## 1. Read kinship (.kin0 only, plink2 output)
kin <- read.table("~/Downloads/Diacarnus0526_king_OZ_nuclear.kin0",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
kin$IID1    <- trimws(as.character(kin$IID1))
kin$IID2    <- trimws(as.character(kin$IID2))
kin$KINSHIP <- as.numeric(kin$KINSHIP)

## 2. Symmetric kinship matrix
ids <- sort(unique(c(kin$IID1, kin$IID2)))
K <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
K[cbind(kin$IID1, kin$IID2)] <- kin$KINSHIP
K[cbind(kin$IID2, kin$IID1)] <- kin$KINSHIP
diag(K) <- NA
stopifnot(sum(!is.na(K) & K != 0) > 0)

## 3. Depth info
depth <- read.table("~/Downloads/pop_info.txt", header = TRUE,
                    stringsAsFactors = FALSE) %>%
  mutate(IID = trimws(IID), pop = trimws(pop)) %>%
  filter(IID %in% ids)
depth_vec <- setNames(depth$pop, depth$IID)

## 4. Relationship class per pair (KING thresholds)
kin$relationship <- cut(kin$KINSHIP,
                        breaks = c(-Inf, 0.0442, 0.0884, 0.177, 0.354, Inf),
                        labels = c("Unrelated","3rd_degree","2nd_degree","1st_degree","Duplicate_or_MZ"))

## 5. Table of cross-depth related pairs
cross_pairs <- kin %>%
  mutate(depth1 = depth_vec[IID1], depth2 = depth_vec[IID2]) %>%
  filter(depth1 != depth2, KINSHIP >= 0.0442) %>%
  transmute(IID1, IID2, depth1, depth2,
            KINSHIP = round(KINSHIP, 4), degree = as.character(relationship)) %>%
  arrange(desc(KINSHIP))
print(cross_pairs)
write.table(cross_pairs, "cross_depth_related_pairs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

## 6. Order by depth: Shallow block then Mesophotic block
depth$pop <- factor(depth$pop, levels = c("Shallow", "Mesophotic"))
ord <- depth %>% arrange(pop, IID) %>% pull(IID)
gap <- sum(depth_vec[ord] == levels(depth$pop)[1])
K_disp <- K[ord, ord]; K_disp[K_disp < 0] <- 0
ub <- max(K_disp, na.rm = TRUE); stopifnot(ub > 0)

## 7. Marker matrix: cross-depth related cells
mark <- matrix("", length(ord), length(ord), dimnames = list(ord, ord))
for (a in ord) for (b in ord) {
  if (a != b && depth_vec[a] != depth_vec[b] &&
      !is.na(K[a, b]) && K[a, b] >= 0.0442) {
    mark[a, b] <- ifelse(K[a, b] >= 0.0884, "**", "*")
  }
}

## 8. Annotation: depth only
annotation_df <- data.frame(
  Depth = depth_vec[ord],
  row.names = ord)
kin_colors <- colorRampPalette(c("gray97", "gray","pink","pink3", "red4"))(100)
ann_colors <- list(
  Depth = c(Shallow = "lightblue", Mesophotic = "purple3"))
p <- pheatmap(
  K_disp,
  cluster_rows = FALSE, cluster_cols = FALSE,
  gaps_row = gap, gaps_col = gap,
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  annotation_names_row = FALSE,   # <– remove "Depth" label on row side
  annotation_names_col = FALSE,   # <– remove "Depth" label on column side
  color = kin_colors,
  breaks = seq(0, ub, length.out = length(kin_colors) + 1),
  legend_breaks = c(0, ub/2, ub),
  legend_labels = c("Least related", "Mid\nrelatedness", "Highest related"),
  na_col = "grey92",
  border_color = "grey85",
  display_numbers = mark,
  number_color = "black",
  main = "",
  fontsize = 16,
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45
)

p

