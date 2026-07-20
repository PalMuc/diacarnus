library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(grid)

## 1. Read kinship (.kin0 only, plink2 output)
kin <- read.table(
  "~/Downloads/Diacarnus0526_king_OZ_nuclear.kin0",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

kin$IID1    <- trimws(as.character(kin$IID1))
kin$IID2    <- trimws(as.character(kin$IID2))
kin$KINSHIP <- as.numeric(kin$KINSHIP)

## 2. Build symmetric kinship matrix
ids <- sort(unique(c(kin$IID1, kin$IID2)))
K <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
K[cbind(kin$IID1, kin$IID2)] <- kin$KINSHIP
K[cbind(kin$IID2, kin$IID1)] <- kin$KINSHIP
diag(K) <- NA

stopifnot(sum(!is.na(K) & K != 0) > 0)

## 3. Read depth info
depth <- read.table(
  "~/Downloads/pop_info.txt",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    IID = trimws(as.character(IID)),
    pop = trimws(as.character(pop))
  ) %>%
  filter(IID %in% ids)

depth$pop <- factor(depth$pop, levels = c("Shallow", "Mesophotic"))
depth_vec <- setNames(as.character(depth$pop), depth$IID)

## 4. Relationship classes (KING thresholds)
kin$relationship <- cut(
  kin$KINSHIP,
  breaks = c(-Inf, 0.0442, 0.0884, 0.177, 0.354, Inf),
  labels = c("Unrelated", "3rd_degree", "2nd_degree", "1st_degree", "Duplicate_or_MZ")
)

## 5. Export cross-depth related pairs
cross_pairs <- kin %>%
  mutate(
    depth1 = depth_vec[IID1],
    depth2 = depth_vec[IID2]
  ) %>%
  filter(
    !is.na(depth1),
    !is.na(depth2),
    depth1 != depth2,
    KINSHIP >= 0.0442
  ) %>%
  transmute(
    IID1,
    IID2,
    depth1,
    depth2,
    KINSHIP = round(KINSHIP, 4),
    degree = as.character(relationship)
  ) %>%
  arrange(desc(KINSHIP))

print(cross_pairs)

write.table(
  cross_pairs,
  "cross_depth_related_pairs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## 6. Order by depth: Shallow first, then Mesophotic
ord <- depth %>%
  arrange(pop, IID) %>%
  pull(IID)

n <- length(ord)
gap <- sum(depth_vec[ord] == "Shallow")

K_disp <- K[ord, ord]
K_disp[K_disp < 0] <- 0
ub <- max(K_disp, na.rm = TRUE)

## 7. Convert matrix to long format
df <- as.data.frame(as.table(K_disp), stringsAsFactors = FALSE) %>%
  rename(Row = Var1, Col = Var2, KINSHIP = Freq) %>%
  mutate(
    Row = as.character(Row),
    Col = as.character(Col),
    depth_row = depth_vec[Row],
    depth_col = depth_vec[Col],
    row_idx = match(Row, ord),
    col_idx = match(Col, ord)
  )

## 8. Remove redundant parts:
##    - remove upper-right cross-depth block
##    - within Shallow/Shallow keep only lower triangle
##    - within Mesophotic/Mesophotic keep only lower triangle
df_plot <- df %>%
  filter(!(depth_row == "Shallow" & depth_col == "Mesophotic")) %>%
  filter(!(depth_row == "Shallow" &
             depth_col == "Shallow" &
             col_idx > row_idx)) %>%
  filter(!(depth_row == "Mesophotic" &
             depth_col == "Mesophotic" &
             col_idx > row_idx))

## 9. Set plotting order and marker labels
df_plot <- df_plot %>%
  mutate(
    Col = factor(Col, levels = ord),
    Row = factor(Row, levels = rev(ord)),
    x = as.numeric(Col),
    y = as.numeric(Row),
    mark = case_when(
      !is.na(KINSHIP) & depth_row != depth_col & KINSHIP >= 0.0884 ~ "**",
      !is.na(KINSHIP) & depth_row != depth_col & KINSHIP >= 0.0442 ~ "*",
      TRUE ~ ""
    )
  )

## 10. Left depth strip
left_ann <- tibble(
  sample = rev(ord),
  Depth = factor(depth_vec[rev(ord)], levels = c("Shallow", "Mesophotic")),
  x = 0,
  y = seq_len(n)
)

## 11. Bottom depth strip
bottom_ann <- tibble(
  sample = ord,
  Depth = factor(depth_vec[ord], levels = c("Shallow", "Mesophotic")),
  x = seq_len(n),
  y = -0.75
)

## 12. Plot
p <- ggplot() +
  geom_tile(
    data = df_plot,
    aes(x = x, y = y, fill = KINSHIP),
    color = "grey85",
    linewidth = 0.3
  ) +
  geom_text(
    data = df_plot %>% filter(mark != ""),
    aes(x = x, y = y, label = mark),
    color = "black",
    size = 6
  ) +
  scale_fill_gradientn(
    colours = c("#f7f7f7", "#f1b6c1", "hotpink3", "brown1"),
    values  = scales::rescale(c(0, 0.08, 0.177, 0.354)),
    limits  = c(0, max(ub, 0.354)),
    breaks  = c(0, 0.08, 0.177, 0.354),
    labels  = c(
      "Least related",
      "Mid relatedness",
      "Higher relatedness",
      "Clonality"
    ),
    name = NULL,
    na.value = "white",
    guide = guide_colorbar(
      order = 1,
      barheight = unit(90, "pt"),
      barwidth  = unit(10, "pt")
    )
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    data = left_ann,
    aes(x = x, y = y, fill = Depth),
    width = 0.85,
    show.legend = TRUE
  ) +
  geom_tile(
    data = bottom_ann,
    aes(x = x, y = y, fill = Depth),
    height = 0.45,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = c(
      Shallow = "lightblue",
      Mesophotic = "purple3"
    ),
    name = "Depth",
    guide = guide_legend(order = 2)
  ) +
  scale_x_continuous(
    breaks = seq_len(n),
    labels = ord,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq_len(n),
    labels = rev(ord),
    expand = expansion(mult = c(0.02, 0.14))
  ) +
  
  coord_fixed(clip = "off") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 18
    ),
    axis.text.y = element_text(size = 18),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    plot.margin = margin(15, 45, 50, 15)
  )+ theme(
    legend.position = c(0.63, 0.70),
    legend.justification = c("left", "center"),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

p

