library(pcadapt)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

### 1. pcadapt and FDR (5%)

# Read PLINK files
x <- read.pcadapt("~/Downloads/Diacarnus_admix_ld.bed", type = "bed")

# Run pcadapt with K = 1
obj <- pcadapt(x, K = 1)

# Extract p-values and FDR-adjusted q-values (BH)
pvals <- obj$pvalues
qvals <- p.adjust(pvals, method = "BH")

# Flag SNP outliers at FDR < 0.05 (5%)
outliers_idx <- which(qvals < 0.05)
length(outliers_idx)

### 2. Merge with BIM

bim <- read_tsv(
  "~/Downloads/Diacarnus_admix_ld.bim",
  col_names = c("chrom", "snp", "cm", "pos", "a1", "a2"),
  col_types = cols(
    chrom = col_character(),
    snp   = col_character(),
    cm    = col_double(),
    pos   = col_double(),
    a1    = col_character(),
    a2    = col_character()
  )
)

res <- tibble(
  snp_index = seq_along(pvals),
  p = pvals,
  q = qvals
) %>%
  left_join(bim %>% mutate(snp_index = row_number()), by = "snp_index")

### 3. Keep nuclear chromosomes

chr_levels <- c(
  "OZ024669.1","OZ024670.1","OZ024671.1","OZ024672.1",
  "OZ024673.1","OZ024674.1","OZ024675.1","OZ024676.1",
  "OZ024677.1","OZ024678.1","OZ024679.1","OZ024680.1",
  "OZ024681.1","OZ024682.1","OZ024683.1","OZ024684.1",
  "OZ024685.1","OZ024686.1"
)

res_main <- res %>%
  filter(chrom %in% chr_levels) %>%
  mutate(
    chrom = factor(chrom, levels = chr_levels),
    logp = -log10(p)
  ) %>%
  arrange(chrom, pos)

### 4. Extract 5% FDR outliers (nuclear chromosomes only)

alpha <- 0.05

outliers <- res_main %>%
  filter(q < alpha) %>%
  arrange(chrom, pos)

nrow(outliers)
head(outliers)

### 5. Read BUSCO table and prepare gene intervals

# first two lines and '#' line already removed manually
Diacarnus_busco <- read.delim("~/Downloads/Diacarnus_busco_full_table_edited.csv")

busco <- Diacarnus_busco %>%
  filter(
    Status %in% c("Complete", "Fragmented"),
    !is.na(Sequence),
    !is.na(Gene.Start),
    !is.na(Gene.End)
  ) %>%
  transmute(
    busco_id    = Busco.id,
    chrom       = as.character(Sequence),
    gene_start  = as.numeric(Gene.Start),
    gene_end    = as.numeric(Gene.End),
    strand      = Strand,
    score       = Score,
    length      = Length,
    description = Description
  )

### 6. Match 5% FDR outliers to BUSCO genes

outliers_busco <- outliers %>%
  mutate(
    chrom = as.character(chrom),
    pos   = as.numeric(pos)
  ) %>%
  inner_join(busco, by = "chrom", relationship = "many-to-many") %>%
  filter(pos >= gene_start, pos <= gene_end)

# how many SNPs per BUSCO gene?
outliers_busco %>%
  count(busco_id, description, sort = TRUE)

# export table of 5% FDR outliers with BUSCO annotation
write.table(
  outliers_busco,
  file = "outliers_with_BUSCO_genes_FDR0.05.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

### 7. Prepare Manhattan-style cumulative coordinates

chr_sizes <- res_main %>%
  group_by(chrom) %>%
  summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
  arrange(chrom) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

res_main <- res_main %>%
  left_join(chr_sizes, by = "chrom") %>%
  mutate(pos_cum = pos + offset)

axis_df <- chr_sizes %>%
  mutate(
    center = offset + chr_len / 2,
    chr_label = gsub("^OZ0*", "", as.character(chrom)),
    chr_label = gsub("\\.1$", "", chr_label)
  )

### 8. Add BUSCO annotation back onto all SNPs

annotated_snps <- outliers_busco %>%
  mutate(
    description_short = case_when(
      grepl("RNA polymerase II-associated protein 1", description, ignore.case = TRUE) ~ "RNAPII-assoc 1",
      grepl("glutaminyl-tRNA synthetase", description, ignore.case = TRUE) ~ "GlnRS",
      grepl("transcription initiation factor tfiid subunit 8", description, ignore.case = TRUE) ~ "TFIID subunit 8",
      grepl("anaphase-promoting complex subunit 2", description, ignore.case = TRUE) ~ "APC2",
      grepl("tbc1 domain family member 23", description, ignore.case = TRUE) ~ "TBC1D23",
      grepl("Tuftelin-interacting protein", description, ignore.case = TRUE) ~ "Tuftelin-int.",
      grepl("^GCN1$", description, ignore.case = TRUE) ~ "GCN1",
      grepl("Nuclear speckle splicing regulatory protein 1", description, ignore.case = TRUE) ~ "NSRP1",
      grepl("Hypoxia up-regulated protein 1", description, ignore.case = TRUE) ~ "HYOU1",
      grepl("intraflagellar transport protein 46", description, ignore.case = TRUE) ~ "IFT46",
      grepl("Ubiquitin-like-conjugating enzyme ATG3", description, ignore.case = TRUE) ~ "ATG3",
      grepl("Exonuclease", description, ignore.case = TRUE) ~ "Exonuclease",
      grepl("phosphatidylinositol", description, ignore.case = TRUE) ~ "Phosphatidylinositol",
      TRUE ~ description
    )
  ) %>%
  distinct(chrom, pos, .keep_all = TRUE)

res_main <- res_main %>%
  left_join(
    annotated_snps %>%
      select(chrom, pos, busco_id, description, description_short),
    by = c("chrom", "pos")
  ) %>%
  mutate(is_annotated_outlier = !is.na(busco_id))

### 9. Keep one label per gene (most significant SNP per BUSCO gene)

label_df <- res_main %>%
  filter(is_annotated_outlier) %>%
  group_by(busco_id, description_short) %>%
  slice_max(order_by = logp, n = 1, with_ties = FALSE) %>%
  ungroup()

### 10. Chromosome colours

chr_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(chr_levels))
names(chr_colors) <- chr_levels

### 11. Manhattan plot with BUSCO gene labels

thr <- -log10(0.01)

p_manhattan <- ggplot(res_main, aes(x = pos_cum, y = logp)) +
  # all SNPs
  geom_point(
    aes(color = chrom),
    alpha = 0.45,
    size = 0.4
  ) +
  # all FDR<0.05 outliers
  geom_point(
    data = res_main %>% filter(q < alpha),
    color = "grey",
    alpha = 0.75,
    size = 0.45
  ) +
  # BUSCO-annotated outlier SNPs
  geom_point(
    data = res_main %>% filter(is_annotated_outlier),
    color = "#D55E00",
    alpha = 0.95,
    size = 3
  ) +
  # one label per gene
  geom_text_repel(
    data = label_df,
    aes(label = description_short),
    size = 4.4,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.color = "grey50",
    segment.size = 0.25,
    min.segment.length = 0,
    seed = 123
  ) +
  geom_hline(
    yintercept = thr,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey35"
  ) +
  scale_color_manual(values = chr_colors, guide = "none") +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$chr_label,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p)))
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.margin = margin(5, 8, 5, 5)
  )

print(p_manhattan)

### 12. Save figure

ggsave(
  "pcadapt_manhattan_BUSCO_annotated_FDR0.05.pdf",
  plot = p_manhattan,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "pcadapt_manhattan_BUSCO_annotated_FDR0.05.png",
  plot = p_manhattan,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

