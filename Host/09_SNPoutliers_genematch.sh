#Run BUSCO to get ortholohue coordinates on the base genome

module load gnu/12
source ~/.bashrc
conda activate busco

IN_FA=/home/lvilleg/Diacarnus/00_Genome/Diacarnus_genomic.fa
OUT_DIR=/home/lvilleg/Diacarnus/09_busco_Diacarnus
DL_DIR=/home/lvilleg/Diacarnus/09_busco_downloads

mkdir -p "$OUT_DIR" "$DL_DIR"
cd "$OUT_DIR"

# Optional: only first time
# busco --download metazoa_odb12 --download_path "$DL_DIR"

busco \
  --in "$IN_FA" \
  --out Diacarnus_busco_metazoa_odb12 \
  --mode genome \
  --lineage_dataset metazoa_odb12 \
  --download_path "$DL_DIR" \
  --cpu "$SLURM_CPUS_PER_TASK"

echo "BUSCO Diacarnus (metazoa_odb12) finished."

#Then on R match the previouls identified outliers from script 08 to the genes. In this case the 08 script has been modified to scan within a flase discovery rate of 5% instead of 1%

library(pcadapt)
library(dplyr)
library(readr)
library(ggplot2)

### 1. pcadapt and FDR (5%)

# Read PLINK files
x <- read.pcadapt("~/Downloads/Diacarnus_admix_ld.bed", type = "bed")

# Run pcadapt with K = 1
obj <- pcadapt(x, K = 1)

# Extract p-values and FDR-adjusted q-values (BH)
pvals <- obj$pvalues
qvals <- p.adjust(pvals, method = "BH")   # BH = FDR

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

# (Optional) If you still want to look at 1% FDR separately:
# res_outliers_1pct <- res %>% filter(q < 0.01)

### 3. Keep nuclear chromosomes and plot

chr_levels <- c("OZ024669.1","OZ024670.1","OZ024671.1","OZ024672.1",
                "OZ024673.1","OZ024674.1","OZ024675.1","OZ024676.1",
                "OZ024677.1","OZ024678.1","OZ024679.1","OZ024680.1",
                "OZ024681.1","OZ024682.1","OZ024683.1","OZ024684.1",
                "OZ024685.1","OZ024686.1")

res_main <- res %>%
  filter(chrom %in% chr_levels) %>%
  mutate(chrom = factor(chrom, levels = chr_levels))

# Simple p-value threshold line at p = 0.01 (for visualization only)
thr <- -log10(0.01)

ggplot(res_main, aes(x = pos, y = -log10(p))) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = thr, linetype = "dashed", colour = "red") +
  labs(x = "Position within chromosome (bp)",
       y = "-log10(p)",
       title = "pcadapt outlier scan (nuclear chromosomes, FDR = 5%)") +
  theme_bw()

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
  filter(Status %in% c("Complete", "Fragmented"),
         !is.na(Sequence),
         !is.na(Gene.Start),
         !is.na(Gene.End)) %>%
  transmute(
    busco_id   = Busco.id,
    chrom      = as.character(Sequence),
    gene_start = as.numeric(Gene.Start),
    gene_end   = as.numeric(Gene.End),
    strand     = Strand,
    score      = Score,
    length     = Length,
    description = Description
  )

### 6. Match 5% FDR outliers to BUSCO genes

outliers <- outliers %>%
  mutate(
    chrom = as.character(chrom),
    pos   = as.numeric(pos)
  )

outliers_busco <- outliers %>%
  inner_join(busco, by = "chrom", relationship = "many-to-many") %>%
  filter(pos >= gene_start, pos <= gene_end)

# how many SNPs per BUSCO gene?
outliers_busco %>%
  count(busco_id, sort = TRUE)

# export table of 5% FDR outliers with BUSCO annotation
write.table(
  outliers_busco,
  file = "outliers_with_BUSCO_genes_FDR0.05.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

