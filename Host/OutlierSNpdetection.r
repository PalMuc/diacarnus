library(pcadapt)

# Read PLINK files
x <- read.pcadapt("~/Downloads/Diacarnus_admix_ld.bed", type = "bed")

# Run pcadapt with K = 1
obj <- pcadapt(x, K = 1)

# Extract p-values and FDR-adjusted q-values (BH)
pvals <- obj$pvalues
qvals <- p.adjust(pvals, method = "BH")   # BH = FDR

# Flag SNP outliers at FDR < 0.01 (1%)
outliers_idx <- which(qvals < 0.01)
length(outliers_idx)

library(dplyr)
library(readr)

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

# Subset to 1% FDR outliers
res_outliers <- res %>% filter(q < 0.01)

library(ggplot2)

# chr IDs to keep
chr_levels <- c("OZ024669.1","OZ024670.1","OZ024671.1","OZ024672.1",
                "OZ024673.1","OZ024674.1","OZ024675.1","OZ024676.1",
                "OZ024677.1","OZ024678.1","OZ024679.1","OZ024680.1",
                "OZ024681.1","OZ024682.1","OZ024683.1","OZ024684.1",
                "OZ024685.1","OZ024686.1")

res_main <- res %>%
  filter(chrom %in% chr_levels) %>%
  mutate(chrom = factor(chrom, levels = chr_levels))

# Plot threshold at p = 0.01
thr <- -log10(0.01)

ggplot(res_main, aes(x = pos, y = -log10(p))) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = thr, linetype = "dashed", colour = "red") +
  labs(x = "Position within chromosome (bp)",
       y = "-log10(p)",
       title = "pcadapt outlier scan (nuclear chromosomes, 1% FDR)") +
  theme_bw()

# Extract 1% FDR outliers (nuclear chromosomes only, if you want)
alpha <- 0.01
outliers <- res_main %>% filter(q < alpha) %>% arrange(chrom, pos)
nrow(outliers)
head(outliers)

library(pcadapt)
library(readr)

# 1) Run pcadapt with K = 2
filename <- read.pcadapt("~/Downloads/Diacarnus_admix_ld.bed", type = "bed")
x <- pcadapt(filename, K = 2)

# 2) Read FAM to confirm individual order
fam <- read_table2(
  "~/Downloads/Diacarnus_admix_ld.fam",
  col_names = c("FID","IID","PID","MID","SEX","PHENO")
)

fam$IID   # should list GW34916 ... GW34934

# 3) Build population vector matching that order
pop <- c(
  rep("Shallow", 12),    # GW34916–GW34927
  rep("Mesophotic", 6)   # GW34928–GW34933
)

stopifnot(length(pop) == nrow(fam))  # safety check

# 4) Scree plot and score plot
plot(x, option = "screeplot")
plot(x, option = "scores", pop = pop)   # PCA coloured by habitat

