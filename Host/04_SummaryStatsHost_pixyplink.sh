module load gnu/12
source ~/.bashrc
conda activate ngs-core   # pixy here

# Paths
WORKDIR=/home/lvilleg/Diacarnus/03_mapping/dedup_q30
VCF=${WORKDIR}/host_freebayes.allsites.filtered.pixy_ready.vcf.gz
POPFILE=/home/lvilleg/Diacarnus/config/pixy_pops.txt
OUTDIR=/home/lvilleg/Diacarnus/04_pixy
PREFIX=Diacarnus_shallow_mesophotic

mkdir -p "${OUTDIR}"
cd "${WORKDIR}"

# Ensure VCF is indexed
if [ ! -f "${VCF}.tbi" ]; then
    echo "Indexing pixy-ready VCF..."
    tabix -p vcf "${VCF}"
fi

# Parameters
WINDOW_SIZE=10000      # 10 kb windows; adjust if desired
N_CORES=12             # match cpus-per-task
CHUNK_SIZE=100000      # bp per chunk

echo "Running pixy (pi, dxy, fst)..."
pixy \
  --stats pi fst dxy \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size ${WINDOW_SIZE} \
  --n_cores ${N_CORES} \
  --chunk_size ${CHUNK_SIZE} \
  --output_folder "${OUTDIR}" \
  --output_prefix "${PREFIX}"

echo "pixy run finished."

module load gnu/12
source ~/.bashrc
conda activate ngs-core

WORKDIR=/home/lvilleg/Diacarnus/03_mapping/dedup_q30
OUTDIR=/home/lvilleg/Diacarnus/04_pixy
VCF=${WORKDIR}/host_freebayes.allsites.filtered.pixy_ready.vcf.gz

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# 1) VCF -> PLINK
plink \
  --vcf "${VCF}" \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out Diacarnus

# 2) Filter individuals by missingness, then PCA
#    Here: remove individuals with >50% missing genotypes
plink \
  --bfile Diacarnus \
  --allow-extra-chr \
  --mind 0.5 \
  --make-bed \
  --out Diacarnus_mind

# 3) PCA on filtered dataset
plink \
  --bfile Diacarnus_mind \
  --allow-extra-chr \
  --pca 20 \
  --out Diacarnus_pca_mind

echo "PCA finished. Results in: ${OUTDIR}"
echo "  Diacarnus_pca_mind.eigenvec  (sample scores)"
echo "  Diacarnus_pca_mind.eigenval  (eigenvalues)"



#Summary Fst checking for CI no mitochondria 

library(boot)

fst <- read.csv("~/Downloads/Diacarnus0526_shallow_mesophotic_fst.txt", sep="")
# --- filtering ---
# keep only OZ nuclear contigs; drop mitochondrion OZ024687.1
fst <- fst[grepl("^OZ", fst$chromosome) & fst$chromosome != "OZ024687.1", ]
# valid Fst values (keep negatives; drop NA)
fst <- fst[!is.na(fst$avg_wc_fst), ]

# --- SNP-count threshold (the consequential choice) ---
# windows with 1-4 SNPs contribute noise-zeros; require a minimum.
min_snps <- 5
f <- fst[fst$no_snps >= min_snps, ]

cat("windows total:", nrow(fst), " kept (>=", min_snps, "SNPs):", nrow(f), "\n")
cat("mean Fst all windows:", round(mean(fst$avg_wc_fst),4),
    " | filtered:", round(mean(f$avg_wc_fst),4), "\n")

# --- window bootstrap ---
mean_fst <- function(d, i) mean(d$avg_wc_fst[i], na.rm=TRUE)
set.seed(1)
bw <- boot(f, mean_fst, R=10000)
cat("\nWindow bootstrap:\n"); print(boot.ci(bw, type=c("perc","bca")))

# --- chromosome block bootstrap (more honest given LD within chroms) ---
chroms <- split(f, f$chromosome)
cat("n chromosomes:", length(chroms), "\n")
boot_chrom <- function() {
  samp <- sample(names(chroms), length(chroms), replace=TRUE)
  mean(sapply(samp, function(c) mean(chroms[[c]]$avg_wc_fst, na.rm=TRUE)))
}
set.seed(1)
dist <- replicate(10000, boot_chrom())
cat("Chromosome block bootstrap:\n")
print(quantile(dist, c(0.025, 0.5, 0.975)))
