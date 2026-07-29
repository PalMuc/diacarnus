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



#Summary Fst only on OZ chromosomes and no mitochondria CI in R

# 18 per-chromosome Fst values (column 3, the one giving 0.069)
fst_chrom <- c(0.076772,0.047276,0.069387,0.078374,0.066238,0.049785,
               0.056196,0.058256,0.088475,0.075211,0.077685,0.065401,
               0.084756,0.068958,0.057955,0.072826,0.062575,0.097436)

cat("n:", length(fst_chrom),
    " mean:", round(mean(fst_chrom),4),
    " sd:", round(sd(fst_chrom),4), "\n")

set.seed(1)
dist <- replicate(10000, mean(sample(fst_chrom, length(fst_chrom), replace=TRUE)))
print(round(quantile(dist, c(0.025, 0.5, 0.975)), 4))

# BCa version via boot, for a more principled interval at n=18
library(boot)
b <- boot(fst_chrom, function(d,i) mean(d[i]), R=10000)
print(boot.ci(b, type=c("perc","bca")))


