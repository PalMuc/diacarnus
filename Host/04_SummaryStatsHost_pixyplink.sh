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


library(dplyr)
library(boot)
fst <- read.table("~/Downloads/Diacarnus0526_shallow_mesophotic_fst.txt",
                  header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  filter(grepl("^OZ", chromosome), chromosome != "OZ024687.1", !is.na(avg_wc_fst)) %>%
  filter(no_snps >= 5)

per_chrom <- fst %>% group_by(chromosome) %>%
  summarise(chrom_fst = mean(avg_wc_fst), .groups="drop")

cat("mean:", round(mean(per_chrom$chrom_fst),4),
    " SD:", round(sd(per_chrom$chrom_fst),4), "\n")

set.seed(1)
b <- boot(per_chrom$chrom_fst, function(d,i) mean(d[i]), R=10000)
print(boot.ci(b, type=c("perc","bca")))
