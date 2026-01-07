module load gnu/12
source ~/.bashrc
conda activate ngs-core    # plink + admixture here

# ---- Paths ----
WORKDIR=/home/lvilleg/Diacarnus/03_mapping/dedup_q30
OUTDIR=/home/lvilleg/Diacarnus/05_admixture
VCF=${WORKDIR}/host_freebayes.allsites.filtered.pixy_ready.vcf.gz

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# ---- Prefixes ----
BASE=Diacarnus_admixbase        # VCF -> PLINK
MIND=Diacarnus_admix_mind       # missingness-filtered
REN=Diacarnus_admix_renamed     # unique SNP IDs
PRUNE=Diacarnus_admix_prune     # LD pruning info
LDSET=Diacarnus_admix_ld        # LD-pruned dataset (orig chrom names)
INT=Diacarnus_admix_ld_int      # LD-pruned dataset with integer chromosomes

echo "### 1) VCF -> PLINK (if needed)"
if [ ! -f ${BASE}.bed ]; then
  plink --vcf "${VCF}" \
        --double-id \
        --allow-extra-chr \
        --make-bed \
        --out ${BASE}
fi

echo "### 2) Filter individuals by missingness (mind 0.5)"
if [ ! -f ${MIND}.bed ]; then
  plink --bfile ${BASE} \
        --allow-extra-chr \
        --mind 0.5 \
        --make-bed \
        --out ${MIND}
fi

echo "### 3) Rename SNPs to unique IDs"
if [ ! -f ${REN}.bed ]; then
  awk 'BEGIN{OFS="\t"} { $2 = "SNP" NR; print }' ${MIND}.bim > ${REN}.bim
  cp ${MIND}.bed ${REN}.bed
  cp ${MIND}.fam ${REN}.fam
fi

echo "### 4) LD pruning on renamed SNPs"
plink --bfile ${REN} \
      --allow-extra-chr \
      --indep-pairwise 50 10 0.2 \
      --out ${PRUNE}


      plink --bfile ${REN} \
      --allow-extra-chr \
      --extract ${PRUNE}.prune.in \
      --make-bed \
      --out ${LDSET}

echo "### 5) Force chromosome codes to integers for ADMIXTURE"
# Here we simply number chromosomes 1..N in the order they appear
awk 'BEGIN{OFS="\t"} { $1 = NR; print }' ${LDSET}.bim > ${INT}.bim
cp ${LDSET}.bed ${INT}.bed
cp ${LDSET}.fam ${INT}.fam

INT_FINAL=${INT}

echo "### 6) Run ADMIXTURE K=1 ^`^s4"
for K in 1 2 3 4; do
  echo "Running ADMIXTURE K=${K}"
  admixture --cv=10 -j${SLURM_CPUS_PER_TASK} ${INT_FINAL}.bed ${K} | tee ${INT_FINAL}.K${K}.log
done

echo "All done. Q/P files are ${INT_FINAL}.K.Q and ${INT_FINAL}.K.P for each K."

