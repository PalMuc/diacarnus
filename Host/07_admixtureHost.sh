

WORKDIR=/home/lvilleg/Diacarnus_0526/03_mapping/dedup_q30
VCF=${WORKDIR}/host_freebayes.allsites.filtered.pixy_ready.vcf.gz

OUT_PLINK=/home/lvilleg/Diacarnus_0526/05_admixture_mind02_intchr_full
OUT_STRUCT=/home/lvilleg/Diacarnus_0526/06_structure_mind02_intchr_full

mkdir -p "${OUT_PLINK}" "${OUT_STRUCT}"

module load gnu/12
source ~/.bashrc
conda activate ngs-core

cd "${WORKDIR}"

echo "## STEP 1: VCF -> PLINK binary (all individuals)"
plink \
  --vcf "${VCF}" \
  --make-bed \
  --double-id \
  --allow-extra-chr \
  --out "${OUT_PLINK}/Diacarnus0526_all"

echo "## STEP 2: basic SNP QC with --mind 0.2 (keep most individuals)"
plink \
  --bfile "${OUT_PLINK}/Diacarnus0526_all" \
  --geno 0.1 \
  --mind 0.2 \
  --maf 0.05 \
  --allow-extra-chr \
  --make-bed \
  --out "${OUT_PLINK}/Diacarnus0526_all_filt_mind02"

# Optional: per-individual missingness report for QC
plink \
  --bfile "${OUT_PLINK}/Diacarnus0526_all" \
  --missing \
  --allow-extra-chr \
  --out "${OUT_PLINK}/Diacarnus0526_all_missing"

echo "Individuals before/after filtering:"
wc -l "${OUT_PLINK}/Diacarnus0526_all.fam"
wc -l "${OUT_PLINK}/Diacarnus0526_all_filt_mind02.fam"

echo "## STEP 3: LD pruning"
plink \
  --bfile "${OUT_PLINK}/Diacarnus0526_all_filt_mind02" \
  --indep-pairwise 50 5 0.2 \
  --allow-extra-chr \
  --out "${OUT_PLINK}/Diacarnus0526_ldprune_mind02"

plink \
  --bfile "${OUT_PLINK}/Diacarnus0526_all_filt_mind02" \
  --extract "${OUT_PLINK}/Diacarnus0526_ldprune_mind02.prune.in" \
  --allow-extra-chr \
  --make-bed \
  --out "${OUT_PLINK}/Diacarnus0526_admix_ld_mind02"

echo "Individuals in LD-pruned dataset:"
wc -l "${OUT_PLINK}/Diacarnus0526_admix_ld_mind02.fam"

###########################################
# 4) Fix chromosome codes for ADMIXTURE   #
###########################################

cd "${OUT_PLINK}"

BASE=Diacarnus0526_admix_ld_mind02
INT=${BASE}.intchr

echo "## STEP 4: Create ADMIXTURE-ready PLINK files with integer chr codes"

# Rewrite chr codes to 1..N for ADMIXTURE (only BIM first column)
awk 'BEGIN{OFS="\t"} { $1 = NR; print }' ${BASE}.bim > ${INT}.bim
cp ${BASE}.bed ${INT}.bed
cp ${BASE}.fam ${INT}.fam

echo "Sanity check on INTCHR dataset:"
wc -l ${INT}.fam
head -n 5 ${INT}.bim

###########################################
# 5) PLINK -> STRUCTURE export (optional) #
###########################################

plink \
  --bfile ${BASE} \
  --allow-extra-chr \
  --recode structure \
  --out Diacarnus0526_for_structure_diploid_mind02_intchrbase

cp Diacarnus0526_for_structure_diploid_mind02_intchrbase.recode.strct_in "${OUT_STRUCT}/"

###########################################
# 6) Run ADMIXTURE K=1–4 on INTCHR set   #
###########################################

echo "## STEP 6: Run ADMIXTURE K=1–4 with ${SLURM_CPUS_PER_TASK} threads on ${INT}"

for K in 1 2 3 4; do
  echo "Running ADMIXTURE K=${K}"
  admixture --cv=10 -j${SLURM_CPUS_PER_TASK} ${INT}.bed ${K} \
    | tee ${INT}.K${K}.log
done

echo "ADMIXTURE finished. Q/P files:"
ls ${INT}.K*.Q ${INT}.K*.P

echo "All steps complete."
echo "PLINK/ADMIXTURE outputs are in: ${OUT_PLINK}"
echo "STRUCTURE input is in: ${OUT_STRUCT}"

