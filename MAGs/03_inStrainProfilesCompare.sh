module load gnu/12
source ~/.bashrc
conda activate instrain

REF_MAGS=/home/lvilleg/Diacarnus/00_Genome/Diacarnus_symbionts_MAGs.uniq.fa
BAM_DIR=/home/lvilleg/Diacarnus/03_mapping_MAGs/dedup_q30

BASE_DIR=/home/lvilleg/Diacarnus/04_inStrain_new
PROFILE_DIR=${BASE_DIR}/MAG_profiles_new
COMPARE_DIR=${BASE_DIR}/MAG_compare_new/MAGs_compare_clustered_new
STB_FILE=${BASE_DIR}/MAGs_new.scaff2mag.stb

mkdir -p "${PROFILE_DIR}"
mkdir -p "${COMPARE_DIR}"

echo "=== Step 1: inStrain profile (new) ==="

for bam in ${BAM_DIR}/GW34*.MAGs.dedup.q30.bam; do
    sample=$(basename "${bam}" .MAGs.dedup.q30.bam)
    outdir="${PROFILE_DIR}/${sample}_MAGs_inStrain_new"
    echo "Profiling ${sample} -> ${outdir}"
    inStrain profile "${bam}" "${REF_MAGS}" \
        -o "${outdir}" \
        -p ${SLURM_CPUS_PER_TASK}
done

echo "=== Step 2: build stb from scaffold_info (new) ==="

# Use one representative profile to get scaffold_info
SI_FILE=${PROFILE_DIR}/GW34916_MAGs_inStrain_new/output/GW34916_MAGs_inStrain_new_scaffold_info.tsv

conda activate r-instrain

Rscript - <<EOF
library(dplyr)
library(readr)

si_file <- "${SI_FILE}"
si <- read_tsv(si_file, show_col_types = FALSE)

stb <- si %>%
  transmute(
    scaffold = scaffold,
    MAG      = sub("\\\\.(primary_contig|scaffold)_.*$", "", scaffold)
  ) %>%
  distinct()

out_stb <- "${STB_FILE}"
write_tsv(stb, out_stb, col_names = FALSE)
cat("Wrote stb to:", out_stb, "\n")
EOF

echo "=== Step 3: clustered inStrain compare (new) ==="

conda activate instrain
cd "${BASE_DIR}"

inStrain compare \
  -i ${PROFILE_DIR}/GW*_MAGs_inStrain_new \
  -o "${COMPARE_DIR}" \
  --stb "${STB_FILE}" \
  -ani 0.99 \
  -cov 0.5 \
  --clusterAlg average \
  -p ${SLURM_CPUS_PER_TASK}

echo "=== Step 4: genome-wide summary for clustered compare (new) ==="

inStrain genome_wide \
  -i "${COMPARE_DIR}" \
  -s "${STB_FILE}" \
  -p ${SLURM_CPUS_PER_TASK}

echo "Done: profiles, stb, clustered compare, and genome-wide outputs written under ${BASE_DIR}"

