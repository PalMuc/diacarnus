module load gnu/12
source ~/.bashrc
conda activate instrain

REF_MAGS=/home/lvilleg/Diacarnus/00_Genome/Diacarnus_symbionts_MAGs.uniq.fa
BAM_DIR=/home/lvilleg/Diacarnus/03_mapping_MAGs/dedup_q30
OUT_DIR=/home/lvilleg/Diacarnus/04_inStrain/MAG_profiles
mkdir -p "${OUT_DIR}"

for bam in ${BAM_DIR}/GW34*.MAGs.dedup.q30.bam; do
    sample=$(basename "${bam}" .MAGs.dedup.q30.bam)
    echo "Profiling ${sample} ..."
    inStrain profile "${bam}" "${REF_MAGS}" \
        -o "${OUT_DIR}/${sample}_MAGs_inStrain" \
        -p ${SLURM_CPUS_PER_TASK}
done


module load gnu/12
source ~/.bashrc
conda activate r-instrain

Rscript - <<'EOF'
library(dplyr)
library(readr)
library(stringr)

base_dir   <- "/home/lvilleg/Diacarnus/04_inStrain"
profile_dir <- file.path(base_dir, "MAG_profiles")

# Use one representative scaffold_info to get scaffold names
si_file <- file.path(profile_dir,
                     "GW34916_MAGs_inStrain",
                     "output",
                     "GW34916_MAGs_inStrain_scaffold_info.tsv")

si <- read_tsv(si_file, show_col_types = FALSE)

# Derive MAG ID from scaffold name, e.g. drop ".primary_contig_..."
stb <- si %>%
  transmute(
    scaffold = scaffold,
    MAG      = sub("\\.primary_contig_.*$", "", scaffold)
  ) %>%
  distinct()

out_stb <- file.path(base_dir, "MAGs.scaff2mag.stb")
write_tsv(stb, out_stb, col_names = FALSE)

cat("Wrote stb to:", out_stb, "\n")
EOF


module load gnu/12
source ~/.bashrc
conda activate instrain

BASE_DIR=/home/lvilleg/Diacarnus/04_inStrain
COMPARE_IN=${BASE_DIR}/MAG_compare/MAGs_compare      # existing compare project
STB_FILE=${BASE_DIR}/MAGs.scaff2mag.stb             # stb you just created

cd "${BASE_DIR}"

inStrain genome_wide \
  -i "${COMPARE_IN}" \
  -s "${STB_FILE}" \
  -p ${SLURM_CPUS_PER_TASK}


module load gnu/12
source ~/.bashrc
conda activate instrain

BASE_DIR=/home/lvilleg/Diacarnus/04_inStrain
PROFILE_DIR=${BASE_DIR}/MAG_profiles
OUT_DIR=${BASE_DIR}/MAG_compare/MAGs_compare_clustered
STB_FILE=${BASE_DIR}/MAGs.scaff2mag.stb

cd "${BASE_DIR}"

# Sanity check
echo "Number of inStrain profiles:"
ls ${PROFILE_DIR}/GW*_MAGs_inStrain | wc -l
echo "Profiles:"
ls ${PROFILE_DIR}/GW*_MAGs_inStrain

inStrain compare \
  -i ${PROFILE_DIR}/GW*_MAGs_inStrain \
  -o "${OUT_DIR}" \
  --stb "${STB_FILE}" \
  -ani 0.99 \
  -cov 0.5 \
  --clusterAlg average \
  -p ${SLURM_CPUS_PER_TASK}


module load gnu/12
source ~/.bashrc
conda activate r-instrain

cd /home/lvilleg/Diacarnus/04_inStrain/MAG_compare/MAGs_compare/output

cat << 'EOF' > mag_popANI_summary.R
library(dplyr)
library(readr)
library(stringr)

# 1) Load inStrain genome-wide compare output
comp <- read_tsv("MAGs_compare_genomeWide_compare.tsv", show_col_types = FALSE)

# 2) Derive SHORT sample IDs from BAM names
comp <- comp %>%
  mutate(
    s1 = str_remove(name1, "[.]MAGs.dedup.q30.bam$"),
    s2 = str_remove(name2, "[.]MAGs.dedup.q30.bam$")
  )

# 3) Sample metadata with SHORT IDs
meta <- tribble(
  ~sample,  ~population,
  "GW34916", "Shallow",
  "GW34917", "Shallow",
  "GW34918", "Shallow",
  "GW34919", "Shallow",
  "GW34920", "Shallow",
  "GW34921", "Shallow",
  "GW34922", "Shallow",
  "GW34923", "Shallow",
  "GW34924", "Shallow",
  "GW34925", "Shallow",
  "GW34926", "Shallow",
  "GW34927", "Shallow",
  "GW34928", "Mesophotic",
  "GW34929", "Mesophotic",
  "GW34930", "Mesophotic",
  "GW34931", "Mesophotic",
  "GW34932", "Mesophotic",
  "GW34933", "Mesophotic",
  "GW34934", "Mesophotic"
)
# 4) Join metadata to both samples
comp2 <- comp %>%
  left_join(meta, by = c("s1" = "sample")) %>%
  rename(pop1 = population) %>%
  left_join(meta, by = c("s2" = "sample")) %>%
  rename(pop2 = population)

# 5) Define pair class (within vs between pops) correctly
comp2 <- comp2 %>%
  mutate(
    pair_class = case_when(
      is.na(pop1) | is.na(pop2) ~ "Unknown",
      pop1 == pop2              ~ pop1,
      TRUE ~ paste(pmin(pop1, pop2), pmax(pop1, pop2), sep = "_vs_")
    )
  )

# 6) Summarise per MAG and pair_class
summary_mag_pop <- comp2 %>%
  group_by(genome, pair_class) %>%
  summarise(
    n_pairs       = n(),
    mean_popANI   = mean(popANI, na.rm = TRUE),
    median_popANI = median(popANI, na.rm = TRUE),
    min_popANI    = min(popANI, na.rm = TRUE),
    mean_conANI   = mean(conANI, na.rm = TRUE),
    median_conANI = median(conANI, na.rm = TRUE),
    min_conANI    = min(conANI, na.rm = TRUE),
    .groups       = "drop"
  )

write_tsv(summary_mag_pop, "MAG_population_popANI_conANI_summary.tsv")
EOF

Rscript mag_popANI_summary.R

