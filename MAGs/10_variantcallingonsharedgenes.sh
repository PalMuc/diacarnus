

set -euo pipefail

# Activate conda env that has tabix/bcftools
source ~/miniconda3/etc/profile.d/conda.sh   # change if your conda lives elsewhere
conda activate tabix

module load gnu/12
module load samtools
module load bcftools

BASE="/home/lvilleg/Diacarnus_0526"

FA_DIR="${BASE}/All_MAGS_46_fasta"
BAM_DIR="${BASE}/03_mapping_MAGs/dedup_q30_byMAG"
SHARED_DIR="${BASE}/07_MAGS_genes/genes_shared_shallow_meso"
VCF_DIR="${BASE}/08_shared_genes_variants"

mkdir -p "$VCF_DIR"

MAGS=(
  "Acidimicrobiales_1.1"
  "Arenicellales_2.1"
  "Latescibacterota_1.1"
  "Poribacteria_1.1"
  "Poribacteria_3.1"
  "Porisulfidus_sp_1.1"
  "Pseudomonadales_3.1"
)

for MAG in "${MAGS[@]}"; do
  REF="${FA_DIR}/odDiaEryt1.${MAG}.primary.renamed.fa"
  REG_BED="${SHARED_DIR}/${MAG}.genes_shared_shallow5_meso5.coords.bed"
  OUT_VCF_GZ="${VCF_DIR}/${MAG}.shared_genes.raw.vcf.gz"

  echo "=== MAG: $MAG ==="
  echo "REF:      $REF"
  echo "BED:      $REG_BED"
  echo "OUT VCF:  $OUT_VCF_GZ"

  if [[ ! -s "$REG_BED" ]]; then
    echo "[WARN] No shared-gene BED for this MAG, skipping."
    echo
    continue
  fi

  if [[ ! -s "$REF" ]]; then
    echo "[WARN] Reference FASTA not found: $REF, skipping."
    echo
    continue
  fi

  BAM_LIST=()
  for BAM in "${BAM_DIR}"/GW*.${MAG}.dedup.q30.bam; do
    [[ -e "$BAM" ]] || continue
    BAM_LIST+=("$BAM")
  done

  if [[ "${#BAM_LIST[@]}" -eq 0 ]]; then
    echo "[WARN] No BAMs for MAG $MAG, skipping."
    echo
    continue
  fi

  echo "  Using ${#BAM_LIST[@]} BAMs:"
  printf "    %s\n" "${BAM_LIST[@]}"

  if [[ ! -s "${REF}.fai" ]]; then
    echo "  Indexing reference..."
    samtools faidx "$REF"
  fi

  echo "  Calling variants on shared genes..."

  bcftools mpileup \
    -f "$REF" \
    -a DP,AD \
    -Q 20 \
    -q 0.01 \
    -b <(printf "%s\n" "${BAM_LIST[@]}") \
    -R "$REG_BED" \
  | bcftools call \
    -mv --ploidy 1 \
    -Oz \
    -o "$OUT_VCF_GZ"

  tabix -p vcf "$OUT_VCF_GZ"

  echo "  Finished MAG $MAG, VCF: $OUT_VCF_GZ"
  echo
done

echo "All MAGs processed."
