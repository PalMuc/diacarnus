#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/home/lvilleg/Diacarnus_0526"
STB="${BASE_DIR}/04_inStrain_46MAGS/MAGs_new.scaff2mag.stb"
OUT_SCAFF="${BASE_DIR}/07_MAGS_genes/scaffolds"

mkdir -p "$OUT_SCAFF"

echo "Using STB: $STB"
echo "Writing scaffold lists to: $OUT_SCAFF"

# Each block: grep MAG ID, strip MAG ID from end of line to get scaffold name

# Acidimicrobiales_1.1
grep "odDiaEryt1.Acidimicrobiales_1.1" "$STB" \
  | sed 's/odDiaEryt1.Acidimicrobiales_1.1$//' \
  > "${OUT_SCAFF}/Acidimicrobiales_1.1.scaffolds.txt"

# Arenicellales_2.1
grep "odDiaEryt1.Arenicellales_2.1" "$STB" \
  | sed 's/odDiaEryt1.Arenicellales_2.1$//' \
  > "${OUT_SCAFF}/Arenicellales_2.1.scaffolds.txt"

# Latescibacterota_1.1
grep "odDiaEryt1.Latescibacterota_1.1" "$STB" \
  | sed 's/odDiaEryt1.Latescibacterota_1.1$//' \
  > "${OUT_SCAFF}/Latescibacterota_1.1.scaffolds.txt"

# Poribacteria_1.1
grep "odDiaEryt1.Poribacteria_1.1" "$STB" \
  | sed 's/odDiaEryt1.Poribacteria_1.1$//' \
  > "${OUT_SCAFF}/Poribacteria_1.1.scaffolds.txt"

# Poribacteria_3.1
grep "odDiaEryt1.Poribacteria_3.1" "$STB" \
  | sed 's/odDiaEryt1.Poribacteria_3.1$//' \
  > "${OUT_SCAFF}/Poribacteria_3.1.scaffolds.txt"

# Porisulfidus_sp_1.1
grep "odDiaEryt1.Porisulfidus_sp_1.1" "$STB" \
  | sed 's/odDiaEryt1.Porisulfidus_sp_1.1$//' \
  > "${OUT_SCAFF}/Porisulfidus_1.1.scaffolds.txt"

# Pseudomonadales_3.1
grep "odDiaEryt1.Pseudomonadales_3.1" "$STB" \
  | sed 's/odDiaEryt1.Pseudomonadales_3.1$//' \
  > "${OUT_SCAFF}/Pseudomonadales_3.1.scaffolds.txt"

echo "Line counts in scaffold files:"
wc -l "${OUT_SCAFF}"/*.scaffolds.txt

#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/home/lvilleg/Diacarnus_0526"
OUT_SCAFF="${BASE_DIR}/07_MAGS_genes/scaffolds"

# Rewrite each existing scaffold file into BED: contig, 0, big_end
for f in "${OUT_SCAFF}"/*.scaffolds.txt; do
  tmp="${f}.bedtmp"
  echo "Rewriting $(basename "$f") to BED format"
  awk 'NF>0{OFS="\t"; print $1, 0, 1000000000000}' "$f" > "$tmp"
  mv "$tmp" "$f"
done

echo "Preview Acidimicrobiales_1.1.scaffolds.txt:"
head -n 3 "${OUT_SCAFF}/Acidimicrobiales_1.1.scaffolds.txt"


module load gnu/12
module load samtools

BASE_DIR="/home/lvilleg/Diacarnus_0526"
BAM_DIR="${BASE_DIR}/03_mapping_MAGs/dedup_q30"
SCAFF_DIR="${BASE_DIR}/07_MAGS_genes/scaffolds"
OUT_DIR="${BASE_DIR}/03_mapping_MAGs/dedup_q30_byMAG"

mkdir -p "$OUT_DIR"

MAGS=("Acidimicrobiales_1.1" \
      "Arenicellales_2.1" \
      "Latescibacterota_1.1" \
      "Poribacteria_1.1" \
      "Poribacteria_3.1" \
      "Porisulfidus_1.1" \
      "Pseudomonadales_3.1")

echo "Input BAM dir:  $BAM_DIR"
echo "Scaffold lists: $SCAFF_DIR"
echo "Output BAM dir: $OUT_DIR"
echo

for bam in "${BAM_DIR}"/GW*.MAGs.dedup.q30.bam; do
  sample=$(basename "$bam" .MAGs.dedup.q30.bam)
  echo "Processing sample: $sample  ($(date))"

  for mag in "${MAGS[@]}"; do
    scaff_file="${SCAFF_DIR}/${mag}.scaffolds.txt"
    out_bam="${OUT_DIR}/${sample}.${mag}.dedup.q30.bam"

    if [ ! -s "$scaff_file" ]; then
      echo "  [WARN] Scaffold file missing or empty for MAG $mag: $scaff_file"
      continue
    fi

    echo "  Subsetting MAG $mag -> $(basename "$out_bam")"

    samtools view -@ "$SLURM_CPUS_PER_TASK" -b -L "$scaff_file" "$bam" > "$out_bam"
    samtools index "$out_bam"
  done

  echo
done

echo "Done at $(date)"
