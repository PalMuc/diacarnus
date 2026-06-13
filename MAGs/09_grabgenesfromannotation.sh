#!/bin/bash
set -euo pipefail

BASE="/home/lvilleg/Diacarnus_0526/07_MAGS_genes"
ANN_DIR="${BASE}/annotation"
BED_DIR="${BASE}/genes_bed"

mkdir -p "$BED_DIR"

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
  GFF="${ANN_DIR}/odDiaEryt1.${MAG}.primary.bakta.gff3"
  BED="${BED_DIR}/${MAG}.genes.bamnames.bed"

  echo "=== ${MAG} ==="
  echo "GFF: ${GFF}"
  echo "BED: ${BED}"

  if [[ ! -s "$GFF" ]]; then
    echo "[WARN] Missing or empty GFF: $GFF"
    continue
  fi

  awk 'BEGIN{FS=OFS="\t"}
    $0 !~ /^#/ && $3=="CDS" {
      contig_num = $1
      sub(/^contig_/, "", contig_num)

      chrom = "odDiaEryt1.'"${MAG}"'.primary_contig_" contig_num
      start = $4 - 1
      end = $5

      id = ""
      locus = ""
      n = split($9, a, ";")
      for (i = 1; i <= n; i++) {
        split(a[i], b, "=")
        if (b[1] == "ID") id = b[2]
        if (b[1] == "locus_tag") locus = b[2]
      }

      gene_id = (locus != "" ? locus : id)
      if (gene_id == "") gene_id = "gene_" NR

      print chrom, start, end, gene_id
    }' "$GFF" > "$BED"

  echo "Wrote $(wc -l < "$BED") CDS entries"
  head -n 3 "$BED"
  echo
done

echo "Done."
