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

#Then do a coverage check, check all the annotated genes and wether they have at least 1X coverage to continue
module load gnu/12
module load samtools

BASE_DIR="/home/lvilleg/Diacarnus_0526"
BAM_DIR="${BASE_DIR}/03_mapping_MAGs/dedup_q30_byMAG"
BED_DIR="${BASE_DIR}/07_MAGS_genes/genes_bed"
OUT_DIR="${BASE_DIR}/07_MAGS_genes/gene_fullcov_1x_from_bakta_v2"

mkdir -p "$OUT_DIR"

MAGS=("Acidimicrobiales_1.1" \
      "Arenicellales_2.1" \
      "Latescibacterota_1.1" \
      "Poribacteria_1.1" \
      "Poribacteria_3.1" \
      "Porisulfidus_sp_1.1" \
      "Pseudomonadales_3.1")

echo "BAM dir:  $BAM_DIR"
echo "BED dir:  $BED_DIR"
echo "Out dir:  $OUT_DIR"
echo

for MAG in "${MAGS[@]}"; do
  BED="${BED_DIR}/${MAG}.genes.bamnames.bed"

  if [ ! -s "$BED" ]; then
    echo "[WARN] Missing or empty BED for MAG $MAG: $BED"
    continue
  fi

  echo "=== MAG: $MAG ==="
  echo "Using BED: $BED"

  OUT="${OUT_DIR}/${MAG}.gene_fullcov_1x.tsv"
  echo -e "sample\tgene_id\tcovered_bases\tlength\tfrac_covered" > "$OUT"

  for BAM in "${BAM_DIR}"/GW*."${MAG}".dedup.q30.bam; do
    [ -e "$BAM" ] || continue
    SAMPLE=$(basename "$BAM" .${MAG}.dedup.q30.bam)

    echo "  Processing sample: $SAMPLE ($(date))"

    DEPTH_TMP=$(mktemp)

    # All depth within all CDS of this MAG for this sample
    samtools depth -a -b "$BED" "$BAM" > "$DEPTH_TMP"

 awk -F'\t' -v OFS="\t" -v sample="$SAMPLE" '
      BEGIN{
        # cd to 1-based internally, map pos ^f^rdepth for current contig
      }
      FNR==NR { next }  # (dummy; we only read one file in this awk)

      # This awk will actually be run with BED via `awk ... BED DEPTH_TMP`
    ' /dev/null "$DEPTH_TMP" > /dev/null

    # Now do the gene-wise scan: loop over BED, query DEPTH_TMP per gene
    # (one awk pass that reads depth into per-contig arrays and checks genes)

    awk -v OFS="\t" -v sample="$SAMPLE" -v depthfile="$DEPTH_TMP" '
      BEGIN{
        FS=OFS="\t"
        # load depth per contig into associative arrays: depth[contig, pos] = depth
        while ((getline line < depthfile) > 0) {
          split(line, a, "\t")
          c   = a[1]
          pos = a[2]
          d   = a[3]
          key = c SUBSEP pos
          depth[key] = d
        }
        close(depthfile)
      }
      {
        # BED: contig, start(0-based), end, gene_id
        contig = $1
        start0 = $2
        end    = $3
        gid    = $4

        len = end - start0
        covered = 0

        # scan all positions in this gene
        for (p = start0 + 1; p <= end; p++) {
          key = contig SUBSEP p
          d = (key in depth ? depth[key] : 0)
          if (d >= 1) covered++
        }

        if (len > 0) {
          frac = covered / len
        } else {
          frac = 0
        }
        # strict: all bases >=1x
        if (covered == len && len > 0) {
          print sample, gid, covered, len, frac
        }

        # to relax to 90% coverage instead, comment the if() above and use:
        # if (len > 0 && frac >= 0.9) { print sample, gid, covered, len, frac }
      }
    ' "$BED" >> "$OUT"

    rm -f "$DEPTH_TMP"
  done

  echo "  Finished MAG $MAG, output: $OUT"
  echo
done

echo "All MAGs processed at $(date)"


