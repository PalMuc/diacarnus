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

      }
    ' "$BED" >> "$OUT"

    rm -f "$DEPTH_TMP"
  done

  echo "  Finished MAG $MAG, output: $OUT"
  echo
done

echo "All MAGs processed at $(date)"

#Then intersect those genes that are fully covered and present at the same time in at least 5 individuals of meso and 5 individuals of shallow

BASE="/home/lvilleg/Diacarnus_0526/07_MAGS_genes"

COV_DIR="${BASE}/gene_fullcov_1x_from_bakta_v2"
BED_DIR="${BASE}/genes_bed"
OUT_DIR="${BASE}/genes_shared_shallow_meso"

mkdir -p "$OUT_DIR"

# Full correct sample->population mapping file
POP_FILE="/home/lvilleg/Diacarnus_0526/sample_pops_shallow_meso_full.txt"

# Thresholds
MIN_SHALLOW=5
MIN_MESO=5

MAGS=(
  "Acidimicrobiales_1.1"
  "Arenicellales_2.1"
  "Latescibacterota_1.1"
  "Poribacteria_1.1"
  "Poribacteria_3.1"
  "Porisulfidus_sp_1.1"
  "Pseudomonadales_3.1"
)

echo "Coverage dir: $COV_DIR"
echo "BED dir:      $BED_DIR"
echo "Output dir:   $OUT_DIR"
echo "Pop file:     $POP_FILE"
echo "Thresholds:   Shallow >= $MIN_SHALLOW ; Mesophotic >= $MIN_MESO"
echo

if [[ ! -s "$POP_FILE" ]]; then
  echo "[ERROR] Population file missing or empty: $POP_FILE"
  exit 1
fi

for MAG in "${MAGS[@]}"; do
  COV="${COV_DIR}/${MAG}.gene_fullcov_1x.tsv"
  BED="${BED_DIR}/${MAG}.genes.bamnames.bed"

  echo "=== MAG: $MAG ==="
  echo "Coverage file: $COV"
  echo "BED file:      $BED"

  if [[ ! -s "$COV" ]]; then
    echo "[WARN] Missing or empty coverage file, skipping."
    echo
    continue
  fi

  if [[ ! -s "$BED" ]]; then
    echo "[WARN] Missing or empty BED, skipping."
    echo
    continue
  fi

  SHARED_GENES_TXT="${OUT_DIR}/${MAG}.genes_shared_shallow${MIN_SHALLOW}_meso${MIN_MESO}.txt"
  SHARED_COUNTS_TSV="${OUT_DIR}/${MAG}.genes_shared_shallow${MIN_SHALLOW}_meso${MIN_MESO}.counts.tsv"
  SHARED_BED="${OUT_DIR}/${MAG}.genes_shared_shallow${MIN_SHALLOW}_meso${MIN_MESO}.coords.bed"

  # Count per gene: #Shallow and #Mesophotic samples with full coverage
  awk -v OFS="\t" -v popfile="$POP_FILE" -v minS="$MIN_SHALLOW" -v minM="$MIN_MESO" '
    BEGIN{
      FS=OFS="\t"
      # load IID -> pop mapping
      while ((getline line < popfile) > 0) {
        if (line ~ /^IID([ \t]+|$)/) continue
        if (line ~ /^[ \t]*$/) continue

        split(line, f, /[ \t]+/)
        s = f[1]
        g = f[2]
        if (s != "" && g != "") {
          pop[s] = g
        }
      }
      close(popfile)
      # header for counts file
      print "gene_id","n_shallow","n_mesophotic" > "'"$SHARED_COUNTS_TSV"'"
    }

    NR==1 { next }  # skip coverage header

    {
      samp = $1
      gid  = $2

      if (!(samp in pop)) next

      key = samp SUBSEP gid
      if (key in seen) next
      seen[key] = 1

      if (pop[samp] == "Shallow") {
        cS[gid]++
      } else if (pop[samp] == "Mesophotic") {
        cM[gid]++
      }
    }

    END{
      # collect all genes seen in either group
      for (g in cS) all[g] = 1
      for (g in cM) all[g] = 1

      for (g in all) {
        s = (g in cS ? cS[g] : 0)
        m = (g in cM ? cM[g] : 0)
        if (s >= minS && m >= minM) {
          print g > "'"$SHARED_GENES_TXT"'"
          print g, s, m >> "'"$SHARED_COUNTS_TSV"'"
        }
      }
    }
  ' "$COV"

  # ensure file exists even if empty
  if [[ ! -f "$SHARED_GENES_TXT" ]]; then
    : > "$SHARED_GENES_TXT"
  fi

  sort -u "$SHARED_GENES_TXT" -o "$SHARED_GENES_TXT"

  NGENES=$(wc -l < "$SHARED_GENES_TXT")
  echo "  Genes with >=${MIN_SHALLOW} Shallow and >=${MIN_MESO} Mesophotic: $NGENES"

  if [[ "$NGENES" -eq 0 ]]; then
    echo "  No genes met thresholds for this MAG."
    echo
    continue
  fi

  # join back to BED for coordinates
  awk 'BEGIN{FS=OFS="\t"}
       FNR==NR {keep[$1]=1; next}
       ($4 in keep) {print}
  ' "$SHARED_GENES_TXT" "$BED" > "$SHARED_BED"

  echo "  Wrote:"
  echo "    Gene IDs : $SHARED_GENES_TXT"
  echo "    Counts   : $SHARED_COUNTS_TSV"
  echo "    BED      : $SHARED_BED"
  echo
done

echo "Done."
