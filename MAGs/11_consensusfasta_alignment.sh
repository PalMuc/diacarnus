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


#Then create consensus sequences per gene for each of the individuals within a MAG

module load gnu/12
module load samtools
module load bcftools

BASE="/home/lvilleg/Diacarnus_0526"

FA_DIR="${BASE}/All_MAGS_46_fasta"
VCF_DIR="${BASE}/08_shared_genes_variants"
SHARED_DIR="${BASE}/07_MAGS_genes/genes_shared_shallow_meso"
COV_DIR="${BASE}/07_MAGS_genes/gene_fullcov_1x_from_bakta_v2"
CONS_DIR="${BASE}/09_gene_consensus_fullcov"

mkdir -p "$CONS_DIR"

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
  VCF_GZ="${VCF_DIR}/${MAG}.shared_genes.raw.vcf.gz"
  REG_BED="${SHARED_DIR}/${MAG}.genes_shared_shallow5_meso5.coords.bed"
  COV_TSV="${COV_DIR}/${MAG}.gene_fullcov_1x.tsv"
  OUT_FASTA="${CONS_DIR}/${MAG}.genes_shared_shallow5_meso5.fullcov.consensus.fasta"

  echo "=== MAG: $MAG ==="
  echo "REF:        $REF"
  echo "VCF:        $VCF_GZ"
  echo "GENE BED:   $REG_BED"
  echo "COV TSV:    $COV_TSV"
  echo "OUT FASTA:  $OUT_FASTA"

  if [[ ! -s "$REF" || ! -s "$VCF_GZ" || ! -s "$REG_BED" || ! -s "$COV_TSV" ]]; then
    echo "[WARN] Missing REF/VCF/BED/COV, skipping."
    echo
    continue
  fi
if [[ ! -s "${VCF_GZ}.tbi" && ! -s "${VCF_GZ}.csi" ]]; then
    echo "  Indexing VCF..."
    tabix -p vcf "$VCF_GZ"
  fi

  if [[ ! -s "${REF}.fai" ]]; then
    echo "  Indexing reference..."
    samtools faidx "$REF"
  fi

  : > "$OUT_FASTA"

  echo "  Pre-loading fully covered sample-gene pairs..."

  # Build a set of (sample,gene) pairs with full coverage
  # Assumes COV_TSV has header; sample in col1, gene_id in col2, covered in col3, length in col4
  FULLCOV_TMP=$(mktemp)
  awk 'NR>1 && $3 == $4 {print $1 "|" $2}' "$COV_TSV" > "$FULLCOV_TMP"

  echo "  Building per-sample consensus for each shared gene..."

  # BED: contig  start(0-based)  end  gene_id
  while read -r CONTIG START END GENE; do
    [[ -z "$CONTIG" ]] && continue

    POS_START=$((START + 1))
    POS_END=$END
    REGION="${CONTIG}:${POS_START}-${POS_END}"

    echo "    Gene $GENE region $REGION"

    TMP_VCF=$(mktemp)
    bcftools view -r "$REGION" "$VCF_GZ" -Oz -o "$TMP_VCF"
    tabix -p vcf "$TMP_VCF"

    # Loop over samples, but only keep those with full coverage for this gene
    bcftools query -l "$TMP_VCF" | while read -r SAMPLE; do
      KEY="${SAMPLE}|${GENE}"
      if ! grep -q "^${KEY}$" "$FULLCOV_TMP"; then
        continue
      fi
    TMP_FA=$(mktemp)

      samtools faidx "$REF" "$REGION" \
        | bcftools consensus \
            -s "$SAMPLE" \
            "$TMP_VCF" > "$TMP_FA"

      SEQ=$(awk 'NR>1 {printf "%s", $0} END{printf "\n"}' "$TMP_FA")

      echo ">${SAMPLE}|${GENE}" >> "$OUT_FASTA"
      echo "$SEQ" >> "$OUT_FASTA"

      rm -f "$TMP_FA"
    done

    rm -f "$TMP_VCF" "$TMP_VCF.tbi" "$TMP_VCF.csi"
  done < "$REG_BED"

  rm -f "$FULLCOV_TMP"

  echo "  Finished full-coverage consensus for MAG $MAG"
  echo
done

echo "All MAGs processed."

#Then split to get one fasta per gene of each of the MAGS

BASE="/home/lvilleg/Diacarnus_0526"

CONS_DIR="${BASE}/09_gene_consensus_fullcov"
SPLIT_DIR="${BASE}/10_gene_fastas_per_MAG"

mkdir -p "$SPLIT_DIR"

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
  CONS_FASTA="${CONS_DIR}/${MAG}.genes_shared_shallow5_meso5.fullcov.consensus.fasta"
  MAG_SPLIT_DIR="${SPLIT_DIR}/${MAG}"

  echo "=== Splitting by gene for MAG: $MAG ==="
  echo "  Consensus FASTA: $CONS_FASTA"
  echo "  Output dir:      $MAG_SPLIT_DIR"

  if [[ ! -s "$CONS_FASTA" ]]; then
    echo "  [WARN] Missing consensus FASTA, skipping $MAG"
    echo
    continue
  fi

  mkdir -p "$MAG_SPLIT_DIR"

  awk -v outdir="$MAG_SPLIT_DIR" '
    /^>/ {
      header = substr($0, 2)
      split(header, a, "|")
      sample = a[1]
      gene   = a[2]
      file = outdir "/" gene ".fasta"
      print ">" sample "|" gene >> file
      next
    }
 /^[^>]/ {
      print >> file
    }
  ' "$CONS_FASTA"

  echo "  Done for MAG $MAG"
  echo
done

echo "All MAGs split by gene."

#ANd finally use MAFFT to get alignments of each gene with all individuals and concatenate it 
BASE="/home/lvilleg/Diacarnus_0526"

CONS_DIR="${BASE}/09_gene_consensus_fullcov"
SPLIT_DIR="${BASE}/10_gene_fastas_per_MAG"
ALIGN_DIR="${BASE}/11_gene_alignments_per_MAG"
CONCAT_DIR="${BASE}/12_concatenated_alignments"

mkdir -p "$SPLIT_DIR" "$ALIGN_DIR" "$CONCAT_DIR"

# MAGs to process
MAGS=(
  "Acidimicrobiales_1.1"
  "Arenicellales_2.1"
  "Latescibacterota_1.1"
  "Poribacteria_1.1"
  "Poribacteria_3.1"
  "Porisulfidus_sp_1.1"
  "Pseudomonadales_3.1"
)

############################################
# 1. Split consensus by gene (per MAG)
############################################

for MAG in "${MAGS[@]}"; do
  CONS_FASTA="${CONS_DIR}/${MAG}.genes_shared_shallow5_meso5.fullcov.consensus.fasta"
  MAG_SPLIT_DIR="${SPLIT_DIR}/${MAG}"

  echo "=== Splitting by gene for MAG: $MAG ==="
  echo "  Consensus FASTA: $CONS_FASTA"
  echo "  Output dir:      $MAG_SPLIT_DIR"

  if [[ ! -s "$CONS_FASTA" ]]; then
    echo "  [WARN] Missing consensus FASTA, skipping MAG $MAG."
    echo
    continue
  fi

  mkdir -p "$MAG_SPLIT_DIR"
 # Headers are >Sample|Gene, sequences are plain A/C/G/T/N
  awk -v outdir="$MAG_SPLIT_DIR" '
    /^>/ {
      header = substr($0, 2)
      split(header, a, "|")
      sample = a[1]
      gene   = a[2]
      file = outdir "/" gene ".fasta"
      print ">" sample "|" gene >> file
      next
    }
    /^[^>]/ {
      print >> file
    }
  ' "$CONS_FASTA"

  echo "  Finished splitting genes for MAG $MAG"
  echo
done

############################################
# 2. Align each gene with MAFFT
############################################

for MAG in "${MAGS[@]}"; do
  MAG_SPLIT_DIR="${SPLIT_DIR}/${MAG}"
  MAG_ALIGN_DIR="${ALIGN_DIR}/${MAG}"

  if [[ ! -d "$MAG_SPLIT_DIR" ]]; then
    echo "[WARN] No split FASTAs for MAG $MAG, skipping alignment."
    continue
  fi

  mkdir -p "$MAG_ALIGN_DIR"

  echo "=== Aligning genes for MAG: $MAG ==="

  for GENE_FASTA in "${MAG_SPLIT_DIR}"/*.fasta; do
    [[ ! -e "$GENE_FASTA" ]] && continue

    GENE_NAME=$(basename "$GENE_FASTA" .fasta)
    OUT_ALN="${MAG_ALIGN_DIR}/${GENE_NAME}.aln.fasta"

    if [[ -s "$OUT_ALN" ]]; then
      echo "  [SKIP] Alignment exists for gene $GENE_NAME"
      continue
    fi
    mafft --auto --thread "$SLURM_CPUS_PER_TASK" "$GENE_FASTA" > "$OUT_ALN"
  done

  echo "  Finished alignments for MAG $MAG"
  echo
done

############################################
# 3. Concatenate all gene alignments (per MAG)
############################################

for MAG in "${MAGS[@]}"; do
  MAG_ALIGN_DIR="${ALIGN_DIR}/${MAG}"
  OUT_CONCAT="${CONCAT_DIR}/${MAG}.concatenated_genes.aln.fasta"

  if [[ ! -d "$MAG_ALIGN_DIR" ]]; then
    echo "[WARN] No alignments for MAG $MAG, skipping concatenation."
    continue
  fi

  echo "=== Concatenating alignments for MAG: $MAG ==="

  ALIGN_FILES=($(ls -1 "${MAG_ALIGN_DIR}"/*.aln.fasta 2>/dev/null || true))
  if [[ ${#ALIGN_FILES[@]} -eq 0 ]]; then
    echo "  [WARN] No .aln.fasta files for MAG $MAG."
    echo
    continue
  fi

  # 3a. List all samples present in any gene alignment
  TMP_SAMPLES=$(mktemp)
  for ALN in "${ALIGN_FILES[@]}"; do
    grep "^>" "$ALN" | sed 's/^>//' | cut -d"|" -f1
  done | sort -u > "$TMP_SAMPLES"

  : > "$OUT_CONCAT"

  # 3b. Concatenate for each sample
  while read -r SAMPLE; do
    echo ">${SAMPLE}" >> "$OUT_CONCAT"
    CONCAT_SEQ=""

    for ALN in "${ALIGN_FILES[@]}"; do
      SEQ=$(awk -v s="$SAMPLE" '
        BEGIN {found = 0}
        /^>/ {
          header = substr($0, 2)
          split(header, a, "|")
          if (a[1] == s) {
            found = 1
          } else {
            found = 0
          }
          next
        }
        found {
          gsub(/\r/, "", $0)
          printf "%s", $0
        }
        END { printf "\n" }
      ' "$ALN")

      if [[ -z "$SEQ" || "$SEQ" == "" ]]; then
        # Sample missing for this gene: fill with gaps of correct length
        GENE_LEN=$(awk '
          /^[^>]/ {gsub(/\r/, "", $0); l+=length($0)}
          END {print l}
        ' "$ALN")
        CONCAT_SEQ="${CONCAT_SEQ}$(perl -e "print '-' x $GENE_LEN")"
      else
        SEQ_CLEAN=$(echo "$SEQ" | tr -d '\n')
        CONCAT_SEQ="${CONCAT_SEQ}${SEQ_CLEAN}"
      fi
    done

    echo "$CONCAT_SEQ" >> "$OUT_CONCAT"
  done < "$TMP_SAMPLES"

  rm -f "$TMP_SAMPLES"

  echo "  Wrote concatenated alignment: $OUT_CONCAT"
  echo
done

echo "All MAGs processed for concatenated alignments."

