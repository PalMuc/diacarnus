#!/usr/bin/env bash
set -euo pipefail

# Run this from: /home/lvilleg/Diacarnus_0526

BASE_DIR="/home/lvilleg/Diacarnus_0526"
STB="$BASE_DIR/04_inStrain_46MAGS/MAGs_new.scaff2mag.stb"
ANN_DIR="$BASE_DIR/07_MAGS_genes/annotation"
OUT_DIR="$BASE_DIR/07_MAGS_genes"

cd "$BASE_DIR"

echo "Using scaffold->MAG mapping: $STB"
echo "Using annotations from:       $ANN_DIR"
echo

mkdir -p "$OUT_DIR/scaffolds" "$OUT_DIR/genes_bed"

# ---- STEP 1: scaffold lists per focal MAG ----
# Keys are short labels, values are the MAG IDs as they appear in the stb file
declare -A MAGS
MAGS["Acidimicrobiales_1.1"]="odDiaEryt1.Acidimicrobiales_1.1.primary"
MAGS["Arenicellales_2.1"]="odDiaEryt1.Arenicellales_2.1.primary"
MAGS["Latescibacterota_1.1"]="odDiaEryt1.Latescibacterota_1.1.primary"
MAGS["Poribacteria_1.1"]="odDiaEryt1.Poribacteria_1.1.primary"
MAGS["Poribacteria_3.1"]="odDiaEryt1.Poribacteria_3.1.primary"
MAGS["Porisulfidus_1.1"]="odDiaEryt1.Porisulfidus_sp_1.1.primary"
MAGS["Porisulfidus_1.1_rep2"]="odDiaEryt1.Porisulfidus_sp_1.1.primary.rep2"

echo "== Step 1: extracting scaffold lists per MAG =="

for key in "${!MAGS[@]}"; do
    mag_id="${MAGS[$key]}"
    out_scaff="$OUT_DIR/scaffolds/${key}.scaffolds.txt"
    echo "  MAG $key  (ID: $mag_id) -> $out_scaff"
    # Assuming stb format: scaffold<TAB>MAG_ID
    grep -F "$mag_id" "$STB" | awk '{print $1}' > "$out_scaff"
    echo "    $(wc -l < "$out_scaff") scaffolds written"
done

echo
echo "== Step 2: generating gene BED files from Bakta GFFs =="

# Function to build BED from a given GFF
build_bed () {
    local gff="$1"
    local label="$2"
    local out_bed="$OUT_DIR/genes_bed/${label}.genes.bed"

    echo "  GFF: $gff"
    echo "  BED: $out_bed"

    grep -v '^#' "$gff" \
      | awk '$3=="CDS"{
               OFS="\t";
               split($9,a,";");
               id=".";
               for(i in a){
                 if(a[i] ~ /^ID=/){
                   sub(/^ID=/,"",a[i]);
                   id=a[i]
                 }
               }
               print $1, $4-1, $5, id
             }' > "$out_bed"

    echo "    $(wc -l < "$out_bed") CDS entries written"
}

build_bed "$ANN_DIR/odDiaEryt1.Acidimicrobiales_1.1.primary.bakta.gff3"       "Acidimicrobiales_1.1"
build_bed "$ANN_DIR/odDiaEryt1.Arenicellales_2.1.primary.bakta.gff3"           "Arenicellales_2.1"
build_bed "$ANN_DIR/odDiaEryt1.Latescibacterota_1.1.primary.bakta.gff3"       "Latescibacterota_1.1"
build_bed "$ANN_DIR/odDiaEryt1.Poribacteria_1.1.primary.bakta.gff3"           "Poribacteria_1.1"
build_bed "$ANN_DIR/odDiaEryt1.Poribacteria_3.1.primary.bakta.gff3"           "Poribacteria_3.1"
build_bed "$ANN_DIR/odDiaEryt1.Porisulfidus_sp_1.1.primary.bakta.gff3"        "Porisulfidus_1.1"
build_bed "$ANN_DIR/odDiaEryt1.Porisulfidus_sp_1.1.primary.rep2.bakta.gff3"   "Porisulfidus_1.1_rep2"

echo
echo "Done."
echo "  Scaffold lists in: $OUT_DIR/scaffolds/"
echo "  Gene BED files in: $OUT_DIR/genes_bed/"
