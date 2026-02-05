
module load gnu/12
source ~/.bashrc
conda activate ngs-core

# paths
TRIM_DIR=/home/lvilleg/Diacarnus/02_trimmed_fastp
MAP_DIR=/home/lvilleg/Diacarnus/08_poribacteria
REF=/home/lvilleg/Diacarnus/00_Genome/odDiaEryt1.Poribacteria_2.1.primary.fa

# index Poribacteria MAG (only needs to be done once for this REF)
bwa-mem2 index "${REF}"   # bwa-mem2 uses the FASTA path/prefix here [web:26]

module load gnu/12
module load samtools

mkdir -p "$MAP_DIR"

cd "$TRIM_DIR"

for SAMPLE_DIR in GW34*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "=== Mapping sample: $SAMPLE ==="

    R1=${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R1.trimmed.fq.gz
    R2=${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R2.trimmed.fq.gz

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "  !! Missing R1 or R2 for $SAMPLE, skipping"
        continue
    fi

    OUT_BAM=${MAP_DIR}/${SAMPLE}.poribacteria.sorted.bam

    RG="@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA\\tPU:1"

    bwa-mem2 mem -t "$SLURM_CPUS_PER_TASK" -R "$RG" "$REF" "$R1" "$R2" \
        | samtools sort -@ 2 -o "$OUT_BAM" -   # samtools sort writes sorted BAM to -o file [web:31][web:35]

    samtools index "$OUT_BAM"

    echo "=== Done: $SAMPLE -> $OUT_BAM ==="
done

echo "All samples mapped."

# Get BUSCO completeness and gene coordinates from the base genome

module load gnu/12
source ~/.bashrc

# Activate BUSCO environment (created beforehand)
conda activate busco

# Paths
IN_FA=/home/lvilleg/Diacarnus/00_Genome/odDiaEryt1.Poribacteria_2.1.primary.fa
OUT_DIR=/home/lvilleg/Diacarnus/09_busco
DL_DIR=/home/lvilleg/Diacarnus/09_busco_downloads

mkdir -p "$OUT_DIR" "$DL_DIR"

cd "$OUT_DIR"

# Optional: pre-download lineage (once). Comment out after first run.
# busco --download bacteria_odb12 --download_path "$DL_DIR"

# Run BUSCO in genome mode on the Poribacteria MAG
busco \
  --in "$IN_FA" \
  --out poribacteria_busco_bacteria_odb12 \
  --mode genome \
  --lineage_dataset bacteria_odb12 \
  --download_path "$DL_DIR" \
  --cpu "$SLURM_CPUS_PER_TASK"

echo "BUSCO run finished."

