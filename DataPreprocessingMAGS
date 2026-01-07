# paths
TRIM_DIR=/home/lvilleg/Diacarnus/02_trimmed_fastp
MAP_DIR=/home/lvilleg/Diacarnus/03_mapping_MAGs
REF=/home/lvilleg/Diacarnus/00_Genome/Diacarnus_symbionts_MAGs.uniq.fa

bwa-mem2 index "${REF}"

module load gnu/12
module load samtools

mkdir -p "$MAP_DIR"

# Optional: sanity check that reference + index exist
if [[ ! -f "${REF}" ]]; then
    echo "ERROR: reference FASTA $REF not found"
    exit 1
fi

cd "$TRIM_DIR"

for SAMPLE_DIR in GW34*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "=== Mapping sample to symbiont MAGs: $SAMPLE ==="

    R1=${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R1.trimmed.fq.gz
    R2=${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R2.trimmed.fq.gz

    # sanity check
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "  !! Missing R1 or R2 for $SAMPLE, skipping"
        continue
    fi

    OUT_BAM=${MAP_DIR}/${SAMPLE}.MAGs.sorted.bam

    # Read group
    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tPU:1"

    # Map to MAG-only reference, sort, index
    bwa-mem2 mem -t "$SLURM_CPUS_PER_TASK" -R "$RG" "$REF" "$R1" "$R2" \
        | samtools sort -@ 2 -o "$OUT_BAM" -

    samtools index "$OUT_BAM"

    echo "=== Done: $SAMPLE -> $OUT_BAM ==="
done

echo "All samples mapped to symbiont MAGs."


module load gnu/12
module load samtools

MAG_DIR=/home/lvilleg/Diacarnus/03_mapping_MAGs
OUT_DIR=${MAG_DIR}/flagstat

mkdir -p "${OUT_DIR}"

for bam in ${MAG_DIR}/GW34*.MAGs.sorted.bam; do
    [[ "${bam}" == *".tmp."* ]] && continue

    base=$(basename "${bam}" .MAGs.sorted.bam)
    echo "Running samtools flagstat on ${bam} ..."
    samtools flagstat -@ ${SLURM_CPUS_PER_TASK} "${bam}" > "${OUT_DIR}/${base}.flagstat.txt"
done

module load gnu/12
module load samtools

IN_DIR=/home/lvilleg/Diacarnus/03_mapping_MAGs/

OUT_DIR=${IN_DIR}/dedup_q30
TMP_DIR=${IN_DIR}/tmp_markdup

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

for bam in ${IN_DIR}/GW34*.sorted.bam; do
    base=$(basename "${bam}" .sorted.bam)
    echo "Processing ${bam} ..."

    # 1) Name-sort for fixmate
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -n \
        -o "${TMP_DIR}/${base}.namesort.bam" "${bam}"

    # 2) fixmate adds ms/MC tags
    samtools fixmate -@ ${SLURM_CPUS_PER_TASK} -m \
        "${TMP_DIR}/${base}.namesort.bam" \
        "${TMP_DIR}/${base}.fixmate.bam"

    # 3) Coordinate-sort (typo fixed: removed extra })
    samtools sort -@ ${SLURM_CPUS_PER_TASK} \
        -o "${TMP_DIR}/${base}.fixmate.possort.bam" \
        "${TMP_DIR}/${base}.fixmate.bam"

    # 4) Mark and remove duplicates
    samtools markdup -@ ${SLURM_CPUS_PER_TASK} -r \
        "${TMP_DIR}/${base}.fixmate.possort.bam" \
        "${OUT_DIR}/${base}.dedup.bam"

    # 5) Filter mapped, MAPQ ≥30
    samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 4 -q 30 \
        "${OUT_DIR}/${base}.dedup.bam" \
        > "${OUT_DIR}/${base}.dedup.q30.bam"

    # 6) Index
    samtools index -@ ${SLURM_CPUS_PER_TASK} "${OUT_DIR}/${base}.dedup.q30.bam"

    echo "Done ${base}"
done

