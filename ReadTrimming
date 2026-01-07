module load gnu/12
module load samtools
source ~/.bashrc
conda activate ngs-core

cd /home/lvilleg/Diacarnus/01_rawreads

OUTBASE=/home/lvilleg/Diacarnus/02_trimmed_fastp
mkdir -p "$OUTBASE"

for SAMPLE_DIR in GW34*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "=== Processing sample: $SAMPLE ==="

    R1=$(ls ${SAMPLE_DIR}/*_1.fq.gz)
    R2=$(ls ${SAMPLE_DIR}/*_2.fq.gz)

    OUTDIR=${OUTBASE}/${SAMPLE}
    mkdir -p "$OUTDIR"

    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "${OUTDIR}/${SAMPLE}_R1.trimmed.fq.gz" \
      -O "${OUTDIR}/${SAMPLE}_R2.trimmed.fq.gz" \
      -w $SLURM_CPUS_PER_TASK \
      -q 20 \
      -u 30 \
      -l 36 \
      -h "${OUTDIR}/${SAMPLE}.fastp.html" \
      -j "${OUTDIR}/${SAMPLE}.fastp.json"
done

echo "All samples processed."

