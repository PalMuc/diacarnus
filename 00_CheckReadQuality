module load gnu/12
source ~/.bashrc
conda activate qc_env   # env containing fastqc + multiqc

cd /home/lvilleg/Diacarnus/01_rawreads

OUTDIR=fastqc_raw
mkdir -p "${OUTDIR}"

# 1) Run FastQC on all raw reads
fastqc -t ${SLURM_CPUS_PER_TASK} -o "${OUTDIR}" \
  $(find . -type f -name "*.fq.gz")

# 2) Run MultiQC summarizing all FastQC reports
cd "${OUTDIR}"
multiqc . --filename multiqc_raw_fastqc_report


