#Run BUSCO to get ortholohue coordinates on the base genome

module load gnu/12
source ~/.bashrc
conda activate busco

IN_FA=/home/lvilleg/Diacarnus/00_Genome/Diacarnus_genomic.fa
OUT_DIR=/home/lvilleg/Diacarnus/09_busco_Diacarnus
DL_DIR=/home/lvilleg/Diacarnus/09_busco_downloads

mkdir -p "$OUT_DIR" "$DL_DIR"
cd "$OUT_DIR"

# Optional: only first time
# busco --download metazoa_odb12 --download_path "$DL_DIR"

busco \
  --in "$IN_FA" \
  --out Diacarnus_busco_metazoa_odb12 \
  --mode genome \
  --lineage_dataset metazoa_odb12 \
  --download_path "$DL_DIR" \
  --cpu "$SLURM_CPUS_PER_TASK"

echo "BUSCO Diacarnus (metazoa_odb12) finished."

#Then on R match the previouls identified outliers from script 08 to the genes 

