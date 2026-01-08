module load gnu/12
source ~/.bashrc
conda activate freebayes-parallel   # env with freebayes, vcflib, samtools, vcftools

BAM_DIR=/home/lvilleg/Diacarnus/03_mapping/dedup_q30
REF=/home/lvilleg/Diacarnus/00_Genome/Diacarnus_genomic.fa

RAW_VCF=${BAM_DIR}/host_freebayes.allsites.raw.vcf
FILT_PREFIX=${BAM_DIR}/host_freebayes.allsites.lowcov.filtered
FILT_VCF=${FILT_PREFIX}.vcf

cd "${BAM_DIR}"

# 1) List BAMs
ls GW34*.dedup.q30.bam > bamlist.txt

echo "Generating regions from FASTA..."
python /home/lvilleg/scripts/make_regions_from_fasta.py \
       "${REF}" 100000 > regions.txt

echo "Starting FreeBayes ALL-SITES calling in parallel on 20 cores..."
freebayes-parallel regions.txt 20 \
    -f "${REF}" \
    --bam-list bamlist.txt \
    --report-monomorphic \
    > "${RAW_VCF}"

echo "FreeBayes all-sites VCF: ${RAW_VCF}"
echo "Starting low-coverage friendly filtering (keeping invariant sites)..."

vcftools --vcf "${RAW_VCF}" \
         --min-meanDP 2 \
         --max-meanDP 30 \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out "${FILT_PREFIX}"

mv "${FILT_PREFIX}.recode.vcf" "${FILT_VCF}"

echo "Filtered all-sites VCF for pixy: ${FILT_VCF}"
echo "Done."

##Filtering the VCf so it works later on with pixy

module load gnu/12
source ~/.bashrc
conda activate ngs-core   # needs vcftools, bcftools, bgzip, tabix

WORKDIR=/home/lvilleg/Diacarnus/03_mapping/dedup_q30
cd "${WORKDIR}"

RAW_ALLSITES=${WORKDIR}/host_freebayes.allsites.raw.vcf.gz

# Make sure raw all-sites VCF is compressed and indexed
if [ ! -f "${RAW_ALLSITES}" ]; then
    echo "Bgzipping raw all-sites VCF..."
    bgzip -c host_freebayes.allsites.raw.vcf > "${RAW_ALLSITES}"
fi
if [ ! -f "${RAW_ALLSITES}.tbi" ]; then
    tabix -p vcf "${RAW_ALLSITES}"
fi

########################################
# 1) FILTER INVARIANT SITES ONLY
########################################
# Invariant = MAF = 0 (pixy doc suggestion).
# Use depth + missingness filters, but avoid minQ which can nuke invariants.
# Tune MIN_MEAN_DP / MAX_MEAN_DP to your 5x data; start conservative.

MIN_MEAN_DP_INVAR=3
MAX_MEAN_DP_INVAR=30
MAX_MISSING_INVAR=0.8

echo "Filtering invariant sites..."
vcftools --gzvcf "${RAW_ALLSITES}" \
         --max-maf 0 \
         --max-missing ${MAX_MISSING_INVAR} \
         --min-meanDP ${MIN_MEAN_DP_INVAR} \
         --max-meanDP ${MAX_MEAN_DP_INVAR} \
         --recode --stdout \
  | bgzip -c > host_freebayes.invariant.filtered.vcf.gz

tabix -p vcf host_freebayes.invariant.filtered.vcf.gz

########################################
# 2) FILTER VARIANT SITES ONLY
########################################
# Variant = MAC >= 1 (pixy doc suggestion).
# You can afford harsher depth/QUAL filters here.

MIN_MEAN_DP_VAR=5
MAX_MEAN_DP_VAR=30
MAX_MISSING_VAR=0.8
MINQ_VAR=30

echo "Filtering variant sites..."
vcftools --gzvcf "${RAW_ALLSITES}" \
         --mac 1 \
         --max-missing ${MAX_MISSING_VAR} \
         --min-meanDP ${MIN_MEAN_DP_VAR} \
         --max-meanDP ${MAX_MEAN_DP_VAR} \
         --minQ ${MINQ_VAR} \
         --recode --stdout \
  | bgzip -c > host_freebayes.variant.filtered.vcf.gz

tabix -p vcf host_freebayes.variant.filtered.vcf.gz

########################################
# 3) CONCATENATE INVARIANT + VARIANT
########################################
# This gives a single all-sites VCF with both classes retained.

echo "Concatenating invariant + variant filtered VCFs..."
bcftools concat \
  --allow-overlaps \
  host_freebayes.variant.filtered.vcf.gz \
  host_freebayes.invariant.filtered.vcf.gz \
  -Oz -o host_freebayes.allsites.filtered.pixy_ready.vcf.gz

tabix -p vcf host_freebayes.allsites.filtered.pixy_ready.vcf.gz

echo "Pixy-ready all-sites VCF: host_freebayes.allsites.filtered.pixy_ready.vcf.gz"
echo "Done."

