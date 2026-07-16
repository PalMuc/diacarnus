 2) LD pruning on the PCA-appropriate marker set   #
#    (same window/step/r2 as the ADMIXTURE pruning: #
#    --indep-pairwise 50 5 0.2)                      #
#####################################################
plink \
  --bfile "${BFILE_MIND}" \
  --allow-extra-chr \
  --indep-pairwise 50 5 0.2 \
  --out "${OUTDIR}/${PREFIX}_pca_ldprune"

plink \
  --bfile "${BFILE_MIND}" \
  --allow-extra-chr \
  --extract "${OUTDIR}/${PREFIX}_pca_ldprune.prune.in" \
  --make-bed \
  --out "${OUTDIR}/${PREFIX}_pca_pruned"

#####################################################
# 3) PCA on the LD-pruned dataset (20 components)   #
#####################################################
plink \
  --bfile "${OUTDIR}/${PREFIX}_pca_pruned" \
  --allow-extra-chr \
  --pca 20 \
  --out "${OUTDIR}/${PREFIX}_pca_pruned_mind"

#####################################################
# 4) KING-robust relatedness (PLINK2, unpruned      #
#    mind-filtered set — pruning is not recommended #
#    for kinship estimation, unlike for PCA)        #
#####################################################
./software/plink2  \
  --bfile "${BFILE_MIND}" \
  --make-king-table \
  --out "${OUTDIR}/${PREFIX}_king"
