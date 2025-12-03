#!/bin/bash

# Script for principal component analyses (PCAs)

# --------------- PCA ON POPGEN WITH CLONES FILE FOR HOST-SYM GENOTYPING --------------- #

VCF_FILE="PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.recode.vcf"
OUTPUT_PREFIX="PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.PCA"

plink --vcf $VCF_FILE --double-id --allow-extra-chr --pca --set-missing-var-ids @:# --out $OUTPUT_PREFIX

# Convert PCA results to CSV for easier visualization
awk '{print $1","$2","$3","$4","$5","$6}' "${OUTPUT_PREFIX}.eigenvec" > "${OUTPUT_PREFIX}.csv"

# continue with /host_symbiont/01.host_symbiont_association.R

# ------------------ PCA ON POPGEN NO CLONES FILE ---------------------- #

VCF_FILE="PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf"
OUTPUT_PREFIX="PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.PCA"

plink --vcf $VCF_FILE --double-id --allow-extra-chr --pca --set-missing-var-ids @:# --out $OUTPUT_PREFIX

# Convert PCA results to CSV for easier visualization
awk '{print $1","$2","$3","$4","$5","$6}' "${OUTPUT_PREFIX}.eigenvec" > "${OUTPUT_PREFIX}.csv"

# continue with 02.PCA_plots.R

# -------------- RANDOM RESAMPING OF SNPs FOR MULTIPLE PCA ------------- # 

VCF_FILE="PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf.gz"
N_SNPS=10000
REPS=20

# Extract SNP IDs
bcftools query -f '%CHROM\t%POS' $VCF_FILE > all_snps.txt

for i in $(seq 1 $REPS); do
  shuf all_snps.txt | head -n $N_SNPS > snps_subset.$i.txt
  bcftools view -R snps_subset.$i.txt $VCF_FILE -Ov -o snps_subset.$i.vcf
  plink --vcf snps_subset.$i.vcf --double-id --allow-extra-chr --pca --set-missing-var-ids @:# --out snps_subset.$i.PCA
  awk '{print $1","$2","$3","$4","$5","$6}' snps_subset.$i.PCA.eigenvec > snps_subset.$i.PCA.csv
done

# continue with 02.PCA_plots.R
