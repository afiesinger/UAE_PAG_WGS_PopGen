#!/bin/bash

# Script for kinship inference, clone mate identification and removal
# VCF files generated here are needed for downstream use

# ---------------- MAKE POPULATION FILES FOR DOWNSTREAM INCLUDING THE CLONE MATES ---------------- #
# bcftools v1.20

# with all samples (one version with and one without header -- for different programs)
bcftools query -l PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.recode.vcf  | awk 'BEGIN{OFS="\t"; print "SampleID","Population"} {split($0,a,"_"); print $0, a[2]}' > PHAR.popfile.all
tail -n +2 PHAR.popfile.all > PHAR.popfile.all.noheader

# PAG
awk '$2=="SA" || $2=="SY" {print $1}' PHAR.popfile.all.noheader > PAG.all.inds

# GO
awk '$2=="SI" {print $1}' PHAR.popfile.all.noheader > GO.all.inds

# -------------- KINSHIP INFERENCE FOR CLONE MATE IDENTIFICATION --------------- #
# PLINK v1.90
# KING v2.2.7

indir="/path/to/kinship"
cd $indir

# make input file for KING (need a .bed file as input)
plink --vcf POPGEN/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.recode.vcf --double-id --allow-extra-chr --make-bed --set-missing-var-ids @:# --out PHAR.OUT.SNPs_filt_kin.ld_pruned

# recode chromosome names all to 1 (i.e., all contigs would lie on CHR1 -- since we do not have chromosome information it is just an approximation for the program)
awk '{$1="1";print $0}' PHAR.OUT.SNPs_filt_kin.ld_pruned.bim > PHAR.OUT.SNPs_filt_kin.ld_pruned.bim.tmp
mv PHAR.OUT.SNPs_filt_kin.ld_pruned.bim.tmp PHAR.OUT.SNPs_filt_kin.ld_pruned.bim

king -b PHAR.OUT.SNPs_filt_kin.ld_pruned.bed --kinship

mv king.kin0 king.kin0_ALL

# ----------------- IBS DENDROGRAM WITH CLONES --------------- #
# PLINK v1.90

# this is for later use in 08.host_symbiont_association.R
plink --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --distance square ibs --out ibs_results

# -------------------  REMOVAL OF CLONES ------------------- #

# keep the file with least amount of missing data for each clone mate

popdir="/path/to/files"
cd $popdir

vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.recode.vcf --remove-indv UAE_SA_Phar_31 --recode --recode-INFO-all --out temp1
vcftools --vcf temp1.recode.vcf --remove-indv UAE_SI_Phar_38 --recode --recode-INFO-all --out PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES

mv PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.recode.vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf
rm temp*

bgzip -k PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf
tabix PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf.gz

# ----------- MAKE POPULATION FILES WITHOUT CLONES ----------- #

bcftools query -l PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf.gz  | awk 'BEGIN{OFS="\t"; print "SampleID","Population"} {split($0,a,"_"); print $0, a[2]}' > PHAR.popfile.noclones
tail -n +2 PHAR.popfile.noclones > PHAR.popfile.noclones.noheader

# PAG
awk '$2=="SA" || $2=="SY" {print $1}' PHAR.popfile.noclones.noheader > PAG.noclones.inds

# GO
awk '$2=="SI" {print $1}' PHAR.popfile.noclones.noheader > GO.noclones.inds

# SA
awk '$2=="SA" {print $1}' PHAR.popfile.noclones.noheader > SA.noclones.inds

# SY
awk '$2=="SY" {print $1}' PHAR.popfile.noclones.noheader > SY.noclones.inds

# SI
awk '$2=="SI" {print $1}' PHAR.popfile.noclones.noheader > SI.noclones.inds

# ----------- SPLIT VCF FILE BY POPULATION FOR FST CALCULATION ---------- #

# PAG and GO VCFs
vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf --keep PAG.noclones.inds --recode --recode-INFO-all --out PAG_noclones
vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf --keep GO.noclones.inds --recode --recode-INFO-all --out GO_noclones

# SA SY SI VCFs
vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf --keep SA.noclones.inds --recode --recode-INFO-all --out SA_noclones
vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf --keep SY.noclones.inds --recode --recode-INFO-all --out SY_noclones
vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf --keep SI.noclones.inds --recode --recode-INFO-all --out SI_noclones

# --------- IBS DENDROGRAM WITHOUT CLONES -------- #
# PLINK v1.90

# to verify the relatedness of remaining individuals

plink --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --distance square ibs --out ibs_results_noclones
