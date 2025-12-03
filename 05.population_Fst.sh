#!/bin/bash

# Script to calculate population-wide Weir & Cockerham's Fst

# ------------------- POPULATION-WIDE Fst ---------------------- #

indir="/path/to/VCF"
outdir="/path/to/output"

cd $outdir

# -- SA vs SI -- #
vcftools --vcf $indir/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf \
    --weir-fst-pop $indir/SA.noclones.inds \
    --weir-fst-pop $indir/SI.noclones.inds \
    --out PHAR_Fst_SA_SI

# -- SA vs SY -- #
vcftools --vcf $indir/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf \
    --weir-fst-pop $indir/SA.noclones.inds \
    --weir-fst-pop $indir/SY.noclones.inds \
    --out PHAR_Fst_SA_SY

# -- SY vs SI -- #
vcftools --vcf $indir/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf \
    --weir-fst-pop $indir/SY.noclones.inds \
    --weir-fst-pop $indir/SI.noclones.inds \
    --out PHAR_Fst_SY_SI

# -- PAG vs GO -- # 
vcftools --vcf $indir/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf \
    --weir-fst-pop $indir/PAG.noclones.inds \
    --weir-fst-pop $indir/GO.noclones.inds \
    --out PHAR_Fst_PAG_GO
