#!/bin/bash

# Script to run ADMIXTURE on VCF file without clone mates

# --------------- ADMIXTURE ON POPGEN NO CLONES FILE --------------- #
PLINK v1.90

# prepare .bed file from VCF
plink --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --make-bed --out Phar_ADMIX

# recode chromosome names so that they match input required for ADMIXTURE
awk '{$1="0";print $0}' Phar_ADMIX.bim > Phar_ADMIX.bim.tmp
mv Phar_ADMIX.bim.tmp Phar_ADMIX.bim

# run ADMIXTURE for K=1 through K=10 with 10 replicates each with random seeds

for K in {1..10}
do
  for run in {1..10}
  do
    SEED=$RANDOM
    admixture -s $SEED --cv Phar_ADMIX.bed $K | tee log${K}_run${run}.out
    echo "Seed used: $SEED" >> log${K}_run${run}.seed.out
    mv Phar_ADMIX.${K}.Q Phar_ADMIX.${K}_run${run}.Q
    mv Phar_ADMIX.${K}.P Phar_ADMIX.${K}_run${run}.P
  done
done

mkdir seeds
mv *.seed.out seeds/

# find optimal K based on cross-validation error
grep -h CV log*.out | sed -E 's/.*K=([0-9]+).*: ([0-9.eE+-]+)/\1 \2/' > cv_error_table.txt

# -------------- CLUMPP TO PERMUTE PROBABILITY OF ASSIGNMENT ------------- #

# make CLUMPP input files from ADMIXTURE output using the R package pophelper in 04.ADMIXTURE_plot.R

# then run CLUMPP like below

CLUMPP="/path/to/CLUMPP_executable"

# change the K here in each instance, then run CLUMPP in each directory
indir="/path/to/CLUMPP/pop_K1"
cd $indir

$CLUMPP paramfile

cd ..
find . -type f -name "*-combined-merged.txt" | while read line; do mv $line $PWD; done

# then back to 04.ADMIXTURE_plot.R for making the plot(s)
