#!/bin/bash

# Script for variant calling and subsequent filtering from preprocessed WGS files & clone mate identification and removal

# ------------------ VARIANT CALLING ------------------

# index all bam files for variant calling
cd /path/to/mapped_files
samtools index -M *.sort.q30.RG.MD.bam

# mount docker image with a volume in target directory to access data
docker run -v $PWD:/gatk -it broadinstitute/gatk

# create .fasta dict file for reference genome & make sure there is a .fai file of the reference genome too
gatk CreateSequenceDictionary -R PAG_UKon_Phar_1.1.fasta

# call variants with HaplotypeCaller
# this takes very long for a high number of samples; consider running on a cluster or splitting files across multiple runs for parallelization
for F in *.sort.q30.RG.MD.bam ; do
    base=$(basename ${F} .sort.q30.RG.MD.bam)
    gatk HaplotypeCaller \
    --java-options "-Xmx4g -XX:ParallelGCThreads=4" \
    -I ${base}.sort.q30.RG.MD.bam \
    -O ${base}.raw.g.vcf \
    -R PAG_UKon_Phar_1.1.fasta \
    --emit-ref-confidence GVCF \
    --native-pair-hmm-threads 16
done

# -------------------- JOINT GENOTYPING ACROSS SAMPLES ------------------ #

# still inside the docker image

prefix="/path/to/files"
output="PHAR.cohort.g.vcf"

# First, create a list of all input gVCF files
input_gvcfs=""
for F in $prefix/*.raw.g.vcf; do
    input_gvcfs="$input_gvcfs --variant $F"
done

# Then, run CombineGVCFs once with all input files
gatk CombineGVCFs -R PAG_UKon_Phar_1.1.fasta $input_gvcfs -O $output

# Now, genotype the combined VCF file
gatk --java-options "-Xmx4g" GenotypeGVCFs -R PAG_UKon_Phar_1.1.fasta -V PHAR.cohort.g.vcf -O PHAR.OUTPUT.GENOTYPE.vcf

# check variant output in a readable format
gatk VariantsToTable -V PHAR.OUTPUT.GENOTYPE.vcf -F CHROM -F TYPE -F QUAL -F NO-CALL -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -O PHAR.variants.table

# ----------------- KEEP ONLY SNPs ------------------ #

# for our analyses we only consider SNPs
gatk SelectVariants -R PAG_UKon_Phar_1.1.fasta -V PHAR.OUTPUT.GENOTYPE.vcf --select-type-to-include SNP -O PHAR.OUT.SNPs.vcf

# how many SNPs remain?
grep -vc '^#' PHAR.OUT.SNPs.vcf

# ----------------- FILTER SNPs ------------------- #

docker run -v $PWD:/gatk/pag -it broadinstitute/gatk

# filter SNPs for DOWNSTREAM
gatk VariantFiltration \
    -R PAG_UKon_Phar_1.1.fasta \
    -V PHAR.OUT.SNPs.vcf \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "SOR > 4.0" --filter-name "StrandBias" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    --genotype-filter-expression "DP < 2 || DP > 238" --genotype-filter-name "DPfilt" \
    --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum" \
    -O PHAR.OUT.SNPs_tmp.vcf

gatk SelectVariants \
    -R PAG_UKon_Phar_1.1.fasta \
    -V PHAR.OUT.SNPs_tmp.vcf \
    --exclude-filtered \
    -O PHAR.OUT.SNPs_GATK.vcf

rm -f PHAR.OUT.SNPs_tmp.vcf*

# how many SNPs remain?
grep -vc '^#' PHAR.OUT.SNPs_GATK.vcf 

# ------- DOWNSTREAM APPLICATION SPECIFIC FILTERING OF VCF FILE ------- #
# vcftools v 0.1.16
# filter VCF file specifically for each downstream application using vcftools

# # # --- POPGEN --- # # #
cd /path/to/vcf_file

vcftools --vcf PHAR.OUT.SNPs_GATK.vcf --max-missing 0.8 --minQ 30 --min-alleles 2 --max-alleles 2 --hwe 0.001 --maf 0.05 --recode --recode-INFO-all --out PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi.vcf
mv PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi.recode.vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi.vcf

vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi.vcf --missing-indv

# list samples with more than 20 or 30% missingness
awk 'NR>1 && $5 > 0.2 {print $1}' out.imiss > remove.inds

# jump to 03.VCF_QC.R to plot missingness & find best threshold to use for filtering below

# remove them from the vcf file
vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi.vcf --remove remove.inds --recode --recode-INFO-all --out PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv

# LD pruning
vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.recode.vcf --geno-r2 --ld-window-bp 200 --min-r2 0.2 --out ld_out

# sites to remove
awk 'NR>1 {print $1"\t"$3}' ld_out.geno.ld | sort | uniq > snps_to_exclude.txt

vcftools --vcf PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.recode.vcf --exclude-positions snps_to_exclude.txt --recode --recode-INFO-all --out PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned

# how many SNPs remain?
grep -vc '^#' PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.recode.vcf
# 663,850

# go to 04.kinship.sh and then 05.kinship_plots.R to identify clone mates and return here to filter VCF file for selective sweep analyses

# # # --- SELECTIVE SWEEPS --- # # #

vcftools --vcf PHAR.OUT.SNPs_GATK.vcf --max-missing 0.8 --minQ 30 --min-alleles 2 --max-alleles 2 --maf 0.002 --recode --recode-INFO-all --out PHAR.OUT.SNPs_filt_sweeps
mv PHAR.OUT.SNPs_filt_sweeps.recode.vcf PHAR.OUT.SNPs_filt_sweeps.vcf

vcftools --vcf PHAR.OUT.SNPs_filt_sweeps.vcf --missing-indv
sort -k5 -n out.imiss > out.sort # look at highest missingness

# list samples with more than 20% missingness
awk 'NR>1 && $5 > 0.2 {print $1}' out.imiss > remove.inds

# add clones to the remove.inds file #

nano remove.inds
# UAE_SA_Phar_31
# UAE_SI_Phar_38

# remove them from the vcf file
vcftools --vcf PHAR.OUT.SNPs_filt_sweeps.vcf --remove remove.inds --recode --recode-INFO-all --out PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES
mv PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.recode.vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.vcf

# how many SNPs remain?
grep -vc '^#' PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.vcf 
# 5,185,934

# LD pruning
vcftools --vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.vcf --geno-r2 --ld-window-bp 200 --min-r2 0.2 --out ld_out

# sites to remove
awk 'NR>1 {print $1"\t"$3}' ld_out.geno.ld | sort | uniq > snps_to_exclude.txt

vcftools --vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.vcf  --exclude-positions snps_to_exclude.txt --recode --recode-INFO-all --out PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned

# how many SNPs remain?
grep -vc '^#' PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf 
# 2,878,739

bgzip -k PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf 
tabix PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf.gz

# make popfiles for SWEEP ANALYSES 

bcftools query -l PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf.gz  | awk 'BEGIN{OFS="\t"; print "SampleID","Population"} {split($0,a,"_"); print $0, a[2]}' > SWEEP.popfile.noclones
tail -n +2 SWEEP.popfile.noclones > SWEEP.popfile.noclones.noheader

# PAG
awk '$2=="SA" || $2=="SY" {print $1}' SWEEP.popfile.noclones.noheader > PAG.noclones.inds

# GO
awk '$2=="SI" {print $1}' SWEEP.popfile.noclones.noheader > GO.noclones.inds

# SA
awk '$2=="SA" {print $1}' SWEEP.popfile.noclones.noheader > SA.noclones.inds

# SY
awk '$2=="SY" {print $1}' SWEEP.popfile.noclones.noheader > SY.noclones.inds

# SI
awk '$2=="SI" {print $1}' SWEEP.popfile.noclones.noheader > SI.noclones.inds

# make separate VCF files for PAG & GO
vcftools --vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf --keep PAG.noclones.inds --recode --recode-INFO-all --out PAG.SWEEP.noclones
vcftools --vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf --keep GO.noclones.inds --recode --recode-INFO-all --out GO.SWEEP.noclones

