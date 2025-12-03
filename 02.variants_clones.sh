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
