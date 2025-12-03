#!/bin/bash

# -------------------------- Fst VCFTOOLS ------------------------ #
vcftools v0.1.16

indir="/path/to/sweeps"

# be careful here; sweeps file has different number of samples than popgen file due to missingness removal !

# -- GLOBAL FST -- #

vcftools --vcf $indir/PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf \
    --weir-fst-pop $indir/PAG.noclones.inds \
    --weir-fst-pop $indir/GO.noclones.inds \
    --out PHAR_Fst_ALL

# -- WINDOWED FST (500 bp windows) -- #

vcftools --vcf $indir/PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf \
    --weir-fst-pop $indir/PAG.noclones.inds \
    --weir-fst-pop $indir/GO.noclones.inds \
    --fst-window-size 500 \
    --fst-window-step 250 \
    --out PHAR_Fst_500bp

# ------------------------------------- H-STATISTICS ------------------------------------------ #

###################################
# BEAGLE PHASE & IMPUTE GENOTYPES #
###################################

# first we have to phase the genotypes

indir="/path/to/sweeps"
workdir="/path/to/sweeps/EHH"
cd $workdir

# first, fix haploid genotypes
bcftools +fixploidy $indir/PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf -o PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.tmp.vcf

bcftools query -f '%CHROM\n' PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.tmp.vcf | sort | uniq -c | awk '$1 == 1 {print $2}'

# iteratively remove these contigs:

chroms=(
Phar_scaff_1048
Phar_scaff_1169
Phar_scaff_1287
Phar_scaff_1595
Phar_scaff_1599
Phar_scaff_1603
Phar_scaff_1607
Phar_scaff_1611
Phar_scaff_1622
Phar_scaff_1644
Phar_scaff_1664
Phar_scaff_1677
Phar_scaff_1692
Phar_scaff_1703
Phar_scaff_1729
Phar_scaff_1806
Phar_scaff_1819
Phar_scaff_1821
Phar_scaff_1835
Phar_scaff_1836
Phar_scaff_1879
)

invcf="PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.tmp.vcf"
outvcf="PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean"

tmpvcf="$invcf"

for chr in "${chroms[@]}"; do
    vcftools --vcf "$tmpvcf" --not-chr "$chr" --recode --recode-INFO-all --out "${outvcf}_$chr"
    tmpvcf="${outvcf}_$chr.recode.vcf"
done

mv PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean_Phar_scaff_1879.recode.vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.vcf
rm PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean_Phar_scaff*
rm PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.tmp.vcf

# phase genotypes with genotype imputation 
# BEAGLE v5.5
java -Xmx16g -jar /home/fiesingera/tools/beagle.27Feb25.75f.jar \
  gt=PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.vcf \
  out=PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd \
  impute=true \
  nthreads=32

gunzip PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.vcf.gz

# TAKE THE IMPUTED FILE 

indir="/path/to/sweeps/EHH"
indir_inds="/path/to/sweeps"
outdir_ihs="/path/to/sweeps/EHH/iHS_OUT"
outdir_xpehh="/path/to/sweeps/XPEHH_OUT"
outdir_nsl="/path/to/sweeps/EHH/nSL_OUT"
outdir_ihh12="/path/to/sweeps/EHH/iHH12_OUT"

cd $indir

# reformat file to add SNP ID column from scaffold name and position -- then make scaffold names numeric
cat PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.vcf | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$3 = $1 "_" $2; print}' > TMP.vcf
cat TMP.vcf | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {sub(/^Phar_scaff_/, "", $1); print}' > PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf

rm TMP.vcf

# make map file for selscan & scaffold file
cat PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf | grep -v "#" | awk '{print $1"\t"$3"\t"$2"\t"$2}' > PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.map
grep -v "^#" PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf | cut -f1 | sort -u > PHAR.OUT_scaffold.list

###################################
# IDENTIFY DOMINANT ALLELE IN PAG #
###################################

# SPLIT VCF FILE BY PAG vs GO  & CALCULATE ALLELE FREQ 

vcftools --vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf --keep $indir_inds/PAG.noclones.inds --freq --out PAG
vcftools --vcf PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf --keep $indir_inds/GO.noclones.inds --freq --out GO

# identify the dominant allele in the OUT population (i.e. PAG). To do this, take the frq output, change : character to tab as allele+freq are in the form a:f. 
# If statement: if frequency of 1st allele is >= freq of 2nd allele, print out 'chr pos 1st_allele' else print 'chr pos 2nd_allele'. Save to file for dominant OUT allele.
cat PAG.frq | sed -E 's/:/\t/g' | while read line ; do awk '{
if ($6 >= $8)
print $1,$2,$5
else
print $1,$2,$7
}'; done > PAG.dom.allele

# print snpid next to OUT dom allele info
cat PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf | grep -v "#" |  awk '{print $3}' | paste - PAG.dom.allele > PAG.dom.allele.newref

# create headerless version of the vcf
grep -v "#" PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf > PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT_nohead.vcf

# change space to tab in OUT dom file, paste with headerless vcf. If dom allele is ref allele in vcf, print out vcf line as is
sed -E 's/\ /\t/g' PAG.dom.allele.newref | paste - PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT_nohead.vcf | while IFS=$'\t' read idd chrd posd dom chr pos id ref alt other; do if [ $dom == $ref ]; then echo -e $chr"\t"$pos"\t"$id"\t"$ref"\t"$alt"\t"$other; fi; done > SNP.dom.match.vcf
 
# change space to tab in OUT dom file, paste with headerless vcf. If dom allele is not ref allele in vcf, print out vcf line with ref and alt alleles switched in info
sed -E 's/\ /\t/g' PAG.dom.allele.newref | paste - PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT_nohead.vcf | while IFS=$'\t' read idd chrd posd dom chr pos id ref alt other; do if [ $dom != $ref ]; then echo -e $chr"\t"$pos"\t"$id"\t"$alt"\t"$ref"\t"$other; fi; done > SNP.dom.mismatch.vcf

# recode genotypes to reflect changed reference allele. Switch genotypes using / instead | to avoid genotypes that have alread been changed. I.e. 0|0 becomes 1/1, 0|1 becomes 1/0, 1|0 becomes 0/1, 1|1 becomes 0/0 - then change / back to | to keep phase. 
cat SNP.dom.mismatch.vcf | sed -E 's/0\|0/1\/1/g' | sed -E 's/0\|1/1\/0/g' | sed -E 's/1\|0/0\/1/g' | sed -E 's/1\|1/0\/0/g' | sed 's/\//\|/g' > SNP.dom.mismatch.final.vcf

# sort vcf and add header
grep "#" PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.vcf > vcf.phased.header
cat SNP.dom.match.vcf SNP.dom.mismatch.final.vcf | sort -n -k1,1 -k2,2 | cat vcf.phased.header - > SNP.dom.phased.vcf
sed 's/ \+/\t/g' SNP.dom.phased.vcf > SNP.dom.phased.corrected.vcf # fix spaces to tabs

# cleanup
rm PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT_nohead.vcf
rm vcf.phased.header

# # REMAKE VCF FILE FOR iHS & XP-EHH CALCULATIONS # #
vcftools --vcf SNP.dom.phased.corrected.vcf --keep $indir_inds/PAG.noclones.inds --recode --recode-INFO-all --out PAG.dom
vcftools --vcf SNP.dom.phased.corrected.vcf --keep $indir_inds/GO.noclones.inds --recode --recode-INFO-all --out GO.dom

################
# iHS & XP-EHH #
################

while read scaffold; do
    # extract scaffold-specific VCFs 
    grep -w "^$scaffold" PAG.dom.recode.vcf > tmp.PAG.vcf
    grep -w "^$scaffold" GO.dom.recode.vcf > tmp.GO.vcf

    # extract corresponding entries from the map file
    grep -w "^$scaffold" PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.map > tmp.map

    # run iHS on PAG
    /path/to/selscan --ihs \
        --vcf tmp.PAG.vcf \
        --map tmp.map \
        --out $outdir_ihs/${scaffold}_PAG \
        --trunc-ok

    # run iHS on GO
    /path/to/selscan --ihs \
        --vcf tmp.GO.vcf \
        --map tmp.map \
        --out $outdir_ihs/${scaffold}_GO \
        --trunc-ok
    
    # run XP-EHH on PAG vs GO
    /path/to/selscan --xpehh \
        --vcf tmp.PAG.vcf \
        --vcf-ref tmp.GO.vcf \
        --map tmp.map \
        --out $outdir_xpehh/${scaffold}"_XPEHH_PAGvGO" \
        --trunc-ok
done < PHAR.OUT_scaffold.list

# normalize data to account for genome-wide differences in haplotype length between PAG & GO #
rm tmp*

cd iHS_OUT/

# # PAG # #
files_pag=$(cat PHAR.OUT_scaffold.list | sed -E 's/$/_PAG.ihs.out/' | while read f; do [ -f "$f" ] && echo "$f"; done | tr '\n' ' ')
/path/to/selscan/src/norm --ihs --files $files_pag

# concatenate all normalized files for outlier calculation
cat *_PAG.ihs.out.100bins.norm > ALL_PAG.ihs.out.norm

# add a header
(echo -e "id\tpos\tdaf\tihh1\tihh0\tuihs\tsihs\tbin"; cat ALL_PAG.ihs.out.norm) > ALL_PAG.ihs.out.norm.header

# # GO # #
files_go=$(cat PHAR.OUT_scaffold.list | sed -E 's/$/_GO.ihs.out/' | while read f; do [ -f "$f" ] && echo "$f"; done | tr '\n' ' ')
/path/to/selscan/src/norm --ihs --files $files_go

# concatenate all normalized files for outlier calculation
cat *_GO.ihs.out.100bins.norm > ALL_GO.ihs.out.norm

# add a header
(echo -e "id\tpos\tdaf\tihh1\tihh0\tuihs\tsihs\tbin"; cat ALL_GO.ihs.out.norm) > ALL_GO.ihs.out.norm.header

# move all individual files to ./ind_files/
mkdir ind_files
mv *.log *.out *.100bins.norm ind_files/

# # PAG vs GO # #
cd XPEHH_OUT/

# prepare normalized files
files_xpehh=$(cat PHAR.OUT_scaffold.list | sed -E 's/$/_XPEHH_PAGvGO.xpehh.out/' | while read f; do [ -f "$f" ] && echo "$f"; done | tr '\n' ' ')
/path/to/selscan/src/norm --xpehh --files $files_xpehh 

# Ccncatenate normalized files

# pick one file to extract the header
header_file="1_XPEHH_PAGvGO.xpehh.out.norm"
output_file="XPEHH_ALL_SCAFF.xpehh.out.norm"

# save header
head -n 1 $header_file > header.tmp

# remove header lines from all files in-place
for f in *.xpehh.out.norm; do
    sed -i '1d' $f
done

# concatenate all headerless files
cat *.xpehh.out.norm > data.tmp

# add the single header back
cat header.tmp data.tmp > $output_file

# clean up
rm -f header.tmp data.tmp

# move all individual files to ind_files/
mkdir ind_files
mv *.log *.out *_PAGvGO.xpehh.out.norm ind_files/

# from Rscript the same table but with cumulative positions along the genome -- for comparison with other metrics ! 
# XPEHH_ALL_SCAFF.xpehh.out.norm_CUM_POSITIONS

# # OUTLIERS # #

# get total number of data lines (excluding header)
total=$(grep -v "pos" XPEHH_ALL_SCAFF.xpehh.out.norm_CUM_POSITIONS | wc -l) # 2878718

# compute the number of top 1% entries
top_n=$(echo "$total * 0.01" | bc | cut -d'.' -f1) # 28787

# absolute value of normxpehh 
grep -v "pos" XPEHH_ALL_SCAFF.xpehh.out.norm_CUM_POSITIONS \
    | awk '{print sqrt($9^2)}' \
    | sort -g | tail -n "$top_n" | head -n1
# 3.09011

# get top 1% of rows directly, sorted by absolute normxpehh
grep -v "pos" XPEHH_ALL_SCAFF.xpehh.out.norm_CUM_POSITIONS \
    | awk '{print sqrt($9^2) "\t" $0}' \
    | sort -k1,1g | tail -n "$top_n" \
    | cut -f2- \
    | awk '$9 > 0' > top.1pc.normxpehh.PAG.outliers

# add header
(echo -e "id\tpos\tgpos\tp1\tihh1\tp2\tihh2\txpehh\tnormxpehh\tcrit\tscaff\tpos2\tscaffnum\tBPcumMb"; cat top.1pc.normxpehh.PAG.outliers) > top.1pc.normxpehh.PAG.outliers.header

# add chromosome info
awk 'BEGIN {OFS="\t"} NR==1 {print "chr", $0} NR>1 {split($1, a, "_"); print a[3], $0}' top.1pc.normxpehh.PAG.outliers.header > top.1pc.normxpehh.PAG.outliers.header.chr

## use this one: top.1pc.xpehh.PAG.outliers.header.chr

# GET SNP_IDs ("Phar_scaff_X_BPcum (in bp)")

awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="BPcumBP"} NR>1{$(NF+1)=$NF*1000000} 1'  top.1pc.normxpehh.PAG.outliers.header.chr > tmp.tsv
awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="snpid"} NR>1{$(NF+1)=$12"_"$NF} 1' tmp.tsv > top.1pc.normxpehh.PAG.outliers.header.chr.SNPids

rm tmp.tsv

###############
# nSL & iHH12 #
###############

# LOOKING FOR SOFT SWEEP RATHER THAN A HARD SWEEP WITH HAPLOTYPE STATISTICS:

# # nSL # #
while read scaffold; do
    # extract scaffold-specific VCFs 
    grep -w "^$scaffold" PAG.dom.recode.vcf > tmp.PAG.vcf
    grep -w "^$scaffold" GO.dom.recode.vcf > tmp.GO.vcf

    # extract corresponding entries from the map file
    grep -w "^$scaffold" PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.map > tmp.map

    # Run nSL on PAG
    /path/to/selscan --nsl \
        --vcf tmp.PAG.vcf \
        --map tmp.map \
        --out $outdir_nsl/${scaffold}_PAG \
        --trunc-ok

    # Run nSL on GO
    /path/to/selscan --nsl \
        --vcf tmp.GO.vcf \
        --map tmp.map \
        --out $outdir_nsl/${scaffold}_GO \
        --trunc-ok
done < PHAR.OUT_scaffold.list

rm tmp*

# normalize data to account for genome-wide differences in haplotype length between PAG & GO #
cd nSL_OUT/

nsl_pag=$(cat PHAR.OUT_scaffold.list | sed -E 's/$/_PAG.nsl.out/' | while read f; do [ -f "$f" ] && echo "$f"; done | tr '\n' ' ')
/home/fiesingera/tools/selscan/src/norm --nsl --files $nsl_pag

# concatenate all normalized files for outlier calculation
cat *_PAG.nsl.out.100bins.norm > ALL_PAG.nsl.out.norm

# add a header
(echo -e "id\tpos\tdaf\tsl1\tsl0\tunsl\tnormnsl\tbin"; cat ALL_PAG.nsl.out.norm) > ALL_PAG.nsl.out.norm.header

nsl_go=$(cat PHAR.OUT_scaffold.list | sed -E 's/$/_GO.nsl.out/' | while read f; do [ -f "$f" ] && echo "$f"; done | tr '\n' ' ')
/path/to/selscan/src/norm --nsl --files $nsl_go

# concatenate all normalized files for outlier calculation
cat *_GO.nsl.out.100bins.norm > ALL_GO.nsl.out.norm

# add a header
(echo -e "id\tpos\tdaf\tsl1\tsl0\tunsl\tnormnsl\tbin"; cat ALL_GO.nsl.out.norm) > ALL_GO.nsl.out.norm.header

mkdir ind_files
mv *.log *.out *.100bins.norm ind_files/

# # iHH12 # #
while read scaffold; do
    # extract scaffold-specific VCFs 
    grep -w "^$scaffold" PAG.dom.recode.vcf > tmp.PAG.vcf
    grep -w "^$scaffold" GO.dom.recode.vcf > tmp.GO.vcf

    # extract corresponding entries from the map file
    grep -w "^$scaffold" PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.clean.phased.imptd.REFORMAT.map > tmp.map
    
    # Run iHH12 on PAG
    /path/to/selscan --ihh12 \
        --vcf tmp.PAG.vcf \
        --map tmp.map \
        --out $outdir_ihh12/${scaffold}_PAG \
        --trunc-ok

    # Run iHH12 on GO
    /path/to/selscan --ihh12 \
        --vcf tmp.GO.vcf \
        --map tmp.map \
        --out $outdir_ihh12/${scaffold}_GO \
        --trunc-ok
done < PHAR.OUT_scaffold.list

rm tmp*

cd iHH12_OUT/

ihh12_pag=$(cat PHAR.OUT_scaffold.list | sed -E 's/$/_PAG.ihh12.out/' | while read f; do [ -f "$f" ] && echo "$f"; done | tr '\n' ' ')
/home/fiesingera/tools/selscan/src/norm --ihh12 --files $ihh12_pag

# concatenate normalized files

# pick one file to extract the header
header_file="1_PAG.ihh12.out.norm"
output_file="ALL_PAG.ihh12.out.norm"

# save header
head -n 1 $header_file > header.tmp

# remove header lines from all files in-place
for f in *_PAG.ihh12.out.norm; do
    sed -i '1d' $f
done

# concatenate all headerless files
cat *_PAG.ihh12.out.norm > data.tmp

# add the single header back
cat header.tmp data.tmp > $output_file

# clean up
rm -f header.tmp data.tmp

ihh12_go=$(cat PHAR.OUT_scaffold.list | sed -E 's/$/_GO.ihh12.out/' | while read f; do [ -f "$f" ] && echo "$f"; done | tr '\n' ' ')
/path/to/selscan/src/norm --ihh12 --files $ihh12_go

# concatenate normalized files

# pick one file to extract the header
header_file="1_GO.ihh12.out.norm"
output_file="ALL_GO.ihh12.out.norm"

# save header
head -n 1 $header_file > header.tmp

# remove header lines from all files in-place
for f in *_GO.ihh12.out.norm; do
    sed -i '1d' $f
done

# concatenate all headerless files
cat *_GO.ihh12.out.norm > data.tmp

# add the single header back
cat header.tmp data.tmp > $output_file

# clean up
rm -f header.tmp data.tmp

mkdir ind_files
mv *.log *.out* ind_files/
mv ind_files/ALL* .

# ----------------------- SNPEFF ------------------------- #

genome_dir="/path/to/reference_genome"
cd $genome_dir

# edit scaffold names to reflect the original assembly submission: Porites_harrisoni_2024.gff3

agat_convert_sp_gff2gtf.pl --gff Porites_harrisoni_2024.gff3 -gtf_version 2.2 -o Porites_harrisoni_2024.gtf

agat_sp_extract_sequences.pl -gff Porites_harrisoni_2024.gff3 -f PAG_UKon_Phar_1.1.fasta -t cds -o Porites_harrisoni_2024.cds.fa
agat_sp_extract_sequences.pl --gff Porites_harrisoni_2024.gff3 -f PAG_UKon_Phar_1.1.fasta -p -o Porites_harrisoni_2024.protein.fa

# remove asterisks from prot file
sed 's/\*//g' Porites_harrisoni_2024.protein.fa > Porites_harrisoni_2024.protein.fa.tmp
mv Porites_harrisoni_2024.protein.fa.tmp Porites_harrisoni_2024.protein.fa

workdir="/path/to/workdir"
cd $workdir

scp $genome_dir/Porites_harrisoni_2024.gtf .
mv Porites_harrisoni_2024.gtf genes.gtf

scp PAG_UKon_Phar_1.1.fasta .
mv PAG_UKon_Phar_1.1.fasta sequences.fa

scp $genome_dir/Porites_harrisoni_2024.protein.fa .
mv Porites_harrisoni_2024.protein.fa protein.fa

scp $genome_dir/Porites_harrisoni_2024.cds.fa .
mv Porites_harrisoni_2024.cds.fa cds.fa

# add the Phar genome to the config file

cd /path/to/snpEff/
vi snpEff.config

# # Porites harrisoni (Phar)
# Phar.genome : Porites_harrisoni

# build database

conda activate wgs
cd /path/to/snpEff/

java -jar snpEff.jar build -gtf22 -v Phar

# annotate SNPs on scaff 29 and 62
INFILE="PHAR.OUT.SNPs_filt_sweeps_indv.NOCLONES.ld_pruned.recode.vcf"
outdir="/path/to/snpEff"
cd $outdir

vcftools --vcf $INFILE --chr Phar_scaff_29 --recode --recode-INFO-all --out SNPEff.INPUT.Phar_scaff_29
vcftools --vcf $INFILE --chr Phar_scaff_62 --recode --recode-INFO-all --out SNPEff.INPUT.Phar_scaff_62

java -Xmx8g -jar snpEff.jar ann Phar $outdir/SNPEff.INPUT.Phar_scaff_29.recode.vcf -o vcf > $outdir/SNPEff.OUTPUT.Phar_scaff_29.vcf
mv snpEff_genes.txt SNPEff.OUTPUT.Phar_scaff_29_genes.txt
mv snpEff_summary.html SNPEff.OUTPUT.Phar_scaff_29_summary.html

java -Xmx8g -jar snpEff.jar ann Phar $outdir/SNPEff.INPUT.Phar_scaff_62.recode.vcf -o vcf > $outdir/SNPEff.OUTPUT.Phar_scaff_62.vcf
mv snpEff_genes.txt SNPEff.OUTPUT.Phar_scaff_62_genes.txt
mv snpEff_summary.html SNPEff.OUTPUT.Phar_scaff_62_summary.html

mv *.txt *.html $outdir
cd $outdir

# --- SNPSift to filter variants --- #

workdir="/path/to/snpEff"

java -jar SnpSift.jar extractFields $workdir/SNPEff.OUTPUT.Phar_scaff_29.vcf CHROM POS REF ALT ANN > $workdir/SNPEff.OUTPUT.Phar_scaff_29.ANN-OUT.txt
java -jar SnpSift.jar extractFields $workdir/SNPEff.OUTPUT.Phar_scaff_62.vcf CHROM POS REF ALT ANN > $workdir/SNPEff.OUTPUT.Phar_scaff_62.ANN-OUT.txt


cd $workdir
