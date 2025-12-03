#!/bin/bash

# Script for calculating nucleotide diversity and H-statistics to identify candidate regions of selection

# ----------------- NUCLEOTIDE DIVERSITY ----------------- #

# to identify regions where loss in nucleotide diversity (π) could be observed in one population but not the other (indicative of directional selection)
# we calculate π on the VCF file output by GATK with an additional flag to include all variant sites (--include-non-variant-sites) with pixy

indir="/path/to/vcf"
workdir="/path/to/sweeps/pi"
cd $workdir

# remove clone mates
vcftools --gzvcf $indir/PHAR.OUTPUT.GENOTYPE_ALL_SITES.vcf.gz --remove-indv UAE_SA_Phar_31 --recode --out temp1
vcftools --vcf temp1.recode.vcf --remove-indv UAE_SI_Phar_38 --recode --out PHAR.OUTPUT.GENOTYPE_ALL_SITES.noclones

# cleanup
rm temp1*
mv PHAR.OUTPUT.GENOTYPE_ALL_SITES.noclones.recode.vcf PHAR.OUTPUT.GENOTYPE_ALL_SITES.noclones.vcf

bgzip -k PHAR.OUTPUT.GENOTYPE_ALL_SITES.noclones.vcf
tabix PHAR.OUTPUT.GENOTYPE_ALL_SITES.noclones.vcf.gz

# 6,677,316 sites

# make a populations file for pixy
bcftools query -l PHAR.OUTPUT.GENOTYPE_ALL_SITES.noclones.vcf.gz  | awk 'BEGIN{OFS="\t"; print "SampleID","Population"} {split($0,a,"_"); print $0, a[2]}' > PI.popfile.noclones
tail -n +2 PI.popfile.noclones > PI.popfile.noclones.noheader

# run pixy on ALL_SITES.vcf file
conda activate pixy
outdir="/path/to/pi_output"

pixy --stats pi fst dxy \
    --vcf PHAR.OUTPUT.GENOTYPE_ALL_SITES.noclones.vcf.gz \
    --populations PI.popfile.noclones.noheader \
    --window_size 500 \
    --n_cores 32 \
    --output_folder $outdir 

# ------------------- ALLELE FREQUENCY ------------------- #

# use frequency files generated in /sweeps/01.selective_sweep_analyses.sh

awk 'NR==1 {print "scaff\tchrom\tpos\tnalleles\tnchr\tfreq1\tfreq2"; next} {print "Phar_scaff_"$1"\t"$0}' PAG.frq > PAG_mod.frq
awk 'NR==1 {print "scaff\tchrom\tpos\tnalleles\tnchr\tfreq1\tfreq2"; next} {print "Phar_scaff_"$1"\t"$0}' GO.frq > GO_mod.frq

vcftools --vcf PAG.SWEEP.noclones.recode.vcf --chr Phar_scaff_29 --recode --recode-INFO-all --out PAG.Phar_scaff_29
vcftools --vcf GO.SWEEP.noclones.recode.vcf --chr Phar_scaff_29 --recode --recode-INFO-all --out GO.Phar_scaff_29

# -------------------------- nSL ---------------------- #
# selscan v2.1

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

# ----------------------- iHH12 ------------------- #
# selscan v2.1

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
