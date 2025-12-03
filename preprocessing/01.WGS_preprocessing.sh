#!/bin/bash

# Script for preprocessing of WGS reads

# ----------------------------------- TRIMMING & ADAPTER REMOVAL  -------------------------------
# using bbduk to trim & remove adapter sequences
# bbmap v39.08

prefix="/path/to/files"
adapters_prefix="/path/to/bbmap"
output="/path/to/output"

for R1 in $prefix/*_R1.fastq.gz; do
    R2="${R1/_R1/_R2}"
    base=$(basename "$R1" _R1.fastq.gz)

    bbduk.sh \
    in1=$prefix/${base}_R1.fastq.gz \
    in2=$prefix/${base}_R2.fastq.gz \
    out1=$output/${base}_R1_trimmed.fastq.gz \
    out2=$output/${base}_R2_trimmed.fastq.gz \
    ref=$adapters_prefix/adapters.fa \
    minlen=100 \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo ;
done

# -------------------------- MAP TO REF GENOME & BASE QUALITY SCORE >= 30 ----------------------------

# samtools v1.21 & bwa v0.7.18

# Before you can align your reads to a reference genome, you need to create an index
# This only needs to be completed once per reference genome. BWA indexes are made from a FASTA genome file using bwa index:
refgenome_dir="/path/to/reference_genome"

# indexing Phar ref genome
bwa index $refgenome_dir/PAG_UKon_Phar_1.1.fasta.gz

prefix="/path/to/trimmed_files"
output="/path/to/mapped_files"

for R1 in $prefix/*_R1_trimmed.fastq.gz; do
    R2=${R1/_R1/_R2}
    base=$(basename $R1 _R1_trimmed.fastq.gz)

    # Align reads
    bwa mem -t 32 $refgenome_dir/PAG_UKon_Phar_1.1.fasta.gz $R1 $R2 > $output/${base}.sam
    
    # Convert SAM to BAM
    samtools view -bS $output/${base}.sam > $output/${base}.bam
    
    # Sort BAM file
    samtools sort -@ 32 $output/${base}.bam -o $output/${base}.sorted.bam
    
    # Optional: Quality filtering (consider removing or adjusting the -q parameter)
    samtools view -q 30 -b $output/${base}.sorted.bam > $output/${base}.sort.q30.bam
    
    # Clean up intermediate files
    rm "$output/${base}.sam" "$output/${base}.bam"

done

# ----------------- COVERAGE OF Q30 MAPPED READS --------------------------

# get the mean coverage of each .sort.q30.bam file and determine whether to use Q30 as a cutoff or to adjust the filtering according to quality

# without q30 filtering
output_file="mean_coverage_results_all.txt"
echo "File,Mean Coverage" > $output_file

for file in *.sorted.bam; do
    mean_coverage=$(samtools coverage $file | awk 'NR>1 {sum+=$6; count++} END {print sum/count}')
    echo "$file,$mean_coverage" >> $output_file
done

# with q30 filtering
output_file="mean_coverage_results_q30.txt"
echo "File,Mean Coverage" > $output_file

for file in *.sort.q30.bam; do
    mean_coverage=$(samtools coverage $file | awk 'NR>1 {sum+=$6; count++} END {print sum/count}')
    echo "$file,$mean_coverage" >> $output_file
done

# MEAN READ COVERAGE FOR FOR ALL SAMPLES #
awk -F, 'NR>1 {sum+=$2; count++} END {if(count>0) print sum/count; else print "No data"}' mean_coverage_results_q30.txt

# --------------------------- ADD READ GROUPS ------------------------------

# PICARD v3.1.1

PICARD="/path/to/picard.jar"
prefix="/path/to/mapped_files"

for F in $prefix/*.sort.q30.bam ; do
    base=$(basename ${F} .sort.q30.bam)
    
    # Extract relevant parts from the filename
    date=$(echo $base | cut -d'_' -f1)
    location=$(echo $base | cut -d'_' -f2,3)
    species=$(echo $base | cut -d'_' -f4)
    colony=$(echo $base | cut -d'_' -f5)
    lib_codes=$(echo $base | cut -d'_' -f10,11)
    platform=$(echo $base | cut -d'_' -f12)
    
    # Construct the read group fields
    RGID="${date}_${location}_${species}_${colony}"
    RGLB="${lib_codes}"
    RGSM="${location}_${species}_${colony}"
    RGPU="${platform}"

    java -jar $PICARD AddOrReplaceReadGroups \
    I=${F} \
    O=${prefix}/${base}.sort.q30.RG.bam \
    RGID=${RGID} \
    RGLB=${RGLB} \
    RGPL=ILLUMINA \
    RGPU=${RGPU} \
    RGSM=${RGSM}
done

# --------------------------- MARK DUPLICATES ------------------------------

# sambamba v1.0.1 

prefix="/path/to/mapped_files"

for F in $prefix/*.sort.q30.RG.bam ; do
    base=$(basename ${F} .sort.q30.RG.bam)

    sambamba markdup -t 32 -p $F ${prefix}/${base}.sort.q30.RG.MD.bam

done
