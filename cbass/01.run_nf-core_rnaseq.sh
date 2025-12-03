#!/bin/bash

# Script to run nf-core/rnaseq pipeline on RNA-Seq raw files
# https://github.com/nf-core/rnaseq

# ------ FOLDER STRUCTURE ------ #

# UAE_PAG_WGS_PopGen/cbass/
#                           |
#                           --- Phar_ref/PAG_UKon_Phar.V2.fasta.gz
#                           --- Phar_ref/Porites_harrisoni_15_09.gtf.gz
#                           |
#                           --- 1_baseline/samplesheet.csv
#                           --- 1_baseline/nf-params.json
#                           --- 1_baseline/nextflow.config
#                                    |
#                                    --- results/
#                                    --- work/ (can be removed after -- will be very large)
#                           |
#                           --- 2_third_temp/samplesheet.csv
#                           --- 2_third_temp/nf-params.json
#                           --- 2_third_temp/nextflow.config
#                                    |
#                                    --- results/
#                                    --- work/ (can be removed after -- will be very large)

# ------ PREPARE REF FILES ------ #

# annotation gtf file
conda activate annotate

cd /home/fiesingera/proj/PAG_Phar_RNASeq_nf/Phar_ref/
scp /home/fiesingera/Phar_latest/*.fasta .

agat_convert_sp_gff2gtf.pl --gff /home/fiesingera/Phar_latest/Porites_harrisoni_15_09.gff3 -o Porites_harrisoni_15_09.gtf 
pigz *.fasta *.gtf 

# ------ 1 BASELINE TEMP PHAR ------ #

# PREPARE input files

echo "sample,fastq_1,fastq_2,strandedness" > samplesheet.csv

for pat in "20220530_UAE_SA_Phar_*_34.1_360" "20220530_UAE_SA_Phar_*_34.1_1080" "20220516_UAE_SY_Phar_*_33.9_360" "20220516_UAE_SY_Phar_*_33.9_1080" "20220523_UAE_SI_Phar_*_32.6_360" "20220523_UAE_SI_Phar_*_32.6_1080"; do
  for r1 in ${pat}_RNA-Seq_*_R1.fastq.gz; do
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    sample=$(echo "$r1" | sed -E 's@.*(UAE_[A-Z]{2}_Phar_[^_]+_[^_]+_(360|1080))_.*@\1@')
    if [[ -f "$r1" && -f "$r2" ]]; then
      echo "$sample,$PWD/$r1,$PWD/$r2,auto" >> samplesheet.csv
    fi
  done
done

# with populations column for DESeq2
echo "sample,fastq_1,fastq_2,strandedness,population" > samplesheet_popl.csv

for pat in "20220530_UAE_SA_Phar_*_34.1_360" "20220530_UAE_SA_Phar_*_34.1_1080" "20220516_UAE_SY_Phar_*_33.9_360" "20220516_UAE_SY_Phar_*_33.9_1080" "20220523_UAE_SI_Phar_*_32.6_360" "20220523_UAE_SI_Phar_*_32.6_1080"; do
  if [[ "$pat" == 20220523_UAE_SI_Phar_* ]]; then
    pop="GO"
  else
    pop="PAG"
  fi

  for r1 in ${pat}_RNA-Seq_*_R1.fastq.gz; do
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    sample=$(echo "$r1" | sed -E 's@.*(UAE_[A-Z]{2}_Phar_[^_]+_[^_]+_(360|1080))_.*@\1@')
    if [[ -f "$r1" && -f "$r2" ]]; then
      echo "$sample,$PWD/$r1,$PWD/$r2,auto,$pop" >> samplesheet_popl.csv
    fi
  done
done

# prepare nextflow config and nf-params.json files

# NEXTFLOW pipeline

conda activate nf_env
nextflow run nf-core/rnaseq -r 3.21.0 -name Phar_baseline -profile singularity -params-file nf-params.json 

# ------ 2 THIRD TEMP PHAR ------ #

# PREPARE input files

echo "sample,fastq_1,fastq_2,strandedness" > samplesheet.csv

for pat in "20220530_UAE_SA_Phar_*_40.1_360" "20220530_UAE_SA_Phar_*_40.1_1080" "20220516_UAE_SY_Phar_*_39.9_360" "20220516_UAE_SY_Phar_*_39.9_1080" "20220523_UAE_SI_Phar_*_38.6_360" "20220523_UAE_SI_Phar_*_38.6_1080"; do
  for r1 in ${pat}_RNA-Seq_*_R1.fastq.gz; do
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    sample=$(echo "$r1" | sed -E 's@.*(UAE_[A-Z]{2}_Phar_[^_]+_[^_]+_(360|1080))_.*@\1@')
    if [[ -f "$r1" && -f "$r2" ]]; then
      echo "$sample,$PWD/$r1,$PWD/$r2,auto" >> samplesheet.csv
    fi
  done
done

# with populations column for DESeq2
echo "sample,fastq_1,fastq_2,strandedness,population" > samplesheet_popl.csv

for pat in "20220530_UAE_SA_Phar_*_40.1_360" "20220530_UAE_SA_Phar_*_40.1_1080" "20220516_UAE_SY_Phar_*_39.9_360" "20220516_UAE_SY_Phar_*_39.9_1080" "20220523_UAE_SI_Phar_*_38.6_360" "20220523_UAE_SI_Phar_*_38.6_1080"; do
  if [[ "$pat" == 20220523_UAE_SI_Phar_* ]]; then
    pop="GO"
  else
    pop="PAG"
  fi

  for r1 in ${pat}_RNA-Seq_*_R1.fastq.gz; do
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    sample=$(echo "$r1" | sed -E 's@.*(UAE_[A-Z]{2}_Phar_[^_]+_[^_]+_(360|1080))_.*@\1@')
    if [[ -f "$r1" && -f "$r2" ]]; then
      echo "$sample,$PWD/$r1,$PWD/$r2,auto,$pop" >> samplesheet_popl.csv
    fi
  done
done

# prepare nextflow config and nf-params.json files

# NEXTFLOW pipeline

conda activate nf_env
nextflow run nf-core/rnaseq -r 3.21.0 -name Phar_third_temp -profile singularity -params-file nf-params.json 
