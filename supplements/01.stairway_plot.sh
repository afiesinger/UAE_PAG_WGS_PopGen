#!/bin/bash

# Script for generation of a stairway plot to infer historical effective population size

# ----------- STAIRWAY PLOT FOR PAG CORALS ---------- #

# Stairway Plot v2
# PAG samples = sites SA & SY 

workdir="/path/to/files"
assembly_dir="/path/to/reference_genome"

cd $workdir

FILTERS_PHAR_SA_SY="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 60 -snp_pval 1e-6 -hwe_pval 0.001"
TODO="-dosaf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -dobcf 1 --ignore-RG 0 -dogeno 1 -docounts 1"
BAMS_PHAR_SA_SY="ANGSD_bam_files_SA_SY.list"

angsd -b $BAMS_PHAR_SA_SY $FILTERS_PHAR_SA_SY $TODO -anc $assembly_dir/PAG_UKon_Phar_1.1.fasta -nThreads 8 -out PHAR_ANGSD_swp_SA_SY 

realSFS PHAR_ANGSD_swp_SA_SY.saf.idx -fold 1 > PHAR_ANGSD_swp_SA_SY_folded.sfs

# count sites for making the blueprint file
NSITES=`zcat PHAR_ANGSD_swp_SA_SY.mafs.gz | wc -l`
echo $NSITES
# 37,960,754

# make blueprint file
nano SA_SY.blueprint

# #input setting
# popid: SASY # id of the population (no white space)
# nseq: 158  # number of sequences
# L: 37960754 # total number of observed nucleic sites, including polymorphic and monomorphic
# whether_folded: true # whether the SFS is folded (true or false)
# SFS: 3571515.848398 4336978.792273 3105644.686619 2473684.738256 2013760.508678 1603327.888976 1360534.141045 1176253.218436 1018261.019333 902053.220204 804150.662814 727126.342810 663309.328481 603481.025213 561588.838497 516784.474952 483778.698401 444870.621199 424949.952682 394114.901818 377878.544736 359961.574485 345135.253941 323824.288639 309949.731688 293577.679133 285741.527264 276371.496247 265693.310677 253721.726003 242527.531140 236671.143913 228999.254693 222462.587956 218291.417988 213615.685255 204279.451018 196212.733395 192564.118715 189794.329491 187424.002459 181593.873994 174632.430041 172446.074076 169857.359378 167984.262973 160939.513113 155789.817815 158036.429849 158414.310885 155816.757260 150144.545115 147062.864875 143139.855706 141505.886457 140494.565925 140122.637686 141456.508288 140007.947578 137390.156538 134384.208745 132961.096415 132896.331273 132369.012443 129906.979603 129154.514510 129729.908224 129685.255024 130723.668258 128555.043005 127061.838642 129405.912051 130388.169359 130047.779267 128166.503382 127518.777717 126412.504635 125408.727646 78273.672033
# #smallest_size_of_SFS_bin_used_for_estimation: 2 # default is 1; to ignore singletons, change this number to 2
# #largest_size_of_SFS_bin_used_for_estimation: 44 # default is nseq/2 for folded SFS
# pct_training: 0.67 # percentage of sites for training
# nrand: 5    10  15  20  25 # number of random break points for each try (separated by white space)
# project_dir: SASY_g34_m483 # project directory
# stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
# ninput: 200 # number of input files to be created for each estimation
# #output setting
# mu: 4.83e-8 # assumed mutation rate per site per generation
# year_per_generation: 34 # assumed generation time (in years)
# #plot setting
# plot_title: SASY_fold_g34_m483_STAIRWAY # title of the plot
# xrange: 1,400 # Time (1k year) range; format: xmin,xmax; "0,0" for default
# yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
# xspacing: 2 # X axis spacing
# yspacing: 2 # Y axis spacing
# fontsize: 14 # Font size

# then run Stairbuilder
java -cp stairway_plot_es Stairbuilder SA_SY.blueprint 

# then run the .sh files to produce the results
bash SA_SY.blueprint.sh
