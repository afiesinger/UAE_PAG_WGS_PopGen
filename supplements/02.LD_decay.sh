#!/bin/bash

# ------------ LINKAGE DISEQUILIBRIUM DECAY ------------- #
# PopLDdecay v3.40

indir="/path/to/vcf/"

grep -vc '^#' PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.recode.vcf
# number of total SNPs: 1,140,169

/home/fiesingera/tools/PopLDdecay/bin/PopLDdecay -InVCF $indir/POPGEN/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.recode.vcf.gz \
    -OutStat PAG_PopLDdecay_out -SubPop $indir/PAG.all.inds

/home/fiesingera/tools/PopLDdecay/bin/PopLDdecay -InVCF $indir/POPGEN/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.recode.vcf.gz \
    -OutStat GO_PopLDdecay_out -SubPop $indir/GO.all.inds

perl /home/fiesingera/tools/PopLDdecay/bin/Plot_OnePop.pl -inFile PAG_PopLDdecay_out.stat.gz -output Fig_PAG_PopLDdecay
perl /home/fiesingera/tools/PopLDdecay/bin/Plot_OnePop.pl -inFile GO_PopLDdecay_out.stat.gz -output Fig_GO_PopLDdecay

nano Pop.ResultPath.list_PAG-GO
# /path/to/file/PAG_PopLDdecay_out.stat.gz  PAG
# /path/to/file/GO_PopLDdecay_out.stat.gz   GO

# the perl script Plot_MultiPop_edit.pl was edited so as to display custom colors for PAG and GO 
perl /home/fiesingera/tools/PopLDdecay/bin/Plot_MultiPop_edit.pl -inList Pop.ResultPath.list_PAG-GO -output Fig_PAG_GO_PopLDdecay
