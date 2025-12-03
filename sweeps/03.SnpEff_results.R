#!/usr/bin/env Rscript

# Script to identify outlier SNPs in annotated SNPs from SNPEff & SNPSift

# ----------- SNPEff & SNPSift --------- #

library(dplyr)
library(tidyr)
library(readr)

setwd("")

# --- CONTIG 29 --- #

ann.29 <- read.table("SNPEff.OUTPUT.Phar_scaff_29.ANN-OUT.txt", header=TRUE, sep="\t")

ann.29$CHR <- as.character(ann.29$CHROM)
ann.29$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", ann.29$CHR))

# load scaffold lengths and cumulative offsets 
scaffold_lengths <- read.table("PAG_UKon_Phar_Genome_scaffold_lengths.txt.header", header=TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

# calculate cumulative position for SNP (in bp and in Mb)
ann.29$BPcumBP <- ann.29$POS + as.numeric(offsets[as.character(ann.29$CHR)])
ann.29$BPcum_Mb <- ann.29$BPcumBP / 1e6

# create unique SNP ID
ann.29$snpid <- paste0(ann.29$CHR, "_", ann.29$BPcumBP)

xpehh.out <- read.table("top.1pc.normxpehh.PAG.outliers.header.chr.SNPids", header = TRUE, sep="", stringsAsFactors=FALSE)
fst.out <- read.table("PHAR_Fst_ALL.OUTLIERS_0.1p_with_SNPIDs.txt", header=TRUE, sep="", stringsAsFactors=FALSE)

# filter for contig 29
xpehh.out.29 <- xpehh.out %>% filter(chr == 29)
fst.out.29 <- fst.out %>% filter(CHROM == "Phar_scaff_29")

ovl_xa_29 <- intersect(as.character(unique(xpehh.out.29$snpid)), as.character(unique(ann.29$snpid)))
ovl_xfa_29 <- intersect(as.character(unique(fst.out.29$snpid)), ovl_xa_29)

print(ovl_xa_29) # 696
print(ovl_xfa_29) # 50

ann_29_OUT <- ann.29[ann.29$snpid %in% ovl_xfa_29, ]
write.table(ann_29_OUT, "SNPEff.OUTPUT.Phar_scaff_29.ANN-OUT_OUTLIER-SNPs.txt", , sep="\t", quote=FALSE, row.names=FALSE)

# --- CONTIG 62 --- #

ann.62 <- read.table("SNPEff.OUTPUT.Phar_scaff_62.ANN-OUT.txt", header=TRUE, sep="\t")

ann.62$CHR <- as.character(ann.62$CHROM)
ann.62$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", ann.62$CHR))

# load scaffold lengths and cumulative offsets (as in your earlier code)
scaffold_lengths <- read.table("PAG_UKon_Phar_Genome_scaffold_lengths.txt.header", header=TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

# calculate cumulative position for SNP (in bp and in Mb)
ann.62$BPcumBP <- ann.62$POS + as.numeric(offsets[as.character(ann.62$CHR)])
ann.62$BPcum_Mb <- ann.62$BPcumBP / 1e6

# create unique SNP ID
ann.62$snpid <- paste0(ann.62$CHR, "_", ann.62$BPcumBP)

xpehh.out <- read.table("top.1pc.normxpehh.PAG.outliers.header.chr.SNPids", header = TRUE, sep="", stringsAsFactors=FALSE)
fst.out <- read.table("PHAR_Fst_ALL.OUTLIERS_0.1p_with_SNPIDs.txt", header=TRUE, sep="", stringsAsFactors=FALSE)

# filter for contig 29
xpehh.out.62 <- xpehh.out %>% filter(chr == 62)
fst.out.62 <- fst.out %>% filter(CHROM == "Phar_scaff_62")

ovl_xa_62 <- intersect(as.character(unique(xpehh.out.62$snpid)), as.character(unique(ann.62$snpid)))
ovl_xfa_62 <- intersect(as.character(unique(fst.out.62$snpid)), ovl_xa_62)

print(ovl_xa_62) # 2531
print(ovl_xfa_62) # 49

ann_62_OUT <- ann.62[ann.62$snpid %in% ovl_xfa_62, ]
write.table(ann_62_OUT, "SNPEff.OUTPUT.Phar_scaff_62.ANN-OUT_OUTLIER-SNPs.txt", , sep="\t", quote=FALSE, row.names=FALSE)

# export as readable

# CONTIG 29
data <- read.delim("SNPEff.OUTPUT.Phar_scaff_29.ANN-OUT_OUTLIER-SNPs.txt", sep = "\t", stringsAsFactors = FALSE)
data_sep <- data %>% separate(ANN, into = paste0("ANN_part", 1:61), sep = "\\|", fill = "right", extra = "merge")

write.table(data_sep, file = "SNPEff.OUTPUT.Phar_scaff_29.ANN-OUT_OUTLIER-SNPs_CLEAN.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# CONTIG 62
data <- read.delim("SNPEff.OUTPUT.Phar_scaff_62.ANN-OUT_OUTLIER-SNPs.txt", sep = "\t", stringsAsFactors = FALSE)
data_sep <- data %>% separate(ANN, into = paste0("ANN_part", 1:80), sep = "\\|", fill = "right", extra = "merge")

write.table(data_sep, file = "SNPEff.OUTPUT.Phar_scaff_62.ANN-OUT_OUTLIER-SNPs_CLEAN.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## make nicer output table with these fields: https://pcingola.github.io/SnpEff/snpeff/inputoutput/
