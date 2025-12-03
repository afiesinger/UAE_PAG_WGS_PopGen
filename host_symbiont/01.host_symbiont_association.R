#!/usr/bin/env Rscript

# Script to plot ITS2 type profiles in order of dendrogram & PCAs by dominant symbiont type per coral colony & conduct Mantel tests for host vs. symbiont genetic distance

# ----------- PACKAGES & WORKDIR ----------- #

library(tidyverse)
library(ggplot2)
library(reshape2)
library(readxl)
library(dplyr)
library(readr)

setwd("/path/to/output_plots")

# -------------- PCA WITH COLORS ACCORDING TO DOMINANT SYMBIONT TYPE --------------- #

symbiont_data <- as.data.frame(read_excel("/path/to/ITS2_abund.xlsx", sheet = 1)) 

symbiont_long <- symbiont_data %>%
  mutate(across(-c(SampleID, Population), ~as.numeric(gsub(",", ".", .x)))) %>%
  pivot_longer(
    cols = -c(SampleID, Population),
    names_to = "Profiles",
    values_to = "Proportion"
  ) %>%
  select(SampleID, Population, Profiles, Proportion)

# extract dominant Symbiont assignment
dominant_symbiont <- symbiont_long %>%
  group_by(SampleID, Population) %>%
  filter(Proportion == max(Proportion)) %>%
  ungroup() %>%
  select(SampleID, Population, DominantSymbiont = Profiles)

# NOW LOAD PCA DATA FROM PLINK # 

pca_indir="/path/to/PCA_data"

# read the PCA results
pca_data <- read_csv(paste0(pca_indir, "PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.PCA.csv"), col_names = c("FID", "IID", "PC1", "PC2", "PC3", "PC4"))

# merge PCA data with dominant_symbiont metadata
pca_merged <- pca_data %>% left_join(dominant_symbiont, by = c("IID" = "SampleID"))

# calculate the proportion of variance explained
eigenvalues <- scan(paste0(pca_indir, "PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.PCA.eigenval"))
total_variance <- sum(eigenvalues)
variance_explained <- eigenvalues / total_variance * 100

pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)

# COLORS FOR PLOTTING #
unique(pca_merged$DominantSymbiont)
# "C3-C3gulf-C3cc" "C3-C3cc-C3gulf-C3ye" "C3-C3gulf-C3ye" "C3/C3c-C3gulf" "A1" "D1/D4-D4c-D1h" "C3yd/C3yc-C3-C3gulf" "C3/C3gulf" "C15/C116" "C3/C3gulf-C115d"  

# CUSTOM COLORS
colors = c(
  "A1" = "#FF6100",
  "C3-C3gulf-C3cc" = "#d75c67",
  "C3-C3cc-C3gulf-C3ye" = "#a01427",
  "C3-C3gulf-C3ye" = "#e0acd0",
  "C3/C3c-C3gulf" = "#e57f68",
  "C3yd/C3yc-C3-C3gulf" = "#d1c5c5", 
  "C3/C3gulf" = "#fd7266",  
  "C3/C3gulf-C115d" = "#ff0000",
  "C15/C116" = "#8a2ae3",
  "D1/D4-D4c-D1h" = "#097d6cd0")

# ensure Population is in desired order and relabel SI -> AA
pca_merged$Population <- fct_recode(pca_merged$Population, "AA" = "SI")
pca_merged$Population <- factor(pca_merged$Population, 
                                levels = c("SA", "SY", "AA"))

# plotting: color by DominantSymbiont, shape by Population
pca_plot <- ggplot(pca_merged, aes(x = PC1, y = PC2, color = DominantSymbiont, shape = Population)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 6) +
  theme(
    legend.position = "right", 
    legend.text  = element_text(size = 5),  
    legend.title = element_text(size = 5),
    axis.text = element_text(size = 8, face = "bold", color = "black"),
    axis.title = element_text(size = 10, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    legend.key.size = unit(0.4, "cm"),        
    legend.spacing.y = unit(0.01, "cm"),     
    legend.box.spacing = unit(0.3, "cm")
    ) +
  labs(title = "", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"), color = "Majority ITS2 type profile", shape = "Population" ) +
  scale_shape_manual(values = c("AA" = 16, "SY" = 17, "SA" = 15)) +
  scale_color_manual(values = colors) 
ggsave("HOST_SYM_PCA.pdf", pca_plot, width = 100, height = 75, units = "mm")

# --------- SEPARATE PCA FOR EACH DOMINANT ITS2 TYPE PROFILE AND ALL SAMPLES THAT HAVE THIS PROFILE --------- # 

# take pca_merged as input (from above)

library(ggplot2)
library(forcats)
library(scales)

# update shapes and colors, relabeling SI -> AA
pop_shapes <- c("SA" = 15, "SY" = 17, "AA" = 16)
#pop_colors <- c("SA" = "#e3403b", "SY" = "#e8744a", "AA" = "#eda531")

# plot with facets for DominantSymbiont

p <- ggplot(pca_merged, aes(x = PC1, y = PC2, shape = Population)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_shape_manual(values = pop_shapes) +
  scale_color_manual(values = "black") +
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_bw(base_size = 12) +   
  labs(
    x = "PC1",
    y = "PC2",
    shape = "Population"
  ) +
  facet_wrap(~ DominantSymbiont, scales = "free", nrow = 3) +
  theme(
    strip.text = element_text(size = 5, face = "bold"),
    axis.text = element_text(size = 4, color = "black"),
    axis.title = element_text(size = 5, color = "black", face = "bold"),
    legend.position = "none"
)
ggsave("PCA_DOM_SYMBIONT_facet.pdf", p, width = 130, height = 100, units = "mm")


# --------------- ORDERED BARPLOT OF ITS2 TYPE PROFILES TO BE PLOTTED IN THE DENDROGRAM OF HOST GENOTYPES ------------- #

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(scales)

cols_2022 <- c(
  "A1" = "#FF6100", 
  "A1-A13a" = "#FA8F11", 
  "A1-A1bv" = "#F5C372", 
  "A1-A1bw" = "#F6E496", 
  "A1-A1bw-A1bf-A1bx-A1eb" = "#fcaf58", 
  "A1-A1bw-A1bx" = "#E7C703", 
  "A1-A1eq-A1ep" = "#FFFF7A", 
  
  "C15" = "#800080",
  "C15h" = "#8a2be2",
  "C15h-C15o-C15k" = "#b66ee8", 
  "C15/C116" = "#cd34b5",
  "C3-C3cc-C3gulf-C3ye" = "#a01427", 
  "C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak" = "#cc2936", 
  "C3-C3gulf-C3cc" = "#c75c67", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak" = "#e5989b", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ak" = "#ffbfc0", 
  "C3-C3gulf-C3d-C3i-C115c-C3ak-C3al" = "#ffe8f0", 
  "C3-C3gulf-C3ye" = "#e0acd0", 
  "C3/C15h-C3gulf" = "#bf69a2", 
  "C3/C3c-C3gulf" = "#e57f68", 
  "C3/C3gulf" = "#fd7266",
  "C3/C3gulf-C115d" = "#ff0000",
  "C39-C39q-C1-C39b-C39p" = "#a13b42",
  "C3yd" = "#bd9092",
  "C3yd/C3yc-C3-C3gulf" = "#d1c5c5",  
  
  "D1" = "#000080", 
  "D1-D4-D4c-D17ao-D17ap-D17aq" = "#114477",
  "D1-D4-D4c-D17ap-D17aq-D17ao-D17ar" = "#4169e1",
  "D1-D4-D4c-D17d" = "#75a4f8",
  "D1-D4-D4c-D17d-D1r-D17c-D17e" = "#1ab2e5b9",  
  "D1-D4-D4c-D17d-D1r-D17c-D17e-D17j" = "#5e8ebf",
  "D1/D2" = "#5f9ea0", 
  "D1/D2-D4-D4c-D1r" = "#57d4d5", 
  "D1/D4-D4c-D1h" = "#097d6cd0", 
  "D1/D4-D4c-D1r" = "#bcf6f5", 
  "D1/D4/D2-D4c-D1c-D1h" = "#b5d0ff", 
  "D5-D5a-D4-D4a-D4b" = "#bcdbdb", 
  "D5-D5a-D5ai-D4-D5f-D5v" = "#81c784", 
  "D5-D5a-D5f" = "#44aa77", 
  "D5-D5c-D4a-D5b-D4-D5i-D5a" = "#558b2f", 
  "D5/D5a-D4-D4a-D2" = "#1ca62f", 
  "D5a-D5-D5ah-D4" = "#006400"
)

setwd("/path/to/output_plots")

# ITS2 type profiles
its2Profs <- read.delim("396_20230721T100004_DBV_20230721T154120.profiles.relative.abund_and_meta.txt", 
                 skip = 6,    
                 header = TRUE,
                 check.names = FALSE)[ , -1]
colnames(its2Profs)[1] = "Genotype"
head(its2Profs)

# add sample names, species and site from metadata file
its2MetaData <- read_excel("its2_profiles_metadata.xlsx")

its2Profs <- its2Profs[match(its2MetaData$sample, its2Profs$Genotype), ]
all(its2Profs$Genotype == its2MetaData$sample)

its2Profs <- cbind(its2MetaData[, c(6,8)], its2Profs) 
colnames(its2Profs)[1:2] <- c("Species", "Site")
head(its2Profs)
            
# filter only Phar (full file has another species in it (Pdae))
its2Profs_Phar <- its2Profs %>% filter(Species == "Phar")

# SET SAMPLE ORDER AS IN DENDROGRAM 

samples <- c(
  "UAE_SY_Phar_27_field_PopGen_MetaB", "UAE_SY_Phar_39_field_PopGen_MetaB", "UAE_SY_Phar_40_field_PopGen_MetaB", 
  "UAE_SY_Phar_3_field_PopGen_MetaB", "UAE_SY_Phar_34_field_PopGen_MetaB", "UAE_SY_Phar_12_field_PopGen_MetaB", 
  "UAE_SY_Phar_17_field_PopGen_MetaB", "UAE_SY_Phar_32_field_PopGen_MetaB", "UAE_SY_Phar_21_field_PopGen_MetaB", 
  "UAE_SY_Phar_25_field_PopGen_MetaB", "UAE_SY_Phar_11_field_PopGen_MetaB", "UAE_SY_Phar_20_field_PopGen_MetaB", 
  "UAE_SY_Phar_22_field_PopGen_MetaB", "UAE_SY_Phar_1_field_PopGen_MetaB", "UAE_SY_Phar_28_field_PopGen_MetaB", 
  "UAE_SY_Phar_35_field_PopGen_MetaB", "UAE_SY_Phar_38_field_PopGen_MetaB", "UAE_SY_Phar_5_field_PopGen_MetaB", 
  "UAE_SY_Phar_4_field_PopGen_MetaB", "UAE_SY_Phar_9_field_PopGen_MetaB", "UAE_SY_Phar_7_field_PopGen_MetaB", 
  "UAE_SY_Phar_18_field_PopGen_MetaB", "UAE_SY_Phar_2_field_PopGen_MetaB", "UAE_SY_Phar_8_field_PopGen_MetaB", 
  "UAE_SY_Phar_13_field_PopGen_MetaB", "UAE_SY_Phar_16_field_PopGen_MetaB", "UAE_SY_Phar_6_field_PopGen_MetaB", 
  "UAE_SY_Phar_33_field_PopGen_MetaB", "UAE_SA_Phar_48_field_PopGen_MetaB", "UAE_SA_Phar_50_field_PopGen_MetaB", 
  "UAE_SA_Phar_43_field_PopGen_MetaB", "UAE_SA_Phar_17_field_PopGen_MetaB", "UAE_SA_Phar_44_field_PopGen_MetaB", 
  "UAE_SA_Phar_42_field_PopGen_MetaB", "UAE_SA_Phar_47_field_PopGen_MetaB", "UAE_SA_Phar_49_field_PopGen_MetaB", 
  "UAE_SA_Phar_23_field_PopGen_MetaB", "UAE_SA_Phar_21_field_PopGen_MetaB", "UAE_SA_Phar_31_field_PopGen_MetaB", 
  "UAE_SA_Phar_41_field_PopGen_MetaB", "UAE_SA_Phar_36_field_PopGen_MetaB", "UAE_SA_Phar_39_field_PopGen_MetaB", 
  "UAE_SA_Phar_12_field_PopGen_MetaB", "UAE_SA_Phar_40_field_PopGen_MetaB", "UAE_SI_Phar_27_field_PopGen_MetaB", 
  "UAE_SA_Phar_15_field_PopGen_MetaB", "UAE_SA_Phar_30_field_PopGen_MetaB", "UAE_SA_Phar_24_field_PopGen_MetaB", 
  "UAE_SA_Phar_35_field_PopGen_MetaB", "UAE_SA_Phar_33_field_PopGen_MetaB", "UAE_SA_Phar_34_field_PopGen_MetaB", 
  "UAE_SA_Phar_13_field_PopGen_MetaB", "UAE_SA_Phar_32_field_PopGen_MetaB", "UAE_SA_Phar_29_field_PopGen_MetaB", 
  "UAE_SA_Phar_25_field_PopGen_MetaB", "UAE_SA_Phar_14_field_PopGen_MetaB", "UAE_SA_Phar_18_field_PopGen_MetaB", 
  "UAE_SA_Phar_16_field_PopGen_MetaB", "UAE_SA_Phar_26_field_PopGen_MetaB", "UAE_SA_Phar_45_field_PopGen_MetaB", 
  "UAE_SI_Phar_40_field_PopGen_MetaB", "UAE_SI_Phar_35_field_PopGen_MetaB", "UAE_SI_Phar_38_field_PopGen_MetaB", 
  "UAE_SI_Phar_39_field_PopGen_MetaB", "UAE_SI_Phar_11_field_PopGen_MetaB", "UAE_SI_Phar_20_field_PopGen_MetaB", 
  "UAE_SI_Phar_19_field_PopGen_MetaB", "UAE_SI_Phar_4_field_PopGen_MetaB", "UAE_SI_Phar_8_field_PopGen_MetaB", 
  "UAE_SI_Phar_25_field_PopGen_MetaB", "UAE_SI_Phar_22_field_PopGen_MetaB", "UAE_SI_Phar_2_field_PopGen_MetaB", 
  "UAE_SI_Phar_3_field_PopGen_MetaB", "UAE_SI_Phar_23_field_PopGen_MetaB", "UAE_SI_Phar_7_field_PopGen_MetaB", 
  "UAE_SI_Phar_32_field_PopGen_MetaB", "UAE_SI_Phar_15_field_PopGen_MetaB", "UAE_SI_Phar_10_field_PopGen_MetaB", 
  "UAE_SI_Phar_14_field_PopGen_MetaB", "UAE_SI_Phar_21_field_PopGen_MetaB", "UAE_SI_Phar_9_field_PopGen_MetaB", 
  "UAE_SI_Phar_12_field_PopGen_MetaB", "UAE_SI_Phar_16_field_PopGen_MetaB", "UAE_SY_Phar_24_field_PopGen_MetaB", 
  "UAE_SY_Phar_15_field_PopGen_MetaB", "UAE_SY_Phar_37_field_PopGen_MetaB", "UAE_SY_Phar_19_field_PopGen_MetaB", 
  "UAE_SY_Phar_30_field_PopGen_MetaB", "UAE_SY_Phar_23_field_PopGen_MetaB", "UAE_SY_Phar_36_field_PopGen_MetaB", 
  "UAE_SY_Phar_29_field_PopGen_MetaB", "UAE_SY_Phar_31_field_PopGen_MetaB", "UAE_SY_Phar_14_field_PopGen_MetaB", 
  "UAE_SY_Phar_26_field_PopGen_MetaB"
)

common_samples <- intersect(its2Profs_Phar$Genotype, samples)
its2Profs_Phar <- its2Profs_Phar[match(common_samples, its2Profs_Phar$Genotype), ]

# pivot to long format
its2Profs_long <- its2Profs_Phar %>%
  pivot_longer(cols = -c(Genotype, Species, Site), 
               names_to = "Profile", 
               values_to = "Abundance")

# ensure numeric
its2Profs_long$Abundance <- as.numeric(its2Profs_long$Abundance)

its2Profs_long$Genotype <- factor(its2Profs_long$Genotype, 
                                   levels = samples)

# FILTER OUT NAs 
profiles_to_remove <- its2Profs_long %>%
  group_by(Profile) %>%
  summarize(all_zero_na = all(is.na(Abundance) | Abundance == 0)) %>%
  filter(all_zero_na) %>%
  pull(Profile)

its2Profs_long_filt <- its2Profs_long %>%
  filter(!Profile %in% profiles_to_remove)

### PLOT ALL TOGETHER ###

ggplot(its2Profs_long_filt, aes(x = Genotype, y = Abundance, fill = Profile)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = cols_2022) +
  labs(y = "", x = "") +
  theme_minimal() +
  theme(
  axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6, color = "black"),
  axis.text.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank())
ggsave("Barplot_for_dendrogram_minimal.pdf", width = 20, height = 6, dpi = 300)

# ---------------- MANTEL TEST HOST VS SYMBIONT PAIRWISE GENETIC DISTANCE --------------- #

library(ggpubr)
library(scales)
library(vegan)
library(vcfR)

setwd("/path/to/files")

# READ HOST GENOTYPE DATA
vcf <- read.vcfR("/path/to/vcf/PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.recode.vcf") 
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

# READ SYMBIONT GENUS-SPECIFIC DIST MATRICES

symb_A_dist <- read.table("/home/fiesingera/proj/PAG_Phar_PopGen/NEW_FILTERS_VCF/Host_Symbiont_Genotyping/20230721T154120_braycurtis_sample_distances_A_no_sqrt.dist", header = FALSE)
symb_C_dist <- read.table("/home/fiesingera/proj/PAG_Phar_PopGen/NEW_FILTERS_VCF/Host_Symbiont_Genotyping/20230721T154120_braycurtis_sample_distances_C_no_sqrt.dist", header = FALSE)
symb_D_dist <- read.table("/home/fiesingera/proj/PAG_Phar_PopGen/NEW_FILTERS_VCF/Host_Symbiont_Genotyping/20230721T154120_braycurtis_sample_distances_D_no_sqrt.dist", header = FALSE)

# remove unnecessary column (your V2) and all Pdae samples
symb_A_dist <- symb_A_dist[,-2]
symb_C_dist <- symb_C_dist[,-2]
symb_D_dist <- symb_D_dist[,-2]

colnames(symb_A_dist) <- c("SampleID", symb_A_dist$V1)
colnames(symb_C_dist) <- c("SampleID", symb_C_dist$V1)
colnames(symb_D_dist) <- c("SampleID", symb_D_dist$V1)

host_samples <- colnames(gt_matrix)

a <- intersect(symb_A_dist$SampleID, host_samples)
filtered_rows <- symb_A_dist[symb_A_dist$SampleID %in% a, ]
rownames(filtered_rows) <- filtered_rows$SampleID
filt_A <- filtered_rows[, c("SampleID", a)]

dim(filt_A)

c <- intersect(symb_C_dist$SampleID, host_samples)
filtered_rows <- symb_C_dist[symb_C_dist$SampleID %in% c, ]
rownames(filtered_rows) <- filtered_rows$SampleID
filt_C <- filtered_rows[, c("SampleID", c)]

dim(filt_C)

d <- intersect(symb_D_dist$SampleID, host_samples)
filtered_rows <- symb_D_dist[symb_D_dist$SampleID %in% d, ]
rownames(filtered_rows) <- filtered_rows$SampleID
filt_D <- filtered_rows[, c("SampleID", d)]

dim(filt_D)

# Convert distance blocks to matrices (drop V1 column, ensure symmetry)
symb_A_mat <- as.matrix(filt_A[,-1])
symb_C_mat <- as.matrix(filt_C[,-1])
symb_D_mat <- as.matrix(filt_D[,-1])

rownames(symb_A_mat) <- a
colnames(symb_A_mat) <- a

rownames(symb_C_mat) <- c
colnames(symb_C_mat) <- c

rownames(symb_D_mat) <- d
colnames(symb_D_mat) <- d

gt_A <- gt_matrix[, colnames(gt_matrix) %in% a]
gt_A <- gt_A[, a]
gt_A_t <- t(gt_A)
host_dist_A <- dist(gt_A_t, method = "euclidean")

symbiont_dist_A <- vegdist(symb_A_mat, method = "bray")

gt_C <- gt_matrix[, colnames(gt_matrix) %in% c]
gt_C <- gt_C[, c]
gt_C_t <- t(gt_C)
host_dist_C <- dist(gt_C_t, method = "euclidean")

symbiont_dist_C <- vegdist(symb_C_mat, method = "bray")

gt_D <- gt_matrix[, colnames(gt_matrix) %in% d]
gt_D <- gt_D[, d]
gt_D_t <- t(gt_D)
host_dist_D <- dist(gt_D_t, method = "euclidean")

symbiont_dist_D <- vegdist(symb_D_mat, method = "bray")

# MANTEL TESTS
mantel_result_A <- mantel(host_dist_A, symbiont_dist_A, method = "pearson", permutations = 9999)
mantel_result_C <- mantel(host_dist_C, symbiont_dist_C, method = "pearson", permutations = 9999)
mantel_result_D <- mantel(host_dist_D, symbiont_dist_D, method = "pearson", permutations = 9999)

print(mantel_result_A)
print(mantel_result_C)
print(mantel_result_D)

