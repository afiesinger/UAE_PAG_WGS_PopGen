#!/usr/bin/env Rscript

# Script to plot PCAs 

# ----------------- PCA NO CLONES ----------------- #

# load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)

# working directory 
setwd("/path/to/vcf/")

# set the output prefix
OUTPUT_PREFIX <- "PHAR.OUT.SNPs_GATK_filt_hwe_maf_bi_indv.ld_pruned.NOCLONES.PCA"

# read the PCA results
pca_data <- read_csv(paste0(OUTPUT_PREFIX, ".csv"), col_names = c("FID", "IID", "PC1", "PC2", "PC3", "PC4"))

# extract population information from IID
pca_data <- pca_data %>% mutate(Population = sub("UAE_(SI|SY|SA)_Phar_.*", "\\1", IID))
pca_data$Population <- factor(pca_data$Population, levels = c("SA", "SY", "SI"))

# calculate the proportion of variance explained
eigenvalues <- scan(paste0(OUTPUT_PREFIX, ".eigenval"))
total_variance <- sum(eigenvalues)
variance_explained <- eigenvalues / total_variance * 100

# round the variance explained to two decimal places
pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)

pc3_var <- round(variance_explained[3], 2)
pc4_var <- round(variance_explained[4], 2)

# PC1 & PC2 #

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Population, label = IID)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw(base_size = 20) +
  theme(
    axis.text = element_text(color = "black", size = 15, face = "bold"),
    axis.title = element_text(color = "black", size = 18, face = "bold"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black")
  ) +
  labs(title = "", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"))
ggsave(paste0(OUTPUT_PREFIX, ".pdf"), pca_plot, width = 10, height = 8)
ggsave(paste0(OUTPUT_PREFIX, ".png"), pca_plot, width = 10, height = 8)

# ELLIPSES COVERING THE SCATTER POINTS
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(type = "norm", lty = 2, position = "identity") +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black")
  ) +
  labs(title = "", 
       x = paste0("PC1 (", pc1_var, "%)"), 
       y = paste0("PC2 (", pc2_var, "%)"))
ggsave(paste0(OUTPUT_PREFIX, "ellipse.pdf"), pca_plot, width = 10, height = 8)
ggsave(paste0(OUTPUT_PREFIX, "ellipse.png"), pca_plot, width = 10, height = 8)

# Display summary of the plot
print(summary(pca_data[, c("PC1", "PC2")]))

# PC1 & PC3 #
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC3, color = Population, label = IID)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw(base_size = 16) +
  labs(title = "", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC3 (", pc3_var, "%)"))
ggsave(paste0(OUTPUT_PREFIX, "_PC1PC3.pdf"), pca_plot, width = 10, height = 8)
ggsave(paste0(OUTPUT_PREFIX, "_PC1PC3.png"), pca_plot, width = 10, height = 8)

# PC1 & PC4 #
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC4, color = Population, label = IID)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw(base_size = 16) +
  labs(title = "", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC4 (", pc4_var, "%)"))
ggsave(paste0(OUTPUT_PREFIX, "_PC1PC4.pdf"), pca_plot, width = 10, height = 8)
ggsave(paste0(OUTPUT_PREFIX, "_PC1PC4.png"), pca_plot, width = 10, height = 8)

# PC2 & PC3 #
pca_plot <- ggplot(pca_data, aes(x = PC2, y = PC3, color = Population, label = IID)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw(base_size = 16) +
  labs(title = "", x = paste0("PC2 (", pc2_var, "%)"), y = paste0("PC3 (", pc3_var, "%)"))
ggsave(paste0(OUTPUT_PREFIX, "_PC2PC3.pdf"), pca_plot, width = 10, height = 8)
ggsave(paste0(OUTPUT_PREFIX, "_PC2PC3.png"), pca_plot, width = 10, height = 8)

# PC2 & PC4 #
pca_plot <- ggplot(pca_data, aes(x = PC2, y = PC4, color = Population, label = IID)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw(base_size = 16) +
  labs(title = "", x = paste0("PC2 (", pc2_var, "%)"), y = paste0("PC4 (", pc4_var, "%)"))
ggsave(paste0(OUTPUT_PREFIX, "_PC2PC4.pdf"), pca_plot, width = 10, height = 8)
ggsave(paste0(OUTPUT_PREFIX, "_PC2PC4.png"), pca_plot, width = 10, height = 8)

# PC3 & PC4 #
pca_plot_PC3PC4 <- ggplot(pca_data, aes(x = PC3, y = PC4, color = Population)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw() +
  labs(title = "", x = paste0("PC3 (", pc3_var, "%)"), y = paste0("PC4 (", pc4_var, "%)"))
ggsave(paste0(OUTPUT_PREFIX, "_PC3PC4.pdf"), pca_plot_PC3PC4, width = 10, height = 8)
ggsave(paste0(OUTPUT_PREFIX, "_PC3PC4.png"), pca_plot_PC3PC4, width = 10, height = 8)

# --------- RANDOM RESAMPLING OF SNPs MULTIPLE PCAs ---------- #

# load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)

setwd("/path/to/files/")

# number of replicates
REPS <- 20

# read all PCA replicates and stack as one data frame
all_pca <- map_dfr(1:REPS, function(i) {
  f <- sprintf("snps_subset.%d.PCA.csv", i)
  dat <- read_csv(f, col_names = c("FID", "IID", "PC1", "PC2", "PC3", "PC4"))
  dat$Replicate <- i
  dat$Population <- sub("UAE_(SI|SY|SA)_Phar_.*", "\\1", dat$IID)
  dat$Population <- factor(dat$Population, levels = c("SA", "SY", "SI"))
  dat
})

# plot all replicates in a facet
p <- ggplot(all_pca, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b")) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black")
  ) +
  labs(x = "PC1", y = "PC2") +
  facet_wrap(~ Replicate)
ggsave("PCA_random_resampling.png", p, width = 20, height = 16, dpi = 300)
