#!/usr/bin/env Rscript

# Script to plot cross-validation error from ADMIXTURE runs and CLUMPP output

# ----------------- CV ERRORS ----------------- #

setwd("/path/to/ADMIXTURE")

# load libraries
library(dplyr)
library(ggplot2)

# read data
cv <- read.table("cv_error_table.txt", header = FALSE, col.names = c("K", "CVerror"))

# summarize
cv_summary <- cv %>% group_by(K) %>% summarise(mean_CV = mean(CVerror), sd_CV = sd(CVerror))

# plot with error bars
plot = ggplot(cv_summary, aes(x = K, y = mean_CV)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_CV - sd_CV, ymax = mean_CV + sd_CV), width = 0.2) +
  scale_x_continuous(breaks = 1:10, labels = as.character(1:10)) +
  theme_bw(base_size = 20) +
  theme(
    axis.text = element_text(color = "black", size = 15, face = "bold"),
    axis.title = element_text(color = "black", size = 18, face = "bold"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black")
  ) +
  labs(x = "Number of clusters (K)", y = "Cross-validation error") 
ggsave("CV_errors_ADMIXTURE.pdf", plot = plot, dpi = 300)

# -------------- MAKE CLUMPP INPUT FILES FROM ADMIXTURE OUTPUT -------------- #

setwd("/path/to/ADMIXTURE")

# load libraries
library(pophelper)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# read replicate Q files (multiple runs)
qfiles <- list.files(path = "/path/to/ADMIXTURE/", pattern = "*.Q", full.names = TRUE)
clumppExport(readQ(qfiles), exportpath="CLUMPP")

# now back to 03.ADMIXTURE.sh
# run CLUMPP, move all "-combined-merged.txt" files for each K into CLUMPP parent directory & return to R for plotting

# ------------------- PLOT CLUMPP ----------------- #

setwd("/path/to/ADMIXTURE")

clumpp_dir <- "/path/to/CLUMPP"

# generate file names for K = 2:4
clumpp_files <- file.path(clumpp_dir, paste0("pop_K", 2:4, "-combined-merged.txt"))

# read the Q-matrices from CLUMPP output files
qlist <- readQ(clumpp_files)

# align clusters within each K
qlist_aligned <- alignK(qlist)

# read in population file (tab-delimited; assigns each sample a population)
grplab = read.table("Phar_ADMIX.popfile")
custom_cols = c("#de8a5a", "#930000", "#F1E8B0", "#B0CCA2", "#70a494")

plot_list <- plotQ(qlist_aligned, clustercol = custom_cols, grplab = grplab, sortind="Cluster1", sharedindlab = FALSE, imgoutput = "join", subsetgrp=c("SA", "SY", "SI"), returnplot = TRUE, showgrplab = TRUE, outputfilename="CLUMPP_NEW_K2-5", exportpath = getwd())

# save as pdf
plot_list <- plotQ(qlist_aligned, clustercol = custom_cols, grplab = grplab, sortind="Cluster1", sharedindlab = FALSE, imgoutput = "join", subsetgrp=c("SA", "SY", "SI"), returnplot = TRUE, showgrplab = TRUE, outputfilename="CLUMPP_NEW_K2-5", imgtype="pdf", exportpath = getwd())

# save as svg
plot_list <- plotQ(qlist_aligned, clustercol = custom_cols, grplab = grplab, sortind="Cluster1", sharedindlab = FALSE, imgoutput = "join", subsetgrp=c("SA", "SY", "SI"), returnplot = TRUE, showgrplab = TRUE, exportplot=FALSE)
g = grid.arrange(plot_list$plot[[1]])
ggsave("CLUMPP_K2-4.svg", g)


