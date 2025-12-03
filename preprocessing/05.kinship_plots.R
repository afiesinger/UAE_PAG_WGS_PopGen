#!/usr/bin/env Rscript

# Script to plot KING output and infer kinship for each pair of samples & plot dendrograms of sample genetic similarity

# -------------- KINSHIP INFERENCE ---------------- #

# load packages
library(ggplot2)
library(readr)
library(dplyr)

setwd("/path/to/files")

# read KING .kin0 output
king_file <- "king.kin0_ALL"
kin <- read_table(king_file, col_names = TRUE)

kin <- kin %>% mutate(Kinship_Clean = pmax(Kinship, 0))

ggplot(kin, aes(x = IBS0, y = Kinship_Clean)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.045, linetype = "dashed", color = "grey40") +   # unrelated
  geom_hline(yintercept = 0.09, linetype = "dashed", color = "grey40") +   # third degree
  geom_hline(yintercept = 0.177, linetype = "dashed", color = "grey40") +   # second degree
  geom_hline(yintercept = 0.35, linetype = "dashed", color = "grey40") +   # full sib
  annotate("text", x = max(kin$IBS0), y = 0.35 + 0.01, label = "Twin", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.35 - 0.01, label = "Full sib", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.177 - 0.01, label = "Second Degree", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.09 - 0.01, label = "Third Degree", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.045 - 0.01, label = "Unrelated", hjust = 1, color = "grey40", size = 3) +
  labs(
    x = "Proportion IBS 0",
    y = "Kinship Coefficient"
  ) +
  theme_bw(base_size = 16)
ggsave("king_kinship_vs_ibs0.png", width = 8, height = 5, dpi = 300)

# identify siblings

library(ggrepel)

# filter outliers (twins and full sibs)
outliers <- kin %>% filter(Kinship_Clean > 0.177)

x_max <- max(kin$IBS0)

ggplot(kin, aes(x = IBS0, y = Kinship_Clean)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0.045, linetype = "dashed", color = "grey40") +   # unrelated
  geom_hline(yintercept = 0.09, linetype = "dashed", color = "grey40") +   # third degree
  geom_hline(yintercept = 0.177, linetype = "dashed", color = "grey40") +   # second degree
  geom_hline(yintercept = 0.35, linetype = "dashed", color = "grey40") +   # full sib
  annotate("text", x = max(kin$IBS0), y = 0.35 + 0.01, label = "Twin", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.35 - 0.01, label = "Full sib", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.177 - 0.01, label = "Second Degree", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.09 - 0.01, label = "Third Degree", hjust = 1, color = "grey40", size = 3) +
  annotate("text", x = max(kin$IBS0), y = 0.045 - 0.01, label = "Unrelated", hjust = 1, color = "grey40", size = 3) +
  # Add labels for outliers
  geom_text_repel(
    data = outliers,
    aes(label = paste(ID1, ID2, sep = "-")),
    color = "purple", size = 3, min.segment.length = 0
  ) +
  labs(
    x = "Proportion IBS 0",
    y = "Kinship Coefficient"
  ) +
  theme_bw(base_size = 16)
ggsave("king_kinship_vs_ibs0_OUTLIERS.png", width = 8, height = 5, dpi = 300)

# --------------- DENDROGRAM ----------------- #

# load packages
library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(svglite)

setwd("/path/to/files")

# load IBS matrix
ibs <- as.matrix(read.table("ibs_results.mibs", header=FALSE))

# load sample names
ids <- read.table("ibs_results.mibs.id", header=FALSE, stringsAsFactors=FALSE)[,1]

# assign row and column names for clarity
rownames(ibs) <- ids
colnames(ibs) <- ids

# convert similarity matrix (identitiy by state = IBS)
dist_mat <- 1 - ibs

# convert to 'dist' object for clustering
dd <- as.dist(dist_mat)

svglite::svglite("Phar_dendrogram_ALL.svg", width = 12, height = 10)
plot(hclust(dd, method="average"), labels=ids, main="Hierarchical clustering dendrogram", cex=0.6)
dev.off()

# get order of sample names in dendrogram for plot
hc$order # number
hc$labels # labels unordered
labels_in_order <- hc$labels[hc$order]

# --------------- DENDROGRAM WITHOUT CLONE MATES ----------------- #

# load packages
library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(svglite)

setwd("/path/to/files")

ibs <- as.matrix(read.table("ibs_results_noclones.mibs", header=FALSE))
ids <- read.table("ibs_results_noclones.mibs.id", header=FALSE, stringsAsFactors=FALSE)[,1]

rownames(ibs) <- ids
colnames(ibs) <- ids

dist_mat <- 1 - ibs

dd <- as.dist(dist_mat)

svglite::svglite("Phar_dendrogram_NOCLONES.svg", width = 12, height = 10)
plot(hclust(dd, method="average"), labels=ids, main="Hierarchical clustering dendrogram", cex=0.6)
dev.off()

