#!/usr/bin/env Rscript

# Script to plot population-wide Weir & Cockerham's Fst as a heatmap

# ----------------- FST HEATMAP ----------------- #

# load libraries
library(ggplot2)
library(reshape2)

# WEIR & COCKERHAM weighted Fst pairwise population Fst

setwd("/path/to/Fst")

# SA vs SI = 0.040977
# SA vs SY = 0.021426
# SY vs SI = 0.027077

# define pairwise Fst values
fst_matrix <- matrix(
  c(NA, 0.021426, 0.040977,
  0.021426, NA, 0.027077,
  0.040977, 0.027077, NA),
  nrow = 3, byrow = TRUE
)

# assign row and column names in desired order
rownames(fst_matrix) <- c("SA", "SY", "SI")
colnames(fst_matrix) <- c("SA", "SY", "SI")

# melt to long format for ggplot
fst_long <- reshape2::melt(fst_matrix, na.rm = TRUE)
fst_long <- fst_long[as.numeric(fst_long$Var1) > as.numeric(fst_long$Var2), ]

# plot heatmap
ggplot(fst_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", value)), color = "black", size = 8) +
  scale_fill_gradient(low = "#a0dde6", high = "#247b7b", name = "Fst") +
  theme_minimal(base_size = 20) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "") +
  coord_fixed()
ggsave("Fst_heatmap.pdf", dpi = 300)

