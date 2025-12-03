#!/usr/bin/env Rscript

# Script to re-plot Stairway Plot 

# ------- STAIRWAY PLOT -- EFFECTIVE POPULATION SIZE ------- #

setwd("/path/to/stairway_output/SASY_g34_m483")
library(ggplot2)

stair = read.table("SASY_fold_g34_m483_STAIRWAY.final.summary", header = TRUE)

g = ggplot(stair, aes(x = year)) +
  geom_line(aes(y = Ne_median), color = "red", lwd = 2) +
  geom_ribbon(aes(ymin = Ne_12.5., ymax = Ne_87.5.), alpha = 0.2) + 
  geom_ribbon(aes(ymin = Ne_2.5., ymax = Ne_97.5.), alpha = 0.2) + 
  scale_x_continuous(trans = 'log10', limits = c(1e+3, 2e+5), breaks = c(2e+3, 4e+3, 6e+3, 8e+3, 1e+4, 2e+4, 4e+4, 6e+4, 8e+4, 1e+5, 2e+5), labels = c("2k", "4k", "6k", "8k", "10k", "20k", "40k", "60k", "80k", "100k", "200k")) +
  scale_y_continuous(trans = 'log10', breaks = c(1e+3, 1e+4, 1e+5, 1e+6), labels = c("1k", "10k", "100k", "1M"))  +
  theme_bw() +
  xlab("Years Before Present") +
  ylab("Effective Population Size")
ggsave("SASY_PAG_STAIRWAY_PLOT_m483_g34_replotted.png", g)
ggsave("SASY_PAG_STAIRWAY_PLOT_m483_g34_replotted.pdf", g)
