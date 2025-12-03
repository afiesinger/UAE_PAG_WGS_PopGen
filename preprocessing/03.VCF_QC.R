#!/usr/bin/env Rscript

# -------------------- QC --------------------- #

# load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)

setwd("/path/to/files/")

# read the .imiss file output from vcftools
imiss <- read.table("out.imiss", header = TRUE)

# extract population code from the INDV name using regex
# this assumes names are formatted like this: UAE_SA_Phar_11, UAE_SY_Phar_12, etc.
imiss <- imiss %>% mutate(Population = str_extract(INDV, "(?<=UAE_)[A-Z]{2}(?=_)"))

# check that populations were correctly extracted
table(imiss$Population)

# plot missingness across all individuals
p1 <- ggplot(imiss, aes(x = reorder(INDV, F_MISS), y = F_MISS, fill = Population)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b"))) +
  labs(x = "Individuals", y = "Proportion of Missing Genotypes (F_MISS)",
       title = "Per-Individual Missingness Across All Samples") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

# plot missingness distribution by population (boxplot)
p2 <- ggplot(imiss, aes(x = Population, y = F_MISS, fill = Population)) +
  geom_boxplot(alpha = 0.8, outlier.color = NA) +
  scale_fill_manual(values = c(c("SI" = "#eda531", "SY" = "#e8744a", "SA" = "#e3403b"))) +
  theme_bw(base_size = 14) +
  labs(x = "Population", y = "Proportion of Missing Genotypes (F_MISS)",
       title = "Distribution of Missingness per Population") +
  theme(legend.position = "none")

g = ggarrange(p1, p2, ncol = 2)
ggsave("Missingness_per_sample.png", g, width = 20, height = 12, dpi = 300)

# missingness with different thresholds (here 20% and 30%) as table
imiss %>% 
  group_by(Population) %>% 
  summarise(
    retained_20 = sum(F_MISS <= 0.2),
    retained_30 = sum(F_MISS <= 0.3),
    total = n()
  ) %>%
  mutate(prop_20 = retained_20/total,
         prop_30 = retained_30/total)

