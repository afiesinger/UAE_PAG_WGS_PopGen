#!/usr/bin/env Rscript 

# ----------- PACKAGES -------------- #

library(ggplot2)
library(dplyr)
library(readr)
library(readxl)
library(tidyr)

# ---------------- BASELINE TEMP PHAR SA SY SI ----------------- #

setwd("UAE_PAG_WGS_PopGen/cbass/1_baseline/")

# read in TPM data
tpm_df <- read_tsv("results/star_salmon/salmon.merged.gene_tpm.tsv", col_types = cols())[,-2]

# read in annotations file (from Porites harrisoni reference genome) and clean-up gene names
# Porites harrisoni reference genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040938025.1/

annotations <- read_excel("PAG_UKon_Phar.V2.annotations.xlsx") %>%
  select(GeneID, Name) %>%
  mutate(Name = ifelse(is.na(Name), "hypoth_prot", Name))

genes_of_interest <- c("ACVR2B", "BMPR2", "CDK2", "COL4A1_2", "COL2A1_1", "ORC4", "ANKS6_1", "ANO10_1", "DCAF11", "DYNC1LI1", "eva-1", "FBXO31", "HNF4A", "KDELC1_1", "PDCD6_1", "PDCD6_2", "PDRG1", "RAB22A", "RBM18_1", "RSPH10B", "SSNA1", "TTLL9", "TUBB4B", "ZDHHC7_1")
keep_ids <- annotations$GeneID[annotations$Name %in% genes_of_interest]

tpm_subset <- tpm_df %>% filter(gene_id %in% keep_ids)

tpm_subset <- tpm_subset %>%
  left_join(annotations, by = c("gene_id" = "GeneID")) %>%
  filter(Name %in% genes_of_interest)

tpm_long <- tpm_subset %>%
  pivot_longer(
    cols = -c(gene_id, Name),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  select(Sample, Gene = Name, TPM)

tpm_long$Population <- ifelse(grepl("_SI_", tpm_long$Sample), "GO", "PAG")

xvalues = c("ACVR2B", "BMPR2", "CDK2", "COL4A1", "COL2A1", "ORC4", "ANKS6", "ANO10", "DCAF11", "DYNC1LI1", "eva-1", "FBXO31", "HNF4A", "KDELC1", "PDCD6", "PDRG1", "RAB22A", "RBM18", "RSPH10B", "SSNA1", "TTLL9", "TUBB4B", "ZDHHC7")

ggplot(tpm_long_filtered, aes(x = Gene, y = TPM, fill = Population, color = Population)) +
  geom_boxplot(alpha = 0.7) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  scale_fill_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  scale_color_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  scale_x_discrete(labels = xvalues) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
  ) +
  labs(x = "Gene", y = "TPM", title = "")
ggsave("Genes_of_interest_TPM_boxplot_temp1.pdf", width = 20, height = 10, dpi = 300)

# SPLIT BY T1 & T2 (360 and 1080) #

# take all files from above

pattern <- paste0("^(", paste(xvalues, collapse = "|"), ")")

# keep only annotation rows whose Name starts with any of xvalues (e.g. PDCD6_1, PDCD6_2)
annotations_sub <- annot_filt %>%
  filter(str_detect(Name, pattern))

tpm_subset <- tpm_df %>%
  filter(gene_id %in% annotations_sub$GeneID) %>%
  left_join(annotations_sub, by = c("gene_id" = "GeneID"))

tpm_subset <- tpm_subset %>%
  mutate(
    Gene_root = str_extract(
      Name,
      "^[^_]+"
    )
  ) %>%
  filter(Gene_root %in% xvalues)

tpm_long <- tpm_subset %>%
  pivot_longer(
    cols = -c(gene_id, Name, Gene_root),
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  select(sample, Gene = Gene_root, TPM)
print(tpm_long)

sample_table <- tpm_long %>%
  distinct(sample) %>%
  mutate(
    population = case_when(
      str_detect(sample, "_SY_") ~ "PAG",
      str_detect(sample, "_SA_") ~ "PAG",
      str_detect(sample, "_SI_") ~ "GO",
      TRUE ~ NA_character_
    ),
    timepoint = case_when(
      str_detect(sample, "_360$")  ~ "T1",
      str_detect(sample, "_1080$") ~ "T2",
      TRUE ~ NA_character_
    )
  )
print(sample_table)

tpm_long_meta <- tpm_long %>%
  left_join(sample_table, by = "sample")
print(tpm_long_meta)

tpm_long_filtered <- tpm_long_meta %>%
  filter(TPM != max(TPM, na.rm = TRUE))

gene_levels <- xvalues

df_plot <- tpm_long_filtered %>%
  mutate(
    timepoint  = factor(timepoint,  levels = c("T1", "T2")),
    population = factor(population, levels = c("PAG", "GO")),
    Gene       = factor(Gene,      levels = gene_levels),
    Gene_tp    = paste(Gene, timepoint, sep = "_"),
    Gene_tp    = factor(
      Gene_tp,
      levels = as.vector(rbind(
        paste0(gene_levels, "_T1"),
        paste0(gene_levels, "_T2")
      ))
    )
  )
print(df_plot)

ggplot(df_plot, aes(x = Gene_tp, y = TPM, fill = population)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_jitter(
    aes(color = population),
    position = position_jitterdodge(jitter.width = 0.2),
    size = 1.5,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  scale_color_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  labs(x = "", y = "TPM") +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(
      size = 14, color = "black", face = "bold",
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title  = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks  = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid   = element_blank()
  ) +
  scale_x_discrete(
    labels = function(x) sub("_(T1|T2)$", "\n\\1", x)
)
ggsave("Genes_of_interest_TPM_boxplot_temp1_BY_TIMEPOINTS.pdf", width = 35, height = 15, dpi = 300)

# find out TPM average by population:
pag_tpm = tpm_long_filtered %>% filter(population == "PAG")
go_tpm = tpm_long_filtered %>% filter(population == "GO")

m_pag = mean(pag_tpm$TPM) # 30.30832
sd_pag = sqrt(m_pag) # 5.505299

m_go = mean(go_tpm$TPM) # 30.96962
sd_go = sqrt(m_go) # 5.565035

pag_tpm_t1 = pag_tpm %>% filter(timepoint == "T1")
pag_tpm_t2 = pag_tpm %>% filter(timepoint == "T2")

m_pag_t1 = mean(pag_tpm_t1$TPM) # 27.39163
m_pag_t2 = mean(pag_tpm_t2$TPM) # 33.22501

sd_pag_t1 = sqrt(m_pag_t1) # 5.233702
sd_pag_t2 = sqrt(m_pag_t2) # 5.764114

go_tpm_t1 = go_tpm %>% filter(timepoint == "T1")
go_tpm_t2 = go_tpm %>% filter(timepoint == "T2")

m_go_t1 = mean(go_tpm_t1$TPM) # 31.91373
m_go_t2 = mean(go_tpm_t2$TPM) # 30.0228

sd_go_t1 = sqrt(m_go_t1) # 5.649224
sd_go_t2 = sqrt(m_go_t2) # 5.479307

# INDIVIDUAL GENES
# make gene list and export for easier use

summary_tbl <- tpm_long_filtered %>%
  group_by(Gene, population, timepoint) %>%
  summarise(
    mean_TPM = mean(TPM),
    sd_TPM   = sqrt(mean(TPM)),
    .groups = "drop"
  ) %>%
  arrange(Gene, population, timepoint)
write.table(summary_tbl, "TEMP1_GENE_SUMMARY_TPM_values_POPxTIMEPOINT.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# --------------- THIRD CBASS TEMP PHAR SA SY SI ----------------- #

setwd("/UAE_PAG_WGS_PopGen/cbass/2_third_temp/")

# read in TPM data
tpm_df <- read_tsv("results/star_salmon/salmon.merged.gene_tpm.tsv", col_types = cols())[,-2]

# read in annotations file (from Porites harrisoni reference genome) and clean-up gene names
# Porites harrisoni reference genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040938025.1/

annotations <- read_excel("PAG_UKon_Phar.V2.annotations.xlsx") %>%
  select(GeneID, Name) %>%
  mutate(Name = ifelse(is.na(Name), "hypoth_prot", Name))

genes_of_interest <- c("ACVR2B", "BMPR2", "CDK2", "COL4A1_2", "COL2A1_1", "ORC4", "ANKS6_1", "ANO10_1", "DCAF11", "DYNC1LI1", "eva-1", "FBXO31", "HNF4A", "KDELC1_1", "PDCD6_1", "PDCD6_2", "PDRG1", "RAB22A", "RBM18_1", "RSPH10B", "SSNA1", "TTLL9", "TUBB4B", "ZDHHC7_1")
keep_ids <- annotations$GeneID[annotations$Name %in% genes_of_interest]

tpm_subset <- tpm_df %>% filter(gene_id %in% keep_ids)

tpm_subset <- tpm_subset %>%
  left_join(annotations, by = c("gene_id" = "GeneID")) %>%
  filter(Name %in% genes_of_interest)

tpm_long <- tpm_subset %>%
  pivot_longer(
    cols = -c(gene_id, Name),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  select(Sample, Gene = Name, TPM)

tpm_long$Population <- ifelse(grepl("_SI_", tpm_long$Sample), "GO", "PAG")

xvalues = c("ACVR2B", "BMPR2", "CDK2", "COL4A1", "COL2A1", "ORC4", "ANKS6", "ANO10", "DCAF11", "DYNC1LI1", "eva-1", "FBXO31", "HNF4A", "KDELC1", "PDCD6", "PDCD6", "PDRG1", "RAB22A", "RBM18", "RSPH10B", "SSNA1", "TTLL9", "TUBB4B", "ZDHHC7")

ggplot(tpm_long_filtered, aes(x = Gene, y = TPM, fill = Population, color = Population)) +
  geom_boxplot(alpha = 0.7) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.8) +
  scale_fill_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  scale_color_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  scale_x_discrete(labels = xvalues) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(size = 14, color = "black", face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
  ) +
  labs(x = "Gene", y = "TPM", title = "")
ggsave("Genes_of_interest_TPM_boxplot_temp3.pdf", width = 20, height = 10, dpi = 300)

# SPLIT BY T1 & T2 (360 and 1080) #

# take files from above

pattern <- paste0("^(", paste(xvalues, collapse = "|"), ")")

# keep only annotation rows whose Name starts with any of xvalues (e.g. PDCD6_1, PDCD6_2)
annotations_sub <- annot_filt %>%
  filter(str_detect(Name, pattern))

tpm_subset <- tpm_df %>%
  filter(gene_id %in% annotations_sub$GeneID) %>%
  left_join(annotations_sub, by = c("gene_id" = "GeneID"))

tpm_subset <- tpm_subset %>%
  mutate(
    Gene_root = str_extract(
      Name,
      "^[^_]+"
    )
  ) %>%
  # Keep only genes that are in xvalues (safety filter)
  filter(Gene_root %in% xvalues)

tpm_long <- tpm_subset %>%
  pivot_longer(
    cols = -c(gene_id, Name, Gene_root),
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  select(sample, Gene = Gene_root, TPM)
print(tpm_long)

sample_table <- tpm_long %>%
  distinct(sample) %>%
  mutate(
    population = case_when(
      str_detect(sample, "_SY_") ~ "PAG",
      str_detect(sample, "_SA_") ~ "PAG",
      str_detect(sample, "_SI_") ~ "GO",
      TRUE ~ NA_character_
    ),
    timepoint = case_when(
      str_detect(sample, "_360$")  ~ "T1",
      str_detect(sample, "_1080$") ~ "T2",
      TRUE ~ NA_character_
    )
  )
print(sample_table)

tpm_long_meta <- tpm_long %>%
  left_join(sample_table, by = "sample")
print(tpm_long_meta)


tpm_long_filtered <- tpm_long_meta %>%
  filter(TPM != max(TPM, na.rm = TRUE))

gene_levels <- xvalues

df_plot <- tpm_long_filtered %>%
  mutate(
    timepoint  = factor(timepoint,  levels = c("T1", "T2")),
    population = factor(population, levels = c("PAG", "GO")),
    Gene       = factor(Gene,      levels = gene_levels),
    Gene_tp    = paste(Gene, timepoint, sep = "_"),
    Gene_tp    = factor(
      Gene_tp,
      levels = as.vector(rbind(
        paste0(gene_levels, "_T1"),
        paste0(gene_levels, "_T2")
      ))
    )
  )
print(df_plot)

ggplot(df_plot, aes(x = Gene_tp, y = TPM, fill = population)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_jitter(
    aes(color = population),
    position = position_jitterdodge(jitter.width = 0.2),
    size = 1.5,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  scale_color_manual(values = c("PAG" = "#930000", "GO" = "#FF8500")) +
  labs(x = "", y = "TPM") +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(
      size = 14, color = "black", face = "bold",
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title  = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks  = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid   = element_blank()
  ) +
  scale_x_discrete(
    labels = function(x) sub("_(T1|T2)$", "\n\\1", x)
)
ggsave("Genes_of_interest_TPM_boxplot_temp3_BY_TIMEPOINTS.pdf", width = 35, height = 15, dpi = 300)

# find out TPM average by population:
pag_tpm = tpm_long_filtered %>% filter(population == "PAG")
go_tpm = tpm_long_filtered %>% filter(population == "GO")

m_pag = mean(pag_tpm$TPM) # 9.548785
sd_pag = sqrt(m_pag) # 3.090111

m_go = mean(go_tpm$TPM) # 8.511072
sd_go = sqrt(m_go) # 2.917374

pag_tpm_t1 = pag_tpm %>% filter(timepoint == "T1")
pag_tpm_t2 = pag_tpm %>% filter(timepoint == "T2")

m_pag_t1 = mean(pag_tpm_t1$TPM) # 9.464368
m_pag_t2 = mean(pag_tpm_t2$TPM) # 9.629001

sd_pag_t1 = sqrt(m_pag_t1) # 3.076421
sd_pag_t2 = sqrt(m_pag_t2) # 3.103063

go_tpm_t1 = go_tpm %>% filter(timepoint == "T1")
go_tpm_t2 = go_tpm %>% filter(timepoint == "T2")

m_go_t1 = mean(go_tpm_t1$TPM) # 7.625869
m_go_t2 = mean(go_tpm_t2$TPM) # 9.396274

sd_go_t1 = sqrt(m_go_t1) # 2.761498
sd_go_t2 = sqrt(m_go_t2) # 3.065334

# INDIVIDUAL GENES
# make gene list and export for easier use

summary_tbl <- tpm_long_filtered %>%
  group_by(Gene, population, timepoint) %>%
  summarise(
    mean_TPM = mean(TPM),
    sd_TPM   = sqrt(mean(TPM)),
    .groups = "drop"
  ) %>%
  arrange(Gene, population, timepoint)
write.table(summary_tbl, "TEMP3_GENE_SUMMARY_TPM_values_POPxTIMEPOINT.txt", sep = "\t", row.names = FALSE, quote = FALSE)

