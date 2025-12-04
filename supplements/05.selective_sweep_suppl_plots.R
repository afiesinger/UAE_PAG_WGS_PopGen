#!/usr/bin/env Rscript

# Script to plot supplementary figures

# --------------- NUCLEOTIDE DIVERSITY (pi) WITH PIXY --------------- #

# Required libraries
library(dplyr)
library(ggplot2)

# Set working directory
setwd("")

# Load Pixy pi data
pi_data <- read.table("pixy_pi.txt", header = TRUE)

# Filter out windows with NA pi
pi_data <- pi_data[!is.na(pi_data$avg_pi), ]

# Load scaffold lengths
scaffold_lengths <- read.table("sweeps/PAG_UKon_Phar_1.1_contig_len.txt", header = TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]
scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

# Add cumulative position in Mb
pi_data$CHR <- as.character(pi_data$chromosome)
pi_data$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", pi_data$CHR))
pi_data <- pi_data[order(pi_data$scaffold_num, pi_data$window_pos_1), ]
pi_data$CHR <- factor(pi_data$CHR, levels = unique(pi_data$CHR))
pi_data$BPcum_Mb <- (pi_data$window_pos_1 + as.numeric(offsets[as.character(pi_data$CHR)])) / 1e6

# ----------------------------
# Genome-wide plots for each population
# ----------------------------

pi_pag <- pi_data %>% filter(pop == c("SA", "SY"))
pi_go  <- pi_data %>% filter(pop == "SI")

png("PIXY_pi.png", width = 1200, height = 1000)
par(mfrow = c(2, 1))

plot(pi_pag$BPcum_Mb, pi_pag$avg_pi, col = "#930000", type = "l",
     xlab = "Genomic position (Mb)", ylab = expression(pi),
     main = "Genome-wide Nucleotide Diversity (π) - PAG")

plot(pi_go$BPcum_Mb, pi_go$avg_pi, col = "#FF8500", type = "l",
     xlab = "Genomic position (Mb)", ylab = expression(pi),
     main = "Genome-wide Nucleotide Diversity (π) - GO")

dev.off()

pi_pag[which.max(pi_pag$avg_pi),]
#        pop     chromosome window_pos_1 window_pos_2 avg_pi no_sites count_diffs count_comparisons count_missing            CHR scaffold_num BPcum_Mb
# 695184  SA Phar_scaff_276       694001       694500      1        1           1                 1          9023 Phar_scaff_276          276 353.4307

pi_go[which.max(pi_go$avg_pi),]
#        pop     chromosome window_pos_1 window_pos_2    avg_pi no_sites count_diffs count_comparisons count_missing            CHR scaffold_num BPcum_Mb
# 691719  SI Phar_scaff_276       694001       694500 0.6666667        1           4                 6          5034 Phar_scaff_276          276 353.4307

# ----------------------------
# Scaffold-specific plots
# ----------------------------

pi_pag$Population <- rep("PAG", nrow(pi_pag))
pi_go$Population <- rep("GO", nrow(pi_go))

plot_scaffold <- function(scaff_id, pop1, pop2, col1 = "#930000", col2 = "#FF8500") {

  df1 <- pop1 %>% filter(CHR == scaff_id) %>% mutate(Population = "PAG")
  df2 <- pop2 %>% filter(CHR == scaff_id) %>% mutate(Population = "GO")
  df <- bind_rows(df1, df2)

  y_min <- 0
  y_max <- max(df$avg_pi, na.rm = TRUE) * 1.05

  cols <- c("PAG" = col1, "GO" = col2)

  p <- ggplot(df, aes(x = BPcum_Mb, y = avg_pi, color = Population)) +
    geom_point(alpha = 0.8) +
    geom_line(linewidth = 0.7) +
    labs(x = "Genomic position (Mb)",
         y = expression(pi),
         title = paste("Nucleotide Diversity (π) -", scaff_id)) +
    ylim(y_min, y_max) +
    scale_color_manual(values = cols) +
    theme_bw(base_size = 20) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 14, color = "black", face = "bold", hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.border = element_rect(color = "black"),
      panel.grid = element_blank())
  ggsave(paste0(scaff_id, "_pi_PIXY.png"), plot = p, width = 25, height = 12, dpi = 300)
}

plot_scaffold("Phar_scaff_29", pi_pag, pi_go)
plot_scaffold("Phar_scaff_62", pi_pag, pi_go)

# -------------------- iHS ---------------------- #

setwd("T")  

ihs_pag <- read.table("ALL_PAG.ihs.out.norm.header", header = TRUE)
ihs_go <- read.table("ALL_GO.ihs.out.norm.header", header = TRUE)

# # PLOT HISTOGRAM OF EACH POP # # 

ihs_pag$sihs <- as.numeric(ihs_pag$sihs)
ihs_go$sihs <- as.numeric(ihs_go$sihs)

# outlier threshold top 1% outliers
ihs_outl_pag <- quantile(ihs_pag$sihs, 0.99, na.rm = TRUE) # 2.670701 
ihs_outl_go <- quantile(ihs_go$sihs, 0.99, na.rm = TRUE) # 1.763913 

ihs_pag[which.max(ihs_pag$sihs), ]

#                           id    pos       daf    ihh1    ihh0    uihs    sihs bin
# 173116 Phar_scaff_191_875493 875493 0.0714286 16935.4 50.1888 2.52819 6.63853   1

ihs_pag[which.min(ihs_pag$sihs), ]

#                           id    pos       daf    ihh1    ihh0     uihs     sihs bin
# 197690 Phar_scaff_204_490473 490473 0.0873016 66.8455 3936.67 -1.77006 -6.38241   1

# create the histogram
ggplot(ihs_pag, aes(x = sihs)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = c(ihs_outl_pag, -ihs_outl_pag), color = "red", linetype = "dotted") +
  labs(x = "iHS", y = "Frequency", title = "Histogram of iHS PAG") +
  theme_bw()
ggsave("iHS_histo_PAG_1pc.png", dpi = 300)

ggplot(ihs_go, aes(x = sihs)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = c(ihs_outl_go, -ihs_outl_go), color = "red", linetype = "dotted") +
  labs(x = "iHS", y = "Frequency", title = "Histogram of iHS Go") +
  theme_bw()
ggsave("iHS_histo_GO_1pc.png", dpi = 300)

scaffold_lengths <- read.table("PAG_UKon_Phar_Genome_scaffold_lengths.txt.header", header = TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

prepare_ihs_data <- function(df, offsets_vec, scaffold_lengths_df) {
  df$CHR <- sub("^(Phar_scaff_\\d+)_\\d+$", "\\1", df$id)
  df$POS <- as.numeric(df$pos)
  df$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", df$CHR))
  df <- df[order(df$scaffold_num, df$POS), ]
  df$CHR <- factor(df$CHR, levels = scaffold_lengths_df$CHR)
  cum_offsets <- as.numeric(offsets_vec[as.character(df$CHR)])
  df$BPcum_Mb <- (df$POS + cum_offsets) / 1000000
  return(df)
}

ihs_pag <- prepare_ihs_data(ihs_pag, offsets, scaffold_lengths)
ihs_go <- prepare_ihs_data(ihs_go, offsets, scaffold_lengths)

assign_colors <- function(df) {
  chroms <- unique(df$CHR)
  nchrom <- length(chroms)
  colors <- rep(c("black", "dodgerblue3"), length.out = nchrom)
  names(colors) <- chroms
  df$color <- colors[as.character(df$CHR)]
  return(df)
}

ihs_pag <- assign_colors(ihs_pag)
ihs_go <- assign_colors(ihs_go)

### PLOT iHS VALUES ###
png("iHS_1pc.png", width = 1200, height = 800)
par(mfrow = c(2,1))

plot(ihs_pag$BPcum_Mb, ihs_pag$sihs, pch = 20, col = ihs_pag$color, cex = 0.6,
     xlab = "Genomic position (Mb)", ylab = "Standardized iHS",
     main = "Genome-wide iHS scan (PAG)")
abline(h = c(ihs_outl_pag, -ihs_outl_pag), col = "red", lty = 2, lwd = 1.5)

plot(ihs_go$BPcum_Mb, ihs_go$sihs, pch = 20, col = ihs_go$color, cex = 0.6,
     xlab = "Genomic position (Mb)", ylab = "Standardized iHS",
     main = "Genome-wide iHS scan (GO)")
abline(h = c(ihs_outl_go, -ihs_outl_go), col = "red", lty = 2, lwd = 1.5)

dev.off()

# # plot only the outliers # # 

# Subset data for PAG and GO iHS
df_PAG_pos <- ihs_pag %>% filter(sihs >= ihs_outl_pag)
df_PAG_neg <- ihs_pag %>% filter(sihs <= -ihs_outl_pag)
df_GO_pos <- ihs_go %>% filter(sihs >= ihs_outl_go)
df_GO_neg <- ihs_go %>% filter(sihs <= -ihs_outl_go)

# For negative values, convert iHS to absolute values for plotting
df_PAG_neg$abs_iHS <- abs(df_PAG_neg$sihs)
df_GO_neg$abs_iHS <- abs(df_GO_neg$sihs)

# Set up a 2x2 plotting layout
png("iHS_OUTLIERS_PAGvGO.png", width = 2000, height = 1000)
par(mfrow = c(2, 2), mar=c(4,4,2,1))

# PAG positive iHS
plot(df_PAG_pos$BPcum_Mb, df_PAG_pos$sihs,
     pch = 20, cex = 0.6, 
     col = df_PAG_pos$color,
     xlab="Genomic position (Mb)",
     ylab="Standardized iHS",
     main="PAG iHS (Positive outliers)"
)
#abline(h=pos_threshold, col="red", lty=2)

# PAG negative iHS (absolute values)
plot(df_PAG_neg$BPcum_Mb, df_PAG_neg$abs_iHS,
     pch = 20, cex = 0.6, 
     col = df_PAG_neg$color,
     xlab="Genomic position (Mb)",
     ylab="Absolute Standardized iHS",
     main="PAG iHS (Negative outliers)"
)
#abline(h=pos_threshold, col="red", lty=2)

# GO positive iHS
plot(df_GO_pos$BPcum_Mb, df_GO_pos$sihs,
     pch = 20, cex = 0.6, 
     col = df_GO_pos$color,
     xlab="Genomic position (Mb)",
     ylab="Standardized iHS",
     main="GO iHS (Positive outliers)"
)
#abline(h=pos_threshold, col="red", lty=2)

# GO negative iHS (absolute values)
plot(df_GO_neg$BPcum_Mb, df_GO_neg$abs_iHS,
     pch = 20, cex = 0.6, 
     col = df_GO_neg$color,
     xlab="Genomic position (Mb)",
     ylab="Absolute Standardized iHS",
     main="GO iHS (Negative outliers)"
)
#abline(h=pos_threshold, col="red", lty=2)

dev.off()

# ------------ nSL --------------- #

setwd("/")

nsl_pag <- read.table("ALL_PAG.nsl.out.norm.header", header = TRUE)
nsl_go <- read.table("ALL_GO.nsl.out.norm.header", header = TRUE)

# # PLOT HISTOGRAM # # 

nsl_pag$normnsl <- as.numeric(nsl_pag$normnsl)
nsl_go$normnsl <- as.numeric(nsl_go$normnsl)

# outlier threshold (99% quantile) -- top 1% outliers
nsl_outl_pag <- quantile(nsl_pag$normnsl, 0.99, na.rm = TRUE)  
nsl_outl_go <- quantile(nsl_go$normnsl, 0.99, na.rm = TRUE) 

# Create the histogram
ggplot(nsl_pag, aes(x = normnsl)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = c(nsl_outl_pag, -nsl_outl_pag), color = "red", linetype = "dotted") +
  labs(x = "nSL", y = "Frequency", title = "Histogram of nSL PAG") +
  theme_bw()
ggsave("nSL_histo_PAG_1pc.png", dpi = 300)

ggplot(nsl_go, aes(x = normnsl)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = c(nsl_outl_go, -nsl_outl_go), color = "red", linetype = "dotted") +
  labs(x = "nSL", y = "Frequency", title = "Histogram of nSL Go") +
  theme_bw()
ggsave("nSL_histo_GO_1pc.png", dpi = 300)

scaffold_lengths <- read.table("sweeps/PAG_UKon_Phar_1.1_contig_len.txt", header = TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

prepare_data <- function(df, offsets_vec, scaffold_lengths_df) {
  df$CHR <- sub("^(Phar_scaff_\\d+)_\\d+$", "\\1", df$id)
  df$POS <- as.numeric(df$pos)
  df$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", df$CHR))
  df <- df[order(df$scaffold_num, df$POS), ]
  df$CHR <- factor(df$CHR, levels = scaffold_lengths_df$CHR)
  cum_offsets <- as.numeric(offsets_vec[as.character(df$CHR)])
  df$BPcum_Mb <- (df$POS + cum_offsets) / 1000000
  return(df)
}

nsl_pag <- prepare_data(nsl_pag, offsets, scaffold_lengths)
nsl_go <- prepare_data(nsl_go, offsets, scaffold_lengths)

assign_colors <- function(df) {
  chroms <- unique(df$CHR)
  nchrom <- length(chroms)
  colors <- rep(c("black", "dodgerblue3"), length.out = nchrom)
  names(colors) <- chroms
  df$color <- colors[as.character(df$CHR)]
  return(df)
}

nsl_pag <- assign_colors(nsl_pag)
nsl_go <- assign_colors(nsl_go)

nsl_pag[which.max(nsl_pag$normnsl), ]

#                          id    pos      daf     sl1     sl0   unsl normnsl bin           CHR    POS scaffold_num BPcum_Mb color
# 634999 Phar_scaff_85_316935 316935 0.253968 106.224 3.76985 1.4499 4.62162   1 Phar_scaff_85 316935           85  171.136 black

nsl_pag[which.min(nsl_pag$normnsl), ]

#                          id    pos      daf     sl1     sl0     unsl  normnsl bin           CHR    POS scaffold_num BPcum_Mb       color
# 541681 Phar_scaff_62_595849 595849 0.246032 4.42796 76.3284 -1.23648 -5.14638   1 Phar_scaff_62 595849           62 137.6717 dodgerblue3

### PLOT nSL VALUES ###
png("nSL_1pc_PAG_GO.png", width = 1200, height = 800)
par(mfrow = c(2,1))

plot(nsl_pag$BPcum_Mb, nsl_pag$normnsl, pch = 20, col = nsl_pag$color, cex = 0.6,
     xlab = "Genomic position (Mb)", ylab = "Normalized nSL",
     main = "Genome-wide nSL scan (PAG)")
abline(h = c(nsl_outl_pag, -nsl_outl_pag), col = "red", lty = 2, lwd = 1.5)

plot(nsl_go$BPcum_Mb, nsl_go$normnsl, pch = 20, col = nsl_go$color, cex = 0.6,
     xlab = "Genomic position (Mb)", ylab = "Normalized nSL",
     main = "Genome-wide nSL scan (GO)")
abline(h = c(nsl_outl_go, -nsl_outl_go), col = "red", lty = 2, lwd = 1.5)

dev.off()

# # plot only the outliers # # 

# Subset data for PAG and GO iHS
df_PAG_pos <- nsl_pag %>% filter(normnsl >= nsl_outl_pag)
df_PAG_neg <- nsl_pag %>% filter(normnsl <= -nsl_outl_pag)
df_GO_pos <- nsl_go %>% filter(normnsl >= nsl_outl_go)
df_GO_neg <- nsl_go %>% filter(normnsl <= -nsl_outl_go)

# For negative values, convert iHS to absolute values for plotting
df_PAG_neg$abs_nSL <- abs(df_PAG_neg$normnsl)
df_GO_neg$abs_nSL <- abs(df_GO_neg$normnsl)

# Set up a 2x2 plotting layout
png("nSL_OUTLIERS_PAG_GO.png", width = 2000, height = 1000)
par(mfrow = c(2, 2), mar=c(4,4,2,1))

plot(df_PAG_pos$BPcum_Mb, df_PAG_pos$normnsl,
     pch = 20, cex = 0.6, 
     col = df_PAG_pos$color,
     xlab="Genomic position (Mb)",
     ylab="Normalized nSL",
     main="PAG nSL (Positive outliers)"
)

plot(df_PAG_neg$BPcum_Mb, df_PAG_neg$abs_nSL,
     pch = 20, cex = 0.6, 
     col = df_PAG_neg$color,
     xlab="Genomic position (Mb)",
     ylab="Absolute Normalized nSL",
     main="PAG nSL (Negative outliers)"
)

plot(df_GO_pos$BPcum_Mb, df_GO_pos$normnsl,
     pch = 20, cex = 0.6, 
     col = df_GO_pos$color,
     xlab="Genomic position (Mb)",
     ylab="Normalized nSL",
     main="GO nSL (Positive outliers)"
)

plot(df_GO_neg$BPcum_Mb, df_GO_neg$abs_nSL,
     pch = 20, cex = 0.6, 
     col = df_GO_neg$color,
     xlab="Genomic position (Mb)",
     ylab="Absolute Normalized nSL",
     main="GO nSL (Negative outliers)"
)

dev.off()

# ------------ iHH12 ------------- #

setwd("/")

ihh12_pag <- read.table("ALL_PAG.ihh12.out.norm", header = TRUE)
ihh12_go  <- read.table("ALL_GO.ihh12.out.norm", header = TRUE)

# # PLOT HISTOGRAM # # 

ihh12_pag$normihh12 <- as.numeric(ihh12_pag$normihh12)
ihh12_go$normihh12  <- as.numeric(ihh12_go$normihh12)

# outlier threshold (99% quantile) -- top 1% outliers
ihh12_outl_pag <- quantile(ihh12_pag$normihh12, 0.99, na.rm = TRUE) 
ihh12_outl_go  <- quantile(ihh12_go$normihh12, 0.99, na.rm = TRUE) 

ihh12_pag[which.max(ihh12_pag$normihh12), ]

#                          id    pos        p1   ihh12 normihh12 crit
# 541718 Phar_scaff_62_628739 628739 0.0555556 62462.5   8.47318    1

ihh12_pag[which.min(ihh12_pag$normihh12), ]

#                          id   pos      p1 ihh12 normihh12 crit
# 33117 Phar_scaff_1139_13529 13529 0.31746     0  -1.10857    0

ggplot(ihh12_pag, aes(x = normihh12)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = ihh12_outl_pag, color = "red", linetype = "dotted") +
  labs(x = "iHH12", y = "Frequency", title = "Histogram of iHH12 PAG") +
  theme_bw()
ggsave("iHH12_histo_PAG_1pc.png", dpi = 300)

ggplot(ihh12_go, aes(x = normihh12)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = ihh12_outl_go, color = "red", linetype = "dotted") +
  labs(x = "iHH12", y = "Frequency", title = "Histogram of iHH12 GO") +
  theme_bw()
ggsave("iHH12_histo_GO_1pc.png", dpi = 300)

scaffold_lengths <- read.table("PAG_UKon_Phar_Genome_scaffold_lengths.txt.header", header = TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

prepare_data <- function(df, offsets_vec, scaffold_lengths_df) {
  df$CHR <- sub("^(Phar_scaff_\\d+)_\\d+$", "\\1", df$id)
  df$POS <- as.numeric(df$pos)
  df$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", df$CHR))
  df <- df[order(df$scaffold_num, df$POS), ]
  df$CHR <- factor(df$CHR, levels = scaffold_lengths_df$CHR)
  cum_offsets <- as.numeric(offsets_vec[as.character(df$CHR)])
  df$BPcum_Mb <- (df$POS + cum_offsets) / 1000000
  return(df)
}

ihh12_pag <- prepare_data(ihh12_pag, offsets, scaffold_lengths)
ihh12_go  <- prepare_data(ihh12_go, offsets, scaffold_lengths)

assign_colors <- function(df) {
  chroms <- unique(df$CHR)
  nchrom <- length(chroms)
  colors <- rep(c("black", "dodgerblue3"), length.out = nchrom)
  names(colors) <- chroms
  df$color <- colors[as.character(df$CHR)]
  return(df)
}

ihh12_pag <- assign_colors(ihh12_pag)
ihh12_go  <- assign_colors(ihh12_go)

### PLOT iHH12 VALUES ###
png("iHH12_1p_PAG_GO.png", width = 1200, height = 800)
par(mfrow = c(2,1))

plot(ihh12_pag$BPcum_Mb, ihh12_pag$normihh12, pch = 20, col = ihh12_pag$color, cex = 0.6,
     xlab = "Genomic position (Mb)", ylab = "Normalized iHH12",
     main = "Genome-wide iHH12 scan (PAG)")
abline(h = ihh12_outl_pag, col = "red", lty = 2, lwd = 1.5)

plot(ihh12_go$BPcum_Mb, ihh12_go$normihh12, pch = 20, col = ihh12_go$color, cex = 0.6,
     xlab = "Genomic position (Mb)", ylab = "Normalized iHH12",
     main = "Genome-wide iHH12 scan (GO)")
abline(h = ihh12_outl_go, col = "red", lty = 2, lwd = 1.5)

dev.off()

# # plot only the outliers # # 

# Subset data for PAG and GO iHH12
df_PAG_pos <- ihh12_pag %>% filter(normihh12 >= ihh12_outl_pag)
df_GO_pos  <- ihh12_go %>% filter(normihh12 >= ihh12_outl_go)

png("iHH12_1p_OUTLIERS_PAG_GO.png", width = 2000, height = 1000)
par(mfrow = c(2, 1))

plot(df_PAG_pos$BPcum_Mb, df_PAG_pos$normihh12,
     pch = 20, cex = 0.6, 
     col = df_PAG_pos$color,
     xlab="Genomic position (Mb)",
     ylab="Normalized iHH12",
     main="PAG iHH12 outliers"
)

plot(df_GO_pos$BPcum_Mb, df_GO_pos$normihh12,
     pch = 20, cex = 0.6, 
     col = df_GO_pos$color,
     xlab="Genomic position (Mb)",
     ylab="Normalized iHH12",
     main="GO iHH12 outliers"
)

dev.off()

# # MANUSCRIPT SUPPLEMENT PLOTS # #

setwd("")

library(ggplot2)
library(patchwork) 

plot_scaffold_all_measures <- function(scaffold, ihs_df, nsl_df, ihh12_df,
                                       outlier_ihs, outlier_nsl, outlier_ihh12,
                                       out_prefix = NULL) {
  df_ihs   <- ihs_df[ihs_df$CHR == scaffold, ]
  df_nsl   <- nsl_df[nsl_df$CHR == scaffold, ]
  df_ihh12 <- ihh12_df[ihh12_df$CHR == scaffold, ]

  out_file <- paste0(out_prefix, "_", scaffold, ".png")

  p1 <- ggplot(df_ihs, aes(x = BPcum_Mb, y = sihs)) +
    geom_point(size = 0.9, color = "black") +
    geom_hline(yintercept = c(outlier_ihs, -outlier_ihs), color = "#626269ff", linetype = "dashed", linewidth = 1.5) +
    labs(x = "",
         y = "Standardized iHS",
         title = paste("Integrated haplotype score (iHS)")) +
    theme_bw(base_size = 20) +
    theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
    )

  p2 <- ggplot(df_nsl, aes(x = BPcum_Mb, y = normnsl)) +
    geom_point(size = 0.9, color = "black") +
    geom_hline(yintercept = c(outlier_nsl, -outlier_nsl), color = "#626269ff", linetype = "dashed", linewidth = 1.5) +
    labs(x = "",
         y = "Normalized nSL",
         title = paste("Number of segregating sites by length (nSL)")) +
    theme_bw(base_size = 20) +
    theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
    )

  p3 <- ggplot(df_ihh12, aes(x = BPcum_Mb, y = normihh12)) +
    geom_point(size = 0.9, color = "black") +
    geom_hline(yintercept = outlier_ihh12, color = "#626269ff", linetype = "dashed", linewidth = 1.5) +
    labs(x = "Genomic Position (Mb)",
         y = "Normalized iHH12",
         title = paste("Integrated haplotype homozygosity pooled (iHH12)")) +
    theme_bw(base_size = 20) +
    theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
    )
  
  combined_plot <- p1 / p2 / p3
  ggsave(out_file, combined_plot, width = 20, height = 26, dpi = 300)
  message(paste("Saved:", out_file))
}

my_scaffolds <- c("Phar_scaff_29", "Phar_scaff_62")

# PAG
for (scaf in my_scaffolds) {
  plot_scaffold_all_measures(scaf,
                             ihs_df   = ihs_pag,
                             nsl_df   = nsl_pag,
                             ihh12_df = ihh12_pag,
                             outlier_ihs   = ihs_outl_pag,
                             outlier_nsl   = nsl_outl_pag,  
                             outlier_ihh12 = ihh12_outl_pag,
                             out_prefix = "H-STATS_PAG")
}

# ----------------- ALLELE FREQUENCY --------------- #

setwd("")
library(ggplot2)

pag <- read.table("PAG_mod.frq", header=TRUE, stringsAsFactors=FALSE)
go <- read.table("GO_mod.frq", header=TRUE, stringsAsFactors=FALSE)

scaffold_lengths <- read.table("sweeps/PAG_UKon_Phar_1.1_contig_len.txt", header = TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

# add cumulative positions in Mb
pag$BPcum_Mb   <- (pag$pos + as.numeric(offsets[as.character(pag$scaff)])) / 1e6
go$BPcum_Mb   <- (go$pos + as.numeric(offsets[as.character(go$scaff)])) / 1e6

contig_of_interest <- "Phar_scaff_29" 
contig_of_interest <- "Phar_scaff_62" 

pag_contig <- pag[pag$scaff == contig_of_interest, ]
go_contig <- go[go$scaff == contig_of_interest, ]

extract_first_af <- function(allele_freq_cols) {
  sapply(allele_freq_cols, function(x) {
    parts <- strsplit(x, ":")[[1]]
    as.numeric(parts[2])
  })
}

pag_contig$AF <- extract_first_af(as.character(pag_contig[, 6]))
go_contig$AF <- extract_first_af(as.character(go_contig[, 6]))

library(ggplot2)
library(ggpubr)

# Plot for PAG population
plot_pag <- ggplot(pag_contig, aes(x = BPcum_Mb, y = AF)) +
  geom_point(color = "#930000", shape = 1, alpha = 0.9) +
  ylim(0, 1) +
  labs(x = "", y = "Allele Frequency") +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)

# Plot for GO population
plot_go <- ggplot(go_contig, aes(x = BPcum_Mb, y = AF)) +
  geom_point(color = "#FF8500", shape = 1, alpha = 0.9) +
  ylim(0, 1) +
  labs(x = "Genomic position (Mb)", y = "Allele Frequency") +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)

pg = ggarrange(plot_pag, plot_go, ncol = 1)
ggsave("Allele_freq_PAG_GO_29.png", pg, width = 20, height = 15, dpi = 300)
ggsave("Allele_freq_PAG_GO_62.png", pg, width = 20, height = 15, dpi = 300)

# plot contig 29 only in the Fst peak range

p_pag <- ggplot(pag_contig, aes(x = BPcum_Mb, y = AF)) +
  geom_point(color = "#930000", shape = 1, alpha = 0.8) +
  ylim(0, 1) +
  xlim(80.4, 80.6) +
  labs(x = "", y = "Allele Frequency") +
  theme_bw(base_size = 20) +
    theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)

p_go <- ggplot(go_contig, aes(x = BPcum_Mb, y = AF)) +
  geom_point(color = "#FF8500", shape = 1, alpha = 0.8) +
  ylim(0, 1) +
  xlim(80.4, 80.6) +
  labs(x = "Genomic position (Mb)", y = "Allele Frequency") +
  theme_bw(base_size = 20) +
    theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)

p = ggarrange(p_pag, p_go, ncol = 1)
ggsave("Allele_freq_PAG_GO_29_zoom.png", p, width = 20, height = 15, dpi = 300)

# plot contig 62 only in the xpehh peak range

g_pag <- ggplot(pag_contig, aes(x = BPcum_Mb, y = AF)) +
  geom_point(color = "#930000", shape = 1, alpha = 0.8) +
  ylim(0, 1) +
  xlim(137.25, 138.0) +
  labs(x = "", y = "Allele Frequency") +
  theme_bw(base_size = 20) +
    theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)

g_go <- ggplot(go_contig, aes(x = BPcum_Mb, y = AF)) +
  geom_point(color = "#FF8500", shape = 1, alpha = 0.8) +
  ylim(0, 1) +
  xlim(137.25, 138.0) +
  labs(x = "Genomic position (Mb)", y = "Allele Frequency") +
  theme_bw(base_size = 20) +
    theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    title = element_text(size = 22, face = "bold", color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)

g = ggarrange(g_pag, g_go, ncol = 1)
ggsave("Allele_freq_PAG_GO_62_zoom.png", g, width = 20, height = 15, dpi = 300)

pag_go = ggarrange(p, g, ncol = 1, nrow = 2)
ggsave("Allele_freq_PAG_GO_29_62_zoom.svg", pag_go, width = 20, height = 24, dpi = 300)

