#!/usr/bin/env Rscript

# --- LIBRARIES --- #

library(ggplot2)
library(dplyr)
library(readr)
library(data.table)

# ---------------- Fst --------------- #

setwd("")

fst <- read.table("PHAR_Fst_ALL.weir.fst", header = TRUE)
fst <- fst[!is.na(fst$WEIR_AND_COCKERHAM_FST), ]

fst$CHR <- as.character(fst$CHROM)
fst$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", fst$CHR))

fst <- fst[order(fst$scaffold_num, fst$POS), ]
fst$CHR <- factor(fst$CHR, levels = unique(fst$CHR))

scaffold_lengths <- read.table("sweeps/PAG_UKon_Phar_1.1_contig_len.txt", header = TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

fst$BPcum_Mb <- (fst$POS + as.numeric(offsets[as.character(fst$CHR)])) / 1000000

outlier_threshold <- quantile(fst$WEIR_AND_COCKERHAM_FST, 0.999, na.rm = TRUE) 

outliers_fst <- fst[fst$WEIR_AND_COCKERHAM_FST > outlier_threshold, ]
write.table(outliers_fst, "PHAR_Fst_ALL.OUTLIERS_0.1p.txt", sep = "\t", row.names = FALSE, quote = FALSE)

max(fst$WEIR_AND_COCKERHAM_FST) # 0.930765

max_fst_scaff = fst %>% filter(fst$WEIR_AND_COCKERHAM_FST > 0.93)
print(max_fst_scaff)

#           CHROM    POS WEIR_AND_COCKERHAM_FST           CHR scaffold_num BPcum_Mb
# 1 Phar_scaff_29 526909               0.930765 Phar_scaff_29           29 80.46393

outliers_fst$BPcumBP <- outliers_fst$BPcum_Mb * 1000000
outliers_fst$snpid <- paste0(outliers_fst$CHR, "_", outliers_fst$BPcumBP)

write.table(outliers_fst, "PHAR_Fst_ALL.OUTLIERS_0.1p_with_SNPIDs.txt", sep="\t", quote=FALSE, row.names=FALSE)

# --------------------- XP-EHH ---------------------- #

library(data.table)
library(ggplot2)

setwd("")  

xpehh <- read.table("XPEHH_ALL_SCAFF.xpehh.out.norm", header = TRUE)
xpehh$CHR <- sub("^(Phar_scaff_\\d+)_\\d+$", "\\1", xpehh$id)
xpehh$POS <- as.numeric(xpehh$pos)
xpehh$normxpehh <- as.numeric(xpehh$normxpehh)

scaffold_lengths <- fread("sweeps/PAG_UKon_Phar_1.1_contig_len.txt", header = TRUE)
scaffold_lengths$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", scaffold_lengths$CHR))
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$scaffold_num), ]

# compute cumulative lengths (starting base of each scaffold)
scaffold_lengths$cumlen <- cumsum(scaffold_lengths$LENGTH)
offsets <- c(0, head(scaffold_lengths$cumlen, -1))
names(offsets) <- scaffold_lengths$CHR

# ensure CHR factor levels match scaffold order
xpehh$scaffold_num <- as.numeric(gsub("Phar_scaff_", "", xpehh$CHR))
xpehh <- xpehh[order(xpehh$scaffold_num, xpehh$POS), ]
xpehh$CHR <- factor(xpehh$CHR, levels = unique(scaffold_lengths$CHR))

# compute cumulative position in Mb
cum_offsets <- as.numeric(offsets[as.character(xpehh$CHR)])
xpehh$BPcum_Mb <- (xpehh$POS + cum_offsets) / 1000000

write.table(xpehh, "XPEHH_ALL_SCAFF.xpehh.out.norm_CUM_POSITIONS", sep = "\t", row.names = FALSE, quote = FALSE)

xpehh_outliers <- fread("top.1pc.normxpehh.PAG.outliers.header.chr", header = TRUE)
xpehh_outliers$normxpehh <- as.numeric(xpehh_outliers$normxpehh)
min(xpehh_outliers$normxpehh) # 3.09011
max(xpehh_outliers$normxpehh) # 10.4701

outl_abs = min(xpehh_outliers$normxpehh)

# assign alternating colors
chroms <- unique(xpehh$CHR)
nchrom <- length(chroms)
colors <- rep(c("black", "dodgerblue3"), length.out = nchrom)
names(colors) <- chroms

xpehh$color <- colors[as.character(xpehh$CHR)]

# create genome-wide XP-EHH plot
png("XPEHH_norm_1p.png", width = 1500, height = 700)
plot(xpehh$BPcum_Mb, xpehh$normxpehh, col = xpehh$color, pch = 20, cex = 0.6, xlab = "Genomic position (Mb)", ylab = "Normalized XP-EHH", main = "Genome-wide XP-EHH scan")
abline(h = outl_abs, col = "red", lty = 2, lwd = 1.5)
dev.off()

# create the histogram
ggplot(xpehh, aes(x = normxpehh)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  geom_vline(xintercept = outl_abs, color = "red", linetype = "dotted") +
  labs(
    title = "Histogram of Normalized XP-EHH Values",
    x = "Normalized XP-EHH (normxpehh)",
    y = "Frequency"
  ) +
  theme_bw()
ggsave("XPEHH_histo_1p.png", dpi = 300)

# find highest scaffold:
xpehh[which.max(xpehh$normxpehh), ]

#                           id    pos   gpos        p1    ihh1   p2    ihh2   xpehh normxpehh crit           CHR    POS scaffold_num BPcum_Mb       color
# 2277084 Phar_scaff_62_597833 597833 597833 0.0714286 23905.6 0.25 516.134 1.66574   10.4701    1 Phar_scaff_62 597833           62 137.6737 dodgerblue3

outl_scaff62 <- xpehh_outliers %>% filter(xpehh_outliers$scaff == "Phar_scaff_62")
write.table(outl_scaff62, "SCAFF62_OUTLIERS_XPEHH_withSNPID.forSNPEff.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# also for scaff 29 (identified through Fst)
outl_scaff29 <- xpehh_outliers %>% filter(xpehh_outliers$scaff == "Phar_scaff_29")
write.table(outl_scaff62, "SCAFF29_OUTLIERS_XPEHH_withSNPID.forSNPEff.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------- FIND OVERLAPPING SCAFFOLDS BETWEEN FST & XP-EHH  --------------------------- #

xpehh.out <- read.table("/path/to/sweeps/EHH/XP-EHH_OUT/top.1pc.normxpehh.PAG.outliers.header.chr.SNPids", header = TRUE, sep="", stringsAsFactors=FALSE)
fst.out <- read.table("path/to/sweeps/Fst/PHAR_Fst_ALL.OUTLIERS_0.1p_with_SNPIDs.txt", header=TRUE, sep="", stringsAsFactors=FALSE)

# find overlapping scaffolds
overlap_fst_xpehh <- intersect(as.character(unique(fst.out$CHROM)), as.character(unique(xpehh.out$scaff)))
print(overlap_fst_xpehh)

# find overlapping SNPs
overlap_fst_xpehh_SNP <- intersect(unique(fst.out$snpid), unique(xpehh.out$snpid))
print(overlap_fst_xpehh_SNP)
# total of 394 SNPs

# how many on either contig 29 and 62 

# CONTIG 29
scaff29 = fst.out %>% filter(fst.out$CHROM == "Phar_scaff_29")
scaff29.x = xpehh.out %>% filter(xpehh.out$chr == 29)
overlap = intersect(unique(scaff29$snpid), unique(scaff29.x$snpid))
# 50 SNPs

# CONTIG 62
scaff62 = fst.out %>% filter(fst.out$CHROM == "Phar_scaff_62")
scaff62.x = xpehh.out %>% filter(xpehh.out$chr == 62)
overlap = intersect(unique(scaff62$snpid), unique(scaff62.x$snpid))
# 49 SNPs

# -------------------- UpSet PLOT ----------------- #

library(UpSetR)

setwd("")

# 0.1 % Fst outliers & 1% XPEHH outliers #
fst.out = read.table("/home/fiesingera/proj/PAG_Phar_PopGen/NEW_FILTERS_VCF/SWEEPS/Fst/PHAR_Fst_ALL.OUTLIERS_0.1p_with_SNPIDs.txt", header=TRUE, sep="", stringsAsFactors=FALSE)
head(fst.out)

xpehh.out = read.table("/home/fiesingera/proj/PAG_Phar_PopGen/NEW_FILTERS_VCF/SWEEPS/EHH/XPEHH_OUT/top.1pc.normxpehh.PAG.outliers.header.chr.SNPids", header = TRUE, sep="", stringsAsFactors=FALSE)
head(xpehh.out)

# SNPs overlap
overlap = intersect(unique(fst.out$snpid), unique(xpehh.out$snpid))
print(overlap) 

outlist = list(fst = fst.out$snpid, xpehh = xpehh.out$snpid)

png("1_UpSet_plot_xpehh_fst_SNPs.png", width = 30, height = 30, unit = "cm", res = 600)
upset(fromList(outlist), order.by = "freq", mainbar.y.label = "Outlier Intersections", sets.x.label = "Outlier SNPs", text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))
dev.off()

# contigs overlap
overlap = intersect(as.character(unique(fst.out$CHROM)), as.character(unique(xpehh.out$scaff)))
print(overlap) 

outlist = list(fst = fst.out$CHROM, xpehh = xpehh.out$scaff)

png("2_UpSet_plot_xpehh_fst_CONTIGS.png", width = 30, height = 30, unit = "cm", res = 600)
upset(fromList(outlist), order.by = "freq", mainbar.y.label = "Outlier Intersections", sets.x.label = "Outlier SNPs", text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))
dev.off()

# ----------- PLOT XPEHH & FST TOGETHER FOR CANDIDATE REGIONS OF SELECTION ----------- #

# OVERVIEW PLOT OF XPEHH WITH OUTLIER SCAFFOLDS HIGHLIGHTED #
setwd("")

highlighted_scaffolds = c("Phar_scaff_29", "Phar_scaff_62")

# ensure CHR column exists as character
xpehh$CHR <- as.character(xpehh$CHR)

# assign alternating grey shades to scaffolds
unique_chr <- unique(xpehh$CHR)
grey_colors <- rep(c("#97979dff", "#323233ff"), length.out = length(unique_chr))
names(grey_colors) <- unique_chr

# override with purple for highlighted scaffolds
custom_colors <- grey_colors
custom_colors[highlighted_scaffolds] <- "#192bc2"

# assign point color to each row
xpehh$color <- custom_colors[xpehh$CHR]

# plot
xpehh_h <- ggplot(xpehh, aes(x = BPcum_Mb, y = normxpehh, color = color)) +
  geom_point(size = 0.7) +
  labs(x = "Genomic position (Mb)",
    y = "Normalized XP-EHH") +
  scale_color_identity() +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
  ) 
ggsave("XPEHH_norm_for_MS_highlighted.png", xpehh_h, width = 15, height = 7, dpi = 300)

# without the x axis label for the combined plot below
xpehh_h <- ggplot(xpehh, aes(x = BPcum_Mb, y = normxpehh, color = color)) +
  geom_point(size = 0.7) +
  labs(x = "",
    y = "Normalized XP-EHH") +
  scale_color_identity() +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
) 

# OVERVIEW PLOT OF FST WITH OUTLIER SCAFFOLDS HIGHLIGHTED #
setwd("")

highlighted_scaffolds = c("Phar_scaff_29", "Phar_scaff_62")

# ensure CHR column exists as character
fst$CHR <- as.character(fst$CHR)

# assign alternating grey shades to scaffolds
unique_chr <- unique(fst$CHR)
grey_colors <- rep(c("#97979dff", "#323233ff"), length.out = length(unique_chr))
names(grey_colors) <- unique_chr

# override with purple for highlighted scaffolds
custom_colors <- grey_colors
custom_colors[highlighted_scaffolds] <- "#192bc2"

# assign point color to each row
fst$color <- custom_colors[fst$CHR]

# pPlot
fst_h <- ggplot(fst, aes(x = BPcum_Mb, y = WEIR_AND_COCKERHAM_FST, color = color)) +
  geom_point(size = 0.7) +
  labs(x = "Genomic position (Mb)",
    y = expression(italic(F)[ST])) +
  scale_color_identity() +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
  )
ggsave("Fst_ALL_for_MS_highlighted.png", fst_h, width = 15, height = 7, dpi = 300)

# CONTIG 29 #

setwd("")  

scaff <- "Phar_scaff_29"

xpehh_scaff <- subset(xpehh, CHR == scaff)
fst_w_scaff <- subset(fst, CHR == scaff)
  
# --- XP-EHH with annotated genes (top panel) --- #

# gene positions in bp
genes <- data.frame(
  name   = c("BMPR2", "CDK2", "ORC4", "COL4A1", "ACVR2B", "COL2A1"),
  start  = c(428518, 445455, 454017, 495764, 588516, 831379),
  end    = c(444279, 450428, 471755, 572314, 615145, 884526),
  strand = c("+", "-", "-", "-", "-", "+")
)

# create a palette of green shades
green_palette <- colorRampPalette(c("#aae9d1ff", "#247b7b"))(nrow(genes))

# create named vector for scale_fill_manual
fill_vector <- setNames(green_palette, genes$name)
fill_vector

# get scaffold offset in bp (for Phar_scaff_29)
offset_bp <- offsets["Phar_scaff_29"]

# convert to Mb for plotting
genes$start_Mb <- (genes$start + offset_bp) / 1e6
genes$end_Mb   <- (genes$end   + offset_bp) / 1e6

# function to make polygon coordinates for genome-style arrow
make_arrow <- function(start, end, y, height, strand) {
  body_height <- height * 0.8
  tip_width <- (end - start) * 0.1
  
  if (strand == "-") {
    data.frame(
      x = c(end, end, start + tip_width, start, start + tip_width, end),
      y = c(y, y + body_height/2, y + body_height/2, y, y - body_height/2, y - body_height/2)
    )
  } else {
    data.frame(
      x = c(start, start, end - tip_width, end, end - tip_width, start),
      y = c(y, y + body_height/2, y + body_height/2, y, y - body_height/2, y - body_height/2)
    )
  }
}

arrow_polygons <- do.call(rbind,
  lapply(1:nrow(genes), function(i) {
    arrow <- make_arrow(genes$start_Mb[i], genes$end_Mb[i],
                        y = max(xpehh_scaff$normxpehh) + 0.5 + (i-1)*0.5,
                        height = 0.4, strand = genes$strand[i])
    arrow$gene <- genes$name[i]
    arrow
  })
)

genes$label_x <- ifelse(genes$strand == "-",
                        genes$start_Mb - 0.02, 
                        genes$end_Mb + 0.02)    
genes$label_y <- max(xpehh_scaff$normxpehh) + 0.5 + (0:(nrow(genes)-1))*0.5

# plot with ggplot2
xpehh_plot <- ggplot(xpehh_scaff, aes(x = BPcum_Mb, y = normxpehh)) +
  geom_point(size = 0.9, color = "black") +
  geom_hline(yintercept = outl_abs, color = "#626269ff", linetype = "dashed") +
  geom_polygon(data = arrow_polygons,
               aes(x = x, y = y, group = gene, fill = gene),
               inherit.aes = FALSE, color = "black") +
  geom_text(data = genes,
            aes(x = label_x,
                y = label_y,
                label = name,
                hjust = ifelse(strand == "-", 1, 0)), 
            inherit.aes = FALSE, color = "black", fontface = "bold", size = 4) +
  labs(x = "", y = "Normalized XP-EHH", title = "CONTIG 29") +
  scale_fill_manual(values = fill_vector) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.25, vjust = 0.5, color = "#192bc2"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)  

# --- FST --- #

fst_plot <- ggplot() +
  geom_point(data = fst_w_scaff, 
             aes(x = BPcum_Mb, y = WEIR_AND_COCKERHAM_FST), 
             color = "black", size = 0.9) +
  geom_hline(yintercept = outlier_threshold, linetype = "dashed", color = "#626269ff") +
  labs(x = "Genomic position (Mb)", 
       y = expression(italic(F)[ST])) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)   

plot29 <- ggarrange(xpehh_plot, fst_plot, ncol = 1, nrow = 2)
ggsave(paste0("ANNOTATED_XPEHH_Fst_", scaff, ".png"), plot29, width = 20, height = 18, dpi = 300)
ggsave(paste0("ANNOTATED_XPEHH_Fst_", scaff, ".pdf"), plot29, width = 20, height = 18, dpi = 300)

## for PLOT BELOW ##
x_pos_sc29 <- mean(range(xpehh_scaff$BPcum_Mb))
y_pos_sc29 <- max(xpehh_scaff$normxpehh) * 0.9

# PHAR SCAFF 62 #

scaff <- "Phar_scaff_62"

xpehh_scaff <- subset(xpehh, CHR == scaff)
fst_w_scaff <- subset(fst, CHR == scaff)
  
# --- XP-EHH with annotated genes (top panel) --- #

# gene positions in bp
genes <- data.frame(
  name   = c("ANKS6", "PDCD6", "KDELC1", "DCAF11", "HNF4A", "ANO10", "ZDHHC7", "RSPH10B", "RAB22A", "TTLL9", "TUBB4B", "PDRG1", "RBM18", "SSNA1", "FBXO31", "DYNC1LI1", "eva-1"),
  start  = c(218235, 247124, 262148, 270716, 290236, 334044, 353102, 361705, 394066, 402751, 429261, 448004, 479501, 490966, 495716, 629465, 695873),
  end    = c(243888, 259508, 263653, 286623, 304431, 347854, 361617, 367467, 402688, 427513, 433184, 463477, 491599, 495384, 517009, 644137, 703413),
  strand = c("+", "-", "+", "+", "-", "+", "-", "+", "-", "+", "-", "-", "+", "-", "+", "+", "+")  
)

# create a palette of 9 green shades
green_palette <- colorRampPalette(c("#aae9d1ff", "#247b7b"))(nrow(genes))

# create named vector for scale_fill_manual
fill_vector <- setNames(green_palette, genes$name)
fill_vector

# get scaffold offset in bp (for Phar_scaff_29)
offset_bp <- offsets["Phar_scaff_62"]

# convert to Mb for plotting
genes$start_Mb <- (genes$start + offset_bp) / 1e6
genes$end_Mb   <- (genes$end   + offset_bp) / 1e6

# function to make polygon coordinates for genome-style arrow
make_arrow <- function(start, end, y, height, strand) {
  body_height <- height * 0.8
  tip_width <- (end - start) * 0.1
  
  if (strand == "-") {
    data.frame(
      x = c(end, end, start + tip_width, start, start + tip_width, end),
      y = c(y, y + body_height/2, y + body_height/2, y, y - body_height/2, y - body_height/2)
    )
  } else {
    data.frame(
      x = c(start, start, end - tip_width, end, end - tip_width, start),
      y = c(y, y + body_height/2, y + body_height/2, y, y - body_height/2, y - body_height/2)
    )
  }
}

arrow_polygons <- do.call(rbind,
  lapply(1:nrow(genes), function(i) {
    arrow <- make_arrow(genes$start_Mb[i], genes$end_Mb[i],
                        y = max(xpehh_scaff$normxpehh) + 0.5 + (i-1)*0.5,
                        height = 0.4, strand = genes$strand[i])
    arrow$gene <- genes$name[i]
    arrow
  })
)

genes$label_x <- ifelse(genes$strand == "-",
                        genes$start_Mb - 0.03, 
                        genes$end_Mb + 0.03)    
genes$label_y <- max(xpehh_scaff$normxpehh) + 0.5 + (0:(nrow(genes)-1))*0.6

# plot with ggplot2
xpehh_plot <- ggplot(xpehh_scaff, aes(x = BPcum_Mb, y = normxpehh)) +
  geom_point(size = 0.9, color = "black") +
  geom_hline(yintercept = outl_abs, color = "#626269ff", linetype = "dashed") +
  geom_polygon(data = arrow_polygons,
               aes(x = x, y = y, group = gene, fill = gene),
               inherit.aes = FALSE, color = "black") +
  geom_text(data = genes,
            aes(x = label_x,
                y = label_y,
                label = name,
                hjust = ifelse(strand == "-", 1, 0)), 
            inherit.aes = FALSE, color = "black", fontface = "bold", size = 4) +
  labs(x = "", y = "Normalized XP-EHH", title = "CONTIG 62") +
  scale_fill_manual(values = fill_vector) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.25, vjust = 0.5, color = "#192bc2"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)  

# --- FST --- #

fst_plot <- ggplot() +
  geom_point(data = fst_w_scaff, 
             aes(x = BPcum_Mb, y = WEIR_AND_COCKERHAM_FST), 
             color = "black", size = 0.9) +
  geom_hline(yintercept = outlier_threshold, linetype = "dashed", color = "#626269ff") +
  labs(x = "Genomic position (Mb)", 
       y = expression(italic(F)[ST])) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black"),
    panel.grid = element_blank()
)

plot62 <- ggarrange(xpehh_plot, fst_plot, ncol = 1, nrow = 2)
ggsave(paste0("ANNOTATED_XPEHH_Fst_", scaff, ".png"), plot62, width = 20, height = 18, dpi = 300)
ggsave(paste0("ANNOTATED_XPEHH_Fst_", scaff, ".pdf"), plot62, width = 20, height = 18, dpi = 300)

## for PLOT BELOW ##
x_pos_sc62 <- mean(range(xpehh_scaff$BPcum_Mb))
y_pos_sc62 <- max(xpehh_scaff$normxpehh) * 0.9

# -------- MANUSCRIPT PLOT HIGHLIGHTED XPEHH & FST TOGETHER WITH PHAR SCAFF 29 --------- #

p1 <- xpehh_h # from above
p2 <- fst_h # from above
p3 <- plot29 # from above, make sure to execute the right subset (Phar_scaff_29)
p4 <- plot62 # from above, make sure to execute the right subset (Phar_scaff_62)

top <- p1 / p2
bottom <- p3 | p4

p <- top / bottom + plot_layout(heights = c(1, 1, 2))

ggsave("1_MS_plot.png", p, dpi = 700, width = 20, height = 25)
