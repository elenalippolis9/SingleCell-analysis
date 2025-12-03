cran_pkgs <- c("dplyr","readr","stringr","tibble","purrr")
to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("SummarizedExperiment")
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc)) BiocManager::install(to_install_bioc, update = TRUE, ask = FALSE)
BiocManager::install(c(
  "recount3",          
  "recount",           
  "edgeR",
  "SummarizedExperiment",
  "MatrixGenerics"
))

library(recount3)
library(recount)
library(edgeR)

rse_brain <- readRDS("./rse_brain (1).RDS")
rse_liver <- readRDS("./rse_liver (1).RDS")
rse_pancreas <- readRDS("./rse_pancreas (1).RDS")

assays(rse_brain)$counts <- transform_counts(rse_brain)
assays(rse_liver)$counts <- transform_counts(rse_liver)
assays(rse_pancreas)$counts <- transform_counts(rse_pancreas)

rse_brain_selected <- rse_brain[,c(30, 40, 41)]
rse_liver_selected <- rse_liver[,c(33,34,36)]
rse_pancreas_selected <- rse_pancreas[,c(30,33, 34)]
counts_brain_selected <- assays(rse_brain_selected)$counts
counts_liver_selected <- assays(rse_liver_selected)$counts
counts_pancreas_selected <- assays(rse_pancreas_selected)$counts

x <- cbind(counts_brain_selected,counts_liver_selected,counts_pancreas_selected)
colnames(x) <- c("Brain30", "Brain40", "Brain41", "Liver33", "Liver34","Liver36","Pancreas30","Pancreas32", "Pancreas33")
rownames(x) <- rowData(rse_brain_selected)$gene_name
y <- DGEList(counts=x)

group <- as.factor(c("Brain","Brain", "Brain", "Liver","Liver","Liver", "Pancreas","Pancreas","Pancreas"))
y$samples$group <- group

y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_liver_selected)$gtex.smrin,colData(rse_pancreas_selected)$gtex.smrin))
y$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,colData(rse_liver_selected)$gtex.smtsd,colData(rse_pancreas_selected)$gtex.smtsd))
y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,colData(rse_liver_selected)$gtex.sex,colData(rse_pancreas_selected)$gtex.sex))
y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,colData(rse_liver_selected)$gtex.age,colData(rse_pancreas_selected)$gtex.age))
y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_liver_selected)$gtex.smrrnart,colData(rse_pancreas_selected)$gtex.smrrnart))
y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(rse_liver_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_pancreas_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))
y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_liver_selected)$"recount_qc.aligned_reads%.chrm",colData(rse_pancreas_selected)$"recount_qc.aligned_reads%.chrm"))
y$samples

table(rowSums(y$counts==0)==9)

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

dim(y)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(viridisLite)   
})

outdir <- "./final"


save_plot <- function(p, name, w=9, h=6, dpi=300) {
  fn_png <- file.path(outdir, paste0(name, ".png"))
  fn_pdf <- file.path(outdir, paste0(name, ".pdf"))
  ggsave(fn_png, p, width=w, height=h, dpi=dpi, bg="white")
  ggsave(fn_pdf, p, width=w, height=h, device=cairo_pdf, bg="white")
  message("Saved: ", fn_png, " & ", fn_pdf)
}

theme_pub <- theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title        = element_text(face="bold", size=13, margin=margin(b=7)),
    axis.title        = element_text(face="bold"),
    axis.text.x       = element_text(angle=40, hjust=1, vjust=1),
    legend.position   = "right",
    legend.title      = element_text(face="bold")
  )
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9", "#000000")
groups <- factor(y$samples$group)
pal    <- setNames(okabe_ito[seq_len(nlevels(groups))], levels(groups))

## Labels
nf_txt     <- sprintf("NF=%.3f", y$samples$norm.factors)
sample_id  <- colnames(y$counts)
sample_lab <- paste0(sample_id, "\n", nf_txt)

mk_box <- function(mat, title) {
  df <- as.data.frame(mat)
  df$gene <- rownames(mat)
  long <- reshape(
    df, direction="long", varying = colnames(df)[colnames(df)!="gene"],
    v.names="logCPM", timevar="sample", times=colnames(mat)
  )
  long$group <- groups[match(long$sample, colnames(y$counts))]
  long$sample <- factor(long$sample, levels=colnames(y$counts), labels = sample_lab)
  
  ggplot(long, aes(sample, logCPM, fill=group, color=group)) +
    geom_boxplot(outlier.shape = NA, alpha=0.7, width=0.65) +
    coord_cartesian(ylim = quantile(long$logCPM, c(0.02, 0.98))) + # trims extreme whiskers for readability
    scale_fill_manual(values=pal) + scale_color_manual(values=pal) +
    labs(x="Samples (TMM factors in labels)", y="log2 CPM", fill="Group", color="Group", title=title) +
    theme_pub
}

p_before <- mk_box(logcpm_before, "logCPM — BEFORE TMM")
p_after  <- mk_box(logcpm_after,  "logCPM — AFTER TMM")

save_plot(p_before, "boxplot_logCPM_before_pretty", w=11, h=6)
save_plot(p_after,  "boxplot_logCPM_after_pretty",  w=11, h=6)
save_plot(p_before + p_after + plot_layout(ncol=2, widths=c(1,1)) & theme(legend.position="none"),
          "boxplots_logCPM_before_after_pretty", w=14, h=6)

#MDS
logcpm <- cpm(y, log=TRUE)
mds <- plotMDS(logcpm, plot=FALSE)
mds_df <- data.frame(
  sample = sample_id,
  x      = mds$x, y = mds$y,
  group  = groups,
  NF     = nf_txt
)

p_mds <- ggplot(mds_df, aes(x, y, color = group, fill = group)) +
  stat_ellipse(type = "t", level = 0.68, alpha = 0.15, show.legend = FALSE) +
  geom_point(size = 3) +
  ## --- LABELS: choose ONE of the blocks below (A or B) ---
  # (B) NO arrows at all  (uncomment this block and remove block A)
  geom_text_repel(
    aes(label = sprintf("%s\nNF=%s", sample, sub("NF=","", NF))),
    size = 3, box.padding = 0.5, point.padding = 0.7,
    force = 2, force_pull = 0.6, max.overlaps = Inf,
    segment.alpha = 0,             # <- hides the arrows
    seed = 42
  ) +
  scale_color_manual(values = pal, name = "Group") +
  scale_fill_manual(values  = pal, guide = "none") +  # <- removes the extra legend
  labs(x = "Dimension 1 (leading logFC)", y = "Dimension 2 (leading logFC)",
       title = "MDS of samples (labels + TMM normalization factors)") +
  theme_pub +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 28, 10, 10))




save_plot(p_mds, "mds_samples_pretty", w=8.5, h=6.5)

# Build MDS data with all metadata + label_all
# Coordinates
mds <- plotMDS(logcpm, plot = FALSE)

# Sample ids in the same order as the columns of y$counts / mds$x
sample_id <- colnames(y$counts)

# Minimal, explicit metadata table
meta <- data.frame(
  sample = sample_id,
  group  = y$samples$group,
  rRNA   = y$samples$rRNA,
  chrm   = y$samples$chrm,
  age    = y$samples$age,
  stringsAsFactors = FALSE
)

# Join coords + meta (preserve order)
mds_df <- data.frame(
  sample = sample_id,
  x = as.numeric(mds$x),
  y = as.numeric(mds$y),
  stringsAsFactors = FALSE
)
mds_df <- merge(mds_df, meta, by = "sample", sort = FALSE)

# Nice number formatter
fmt_num <- function(x) format(round(as.numeric(x), 3), nsmall = 3, trim = TRUE)

# One clean, multi-line label per point
mds_df$label_all <- sprintf(
  "%s\nrRNA=%s  chrM=%s  age=%s",
  mds_df$sample, fmt_num(mds_df$rRNA), fmt_num(mds_df$chrm), as.character(mds_df$age)
)
# gene annotation access 
nG <- nrow(y$counts)

get_col <- function(df, candidates, default = NULL) {
  if (is.null(df)) return(if (is.null(default)) rep(NA, nG) else default)
  for (nm in candidates) if (nm %in% colnames(df)) return(df[[nm]])
  if (is.null(default)) rep(NA, nG) else default
}

genes <- if (!is.null(y$genes)) y$genes else NULL

sym <- get_col(genes, c("gene_name","symbol","hgnc_symbol"), rownames(y$counts))
chr <- tolower(get_col(genes, c("chromosome","seqnames","chr","chromosome_name"),
                       rep(NA_character_, nG)))
bio <- tolower(get_col(genes, c("gene_biotype","biotype","feature_types"),
                       rep(NA_character_, nG)))

# flags (work even if sym/chr/bio contain duplicates)
is_mito <- (!is.na(sym) & startsWith(sym, "MT-")) |
  (!is.na(chr) & chr %in% c("mt","chrmt","m","chrm","mitochondrion"))
is_rrna <- (!is.na(bio) & grepl("rrna", bio)) |
  (!is.na(sym) & grepl("^RNR[0-9]|^MT-RNR", sym, ignore.case = TRUE))

# compute percentages from RAW counts
lib        <- colSums(y$counts)
rrna_counts <- if (any(is_rrna)) colSums(y$counts[is_rrna, , drop = FALSE]) else rep(0, ncol(y$counts))
mito_counts <- if (any(is_mito)) colSums(y$counts[is_mito, , drop = FALSE]) else rep(0, ncol(y$counts))

y$samples$rRNA <- round(100 * rrna_counts / lib, 3)
y$samples$chrm <- round(100 * mito_counts / lib, 2)


# ensure types
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# minimal theme
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(),
    legend.position = "right"
  )

pal <- c(Brain = "#0072B2", Liver = "#E69F00", Pancreas = "#009E73")

# compute MDS coords once
mds <- plotMDS(logcpm, plot = FALSE)
sample_id <- colnames(y$counts)

mds_df <- data.frame(
  sample = sample_id,
  x = as.numeric(mds$x),
  y = as.numeric(mds$y),
  group = y$samples$group,
  age   = factor(y$samples$age),
  rRNA  = as.numeric(y$samples$rRNA),
  chrm  = as.numeric(y$samples$chrm),
  stringsAsFactors = FALSE
)

# label formatters
fmt_txt <- function(x) as.character(x)
fmt_num <- function(x) format(round(as.numeric(x), 1), nsmall = 1, trim = TRUE)

# one nice MDS maker
mds_png_nice <- function(label_vec, title, filename,
                         label_formatter = fmt_txt,
                         width = 1600, height = 1200, res = 200) {

  df <- mds_df
  df$label <- label_formatter(label_vec)

  # pick filled shapes for age so outline&fill show; up to 5 bins
  age_lvls <- levels(df$age)
  age_shapes <- setNames(21:(21 + length(age_lvls) - 1), age_lvls)

  p <- ggplot(df, aes(x, y)) +
    stat_ellipse(aes(color = group), type = "t", level = 0.68,
                 alpha = 0.10, linewidth = 0.6, show.legend = FALSE) +
    geom_point(aes(color = group, shape = age), size = 3.3, stroke = 0.9) +
    geom_text_repel(
      aes(label = label, color = group),
      size = 3.5,
      direction = "y",
      nudge_y = 0.18,
      box.padding = 0.18,
      point.padding = 0.4,
      force = 0.3, force_pull = 0,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.3,
      segment.alpha = 0.7,
      seed = 42,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, name = "Group") +
    scale_shape_manual(values = age_shapes, name = "Age") +
    scale_x_continuous(expand = expansion(mult = c(.06, .10))) +
    scale_y_continuous(expand = expansion(mult = c(.06, .12))) +
    labs(x = "Dimension 1 (leading logFC)",
         y = "Dimension 2 (leading logFC)",
         title = title) +
    theme_pub +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 26, 10, 10))

  png(file.path(outdir, filename), width = width, height = height, res = res)
  print(p)
  dev.off()
}

library(edgeR)
library(ggplot2)
library(ggrepel)

# Align rows of y$samples to columns of y$counts
y$samples <- y$samples[colnames(y$counts), , drop = FALSE]

# Safe factor→numeric helper
numize <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Build plotting frame with CORRECT numeric values
logcpm <- cpm(y, log = TRUE)
mds    <- plotMDS(logcpm, plot = FALSE)

smeta <- y$samples
rRNA_raw <- numize(smeta$rRNA)
chrM_raw <- numize(smeta$chrm)

mds_df <- data.frame(
  sample = colnames(y$counts),
  x      = as.numeric(mds$x),
  y      = as.numeric(mds$y),
  group  = smeta$group,
  NF     = numize(smeta$norm.factors),
  rRNA   = ifelse(!is.na(rRNA_raw) & rRNA_raw <= 1, 100 * rRNA_raw, rRNA_raw),
  chrM   = chrM_raw,                 
  age    = smeta$age,
  sex    = smeta$sex,
  stringsAsFactors = FALSE
)

#Plots
fmt1 <- function(x) format(round(as.numeric(x), 1), nsmall = 1, trim = TRUE)
pal  <- c(Brain="#0072B2", Liver="#E69F00", Pancreas="#009E73")

theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face="bold", size=13, margin=margin(b=6)),
        axis.title = element_text(face="bold"),
        legend.position = "right")

# Labels = NF; shape = sex
p_mds_nf <- ggplot(mds_df, aes(x, y)) +
  stat_ellipse(aes(color = group), type = "t", level = 0.68,
               alpha = 0.12, linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(color = group, shape = sex), size = 3.2, stroke = 0.9) +
  geom_text_repel(
    aes(label = sprintf("%s\nNF=%.3f", sample, NF), color = group),
    size = 3.2, box.padding = 0.25, point.padding = 0.45,
    max.overlaps = Inf, segment.size = 0.25, segment.alpha = 0.6, seed = 42
  ) +
  scale_color_manual(values = pal, name = "Group") +
  scale_shape_discrete(name = "Sex") +
  labs(x = "Dimension 1 (leading logFC)", y = "Dimension 2 (leading logFC)",
       title = "MDS of samples (TMM; labels = NF; shape = sex)") +
  theme_pub

# Labels = rRNA/chrM/age/sex; shape = sex
p_mds_meta <- ggplot(mds_df, aes(x, y)) +
  stat_ellipse(aes(color = group), type = "t", level = 0.68,
               alpha = 0.12, linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(color = group, shape = sex), size = 3.2, stroke = 0.9) +
  geom_text_repel(
    aes(label = sprintf("%s\nrRNA=%s%%  chrM=%s%%  age=%s  sex=%s",
                        sample, fmt1(rRNA), fmt1(chrM), as.character(age), as.character(sex)),
        color = group),
    size = 3.1, box.padding = 0.22, point.padding = 0.42,
    max.overlaps = Inf, segment.size = 0.25, segment.alpha = 0.6, seed = 42
  ) +
  scale_color_manual(values = pal, name = "Group") +
  scale_shape_discrete(name = "Sex") +
  labs(x = "Dimension 1 (leading logFC)", y = "Dimension 2 (leading logFC)",
       title = "MDS of samples with rRNA/chrM/age/sex") +
  theme_pub

# Save 
if (exists("save_plot")) {
  save_plot(p_mds_nf,   "mds_samples_w_sex_nf",   w=8.5, h=6.5, dpi=300)
  save_plot(p_mds_meta, "mds_samples_w_metadata", w=9.0, h=7.0,  dpi=300)
} else {
  ggsave("mds_samples_w_sex_nf.png",   p_mds_nf,   width=8.5, height=6.5, dpi=300, bg="white")
  ggsave("mds_samples_w_metadata.png", p_mds_meta, width=9.0, height=7.0,  dpi=300, bg="white")
}



y <- estimateDisp(y, design)
png(file.path(outdir, "bcv.png"), width = 1600, height = 1200, res = 200)
plotBCV(y, main = "Biological Coefficient of Variation")
dev.off()
pdf(file.path(outdir, "bcv2.pdf"), width = 8, height = 6)
plotBCV(y, main = "Biological Coefficient of Variation")
dev.off()


#BCV plot 
if (is.null(y$trended.dispersion) || is.null(y$tagwise.dispersion)) {
  y <- estimateDisp(y, model.matrix(~0+group, data=y$samples))
}

bcv_df <- data.frame(
  Amean = y$AveLogCPM,
  BCV_tag = sqrt(y$tagwise.dispersion),
  BCV_trend = sqrt(y$trended.dispersion)
)

p_bcv <- ggplot(bcv_df, aes(Amean, BCV_tag)) +
  geom_point(alpha=0.15, size=0.8) +
  geom_line(aes(y=BCV_trend), linewidth=0.9) +
  labs(x="Average log2 CPM", y="Biological Coefficient of Variation",
       title="BCV vs expression (points: tagwise; line: trended)") +
  theme_pub

save_plot(p_bcv, "bcv_pretty", w=8, h=6)

nf_table <- data.frame(
  sample       = sample_id,
  group        = as.character(groups),
  lib.size     = y$samples$lib.size,
  norm.factor  = y$samples$norm.factors,
  eff.lib.size = y$samples$lib.size * y$samples$norm.factors,
  stringsAsFactors = FALSE
)
nf_table
write.csv(nf_table, file.path(outdir, "normalization_factors.csv"), row.names = FALSE)
message("Exported normalization_factors.csv")

# RIPRENDI DI QUA

#liver (top) vs brain (bottom)
qlfLB <- glmQLFTest(fit, contrast=c(-1,1,0))
#pancreas (top) vs brain (bottom)
qlfPB <- glmQLFTest(fit, contrast=c(-1,0,1))
#pancreas (top) vs liver (bottom)
qlfPL <- glmQLFTest(fit, contrast=c(0,-1,1))

qlfLB

head(qlfLB$table)

# top DE genes between Brain vs Liver
topTags(qlfLB, n=10,adjust.method = "BH", sort.by = "PValue")

resultsLB <- topTags(qlfHB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

resultsLB


summary(decideTests(qlfLB, p.value=0.05, adjust.method = "BH", lfc=0))
summary(decideTests(qlfLB, p.value=0.05, adjust.method = "BH", lfc=1))
summary(decideTests(qlfLB, p.value=0.01, adjust.method = "BH", lfc=0))
summary(decideTests(qlfLB, p.value=0.01, adjust.method = "BH",lfc=1))
resultsPB <- topTags(qlfPB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsPL <- topTags(qlfPL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)

resultsPB
resultsPL

assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_liver)$TPM <- recount::getTPM(rse_liver)
assays(rse_pancreas)$TPM <- recount::getTPM(rse_pancreas)
which(rowData(rse_brain)$gene_name == "ARG1")

boxplot(assays(rse_brain)$TPM[44926,],assays(rse_liver)$TPM[44926,], assays(rse_pancreas)$TPM[44926,], outline=F , names=c("Brain","Liver","Pancreas"))

safe_boxplot_save <- function(brain, liver, panc, outdir = "plots",
                              fname = "gene_44926_boxplot", dpi = 200) {
  # Ensure folder exists
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  png_path <- file.path(outdir, paste0(fname, ".png"))
  pdf_path <- file.path(outdir, paste0(fname, ".pdf"))
  
  # ---- PNG ----
  grDevices::png(png_path, width = 1600, height = 1200, res = dpi)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  boxplot(
    brain, liver, panc,
    outline = FALSE,
    names = c("Brain", "Liver", "Pancreas"),
    ylab = "TPM",
    main = "Expression of gene ARG1"
  )
  
  grDevices::dev.off()  # close PNG device
  
  grDevices::pdf(pdf_path, width = 8, height = 6)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  boxplot(
    brain, liver, panc,
    outline = FALSE,
    names = c("Brain", "Liver", "Pancreas"),
    ylab = "TPM",
    main = "Expression of gene ARG1"
  )
  
  grDevices::dev.off()  # close PDF device
  
  message("Saved:\n  ", normalizePath(png_path), "\n  ", normalizePath(pdf_path))
  invisible(list(png = png_path, pdf = pdf_path))
}

# Call it with your data:
safe_boxplot_save(
  brain = assays(rse_brain)$TPM[44926,],
  liver = assays(rse_liver)$TPM[44926,],
  panc  = assays(rse_pancreas)$TPM[44926,],
  outdir = outdir,
  fname = "gene_31286_boxplot"
)

file.exists(file.path(outdir, "gene_31286_boxplot.png"))
file.exists(file.path(outdir, "gene_31286_boxplot.pdf"))

save_colorful_boxplot <- function(brain, liver, panc,
                                  outdir = "plots",
                                  fname  = "gene_ARG1_boxplot_color",
                                  dpi = 300) {
  
  # Make sure the output folder exists
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Colors
  bg_col   <- rgb(255, 255, 255, maxColorValue = 255)  
  box_cols <- c(Brain = "#0072B2", Liver = "#E69F00", Pancreas = "#009E73")  # colorblind-safe
  border_cols <- adjustcolor(box_cols, offset = c(-0.1, -0.1, -0.1, 0))     # slightly darker borders
  
  # Small helper to draw the same plot on any device
  draw_plot <- function() {
    par(bg = bg_col, mar = c(4.5, 5, 3, 2) + 0.1, family = "sans")
    # light grid
    yr <- range(c(brain, liver, panc), na.rm = TRUE)
    plot.new(); plot.window(xlim = c(0.5, 3.5), ylim = yr)
    abline(h = pretty(yr), col = adjustcolor("grey30", 0.25), lwd = 0.6)
    # actual boxplot
    boxplot(
      list(Brain = brain, Liver = liver, Pancreas = panc),
      add = TRUE, outline = FALSE,
      col = box_cols, border = border_cols,
      boxwex = 0.5, lwd = 1.1, whisklty = 1, staplelty = 1,
      yaxt = "n", xaxt = "n"
    )
    axis(2, las = 1)
    axis(1, at = 1:3, labels = names(box_cols), tick = FALSE, cex.axis = 1.0)
    title(main = "Expression of gene arg1", xlab = "", ylab = "TPM", font.main = 2)
    stripchart(
      list(brain, liver, panc), method = "jitter", vertical = TRUE, add = TRUE,
      pch = 16, cex = 0.6, col = adjustcolor("black", 0.35), jitter = 0.08
    )
    box(lwd = 1)  # frame
  }
  
  # PNG 
  png_path <- file.path(outdir, paste0(fname, ".png"))
  png(png_path, width = 1600, height = 1200, res = dpi, bg = bg_col)
  draw_plot()
  dev.off()
  
  # pdf
  pdf_path <- file.path(outdir, paste0(fname, ".pdf"))
  pdf(pdf_path, width = 8, height = 6, bg = bg_col, family = "sans")
  draw_plot()
  dev.off()
  
  message("Saved:\n  ", normalizePath(png_path), "\n  ", normalizePath(pdf_path))
  invisible(list(png = png_path, pdf = pdf_path))
}

# Call it with your data:
save_colorful_boxplot(
  brain = assays(rse_brain)$TPM[44926,],
  liver = assays(rse_liver)$TPM[44926,],
  panc  = assays(rse_pancreas)$TPM[44926,],
  outdir = outdir,
  fname  = "gene_arg1_boxplot_color"
)
install.packages(c("ggtext","glue"))

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggtext)
})

# parameters 
alpha   <- 0.01   # FDR threshold
lfc_min <- 1      # absolute log2FC threshold
cpm_min <- 0      # mean logCPM threshold
bg_col  <- rgb(255,255,255, maxColorValue = 255)  # your RGB background

#helpers 

# Extract a tidy results table from a DGELRT (qlf*) object
tidy_from_qlf <- function(qlf) {
  tt <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
  tibble::rownames_to_column(as.data.frame(tt), var = "gene")
}

# Call Up/Down/NotSig given thresholds
call_de <- function(df, alpha, lfc_min, cpm_min) {
  df %>%
    mutate(
      sig  = PValue,  # already used to compute FDR in topTags
      FDR  = p.adjust(PValue, method="BH"),
      call = case_when(
        FDR < alpha & logCPM > cpm_min & logFC >=  lfc_min ~ "Up",
        FDR < alpha & logCPM > cpm_min & logFC <= -lfc_min ~ "Down",
        TRUE ~ "NotSig"
      )
    )
}

#  theme
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = bg_col, color = NA),
    panel.background = element_rect(fill = bg_col, color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=8)),
    axis.title = element_text(face = "bold")
  )

pal <- c(Brain = "#1F77B4", Liver = "#E69F00", Pancreas = "#2CA02C")

#build DE calls for each contrast 
res_LB <- tidy_from_qlf(qlfLB) %>% call_de(alpha, lfc_min, cpm_min) # Liver vs Brain (positive= L> B)
res_PB <- tidy_from_qlf(qlfPB) %>% call_de(alpha, lfc_min, cpm_min) # Pancreas vs Brain (positive= P> B)
res_PL <- tidy_from_qlf(qlfPL) %>% call_de(alpha, lfc_min, cpm_min) # Pancreas vs Liver (positive= P> L)

#PLOT 1: three Up/Down tiles per pair


pair_counts <- bind_rows(
  "Liver vs Brain"    = res_LB,
  "Pancreas vs Brain" = res_PB,
  "Pancreas vs Liver" = res_PL,
  .id = "pair"
) %>%
  dplyr::count(pair, call, name = "n") %>%          
  filter(call %in% c("Up", "Down")) %>%
  mutate(call = factor(call, levels = c("Up", "Down")))


# tile plot (compact table style)
p_pairs <- ggplot(pair_counts, aes(x = call, y = pair, fill = call, label = scales::comma(n))) +
  geom_tile(color = "white", height = 0.8, width = 0.95) +
  geom_text(fontface = "bold", size = 5) +
  scale_fill_manual(values = c(Up = "#2CA02C", Down = "#A4161A")) +
  labs(title = glue::glue("Pairwise DE calls  •  FDR<{alpha}, |log2FC|≥{lfc_min}, mean logCPM>{cpm_min}"),
       x = NULL, y = NULL) +
  theme_pub +
  theme(legend.position="none")

ggsave(file.path(outdir, "pairwise_de_tiles.png"), p_pairs, width=10, height=3.2, dpi=300, bg=bg_col)
ggsave(file.path(outdir, "pairwise_de_tiles.pdf"), p_pairs, width=10, height=3.2, device="pdf", bg=bg_col)

## tissue-specific UP (vs both others)

# define UP direction for each contrast relative to the first tissue named
up_L_over_B <- res_LB %>% filter(call=="Up")   %>% select(gene) %>% mutate(tissue="Liver")
up_P_over_B <- res_PB %>% filter(call=="Up")   %>% select(gene) %>% mutate(tissue="Pancreas")
up_P_over_L <- res_PL %>% filter(call=="Up")   %>% select(gene) %>% mutate(tissue="Pancreas")

# Tissue-specific up genes:
# Brain-specific up = (Brain > Liver) AND (Brain > Pancreas)
brain_up   <- intersect( res_LB %>% filter(call=="Down") %>% pull(gene),   # Brain > Liver
                         res_PB %>% filter(call=="Down") %>% pull(gene) )  # Brain > Pancreas

# Liver-specific up = (Liver > Brain) AND (Liver > Pancreas)
liver_up   <- intersect( res_LB %>% filter(call=="Up") %>% pull(gene),     # Liver > Brain
                         res_PL %>% filter(call=="Down") %>% pull(gene) )  # Liver > Pancreas

# Pancreas-specific up = (Pancreas > Brain) AND (Pancreas > Liver)
pancreas_up<- intersect( res_PB %>% filter(call=="Up") %>% pull(gene),     # Pancreas > Brain
                         res_PL %>% filter(call=="Up") %>% pull(gene) )    # Pancreas > Liver

tissue_counts <- tibble::tibble(
  Tissue = c("Brain","Liver","Pancreas"),
  n_up   = c(length(brain_up), length(liver_up), length(pancreas_up))
)

p_tissue <- ggplot(tissue_counts, aes(Tissue, n_up, fill = Tissue)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = n_up), vjust = -0.4, fontface = "bold") +
  scale_fill_manual(values = pal) +
  labs(
    title = "Number of upregulated genes per tissue w.r.t. both other tissues",
    subtitle = glue::glue("Selection: FDR<{alpha}  •  |log2FC|≥{lfc_min}  •  mean logCPM>{cpm_min}   (total genes: {nrow(res_LB)})"),
    x = "Tissue", y = "Number of genes"
  ) +
  theme_pub

ggsave(file.path(outdir, "tissue_specific_up_bar.png"), p_tissue, width=9, height=6, dpi=300, bg=bg_col)
ggsave(file.path(outdir, "tissue_specific_up_bar.pdf"), p_tissue, width=9, height=6, device="pdf", bg=bg_col)

#export the actual gene lists
writeLines(brain_up,    file.path(outdir, "brain_specific_up.txt"))
writeLines(liver_up,    file.path(outdir, "liver_specific_up.txt"))
writeLines(pancreas_up, file.path(outdir, "pancreas_specific_up.txt"))


wilcox.test(assays(rse_brain)$TPM[44926,], assays(rse_heart)$TPM[44926,])
wilcox.test(assays(rse_brain)$TPM[44926,], assays(rse_pancreas)$TPM[44926,])

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)          # human gene IDs
  library(ReactomePA)
  library(enrichplot)
})

alpha   <- 0.01   # FDR cutoff (stringent)
lfc_min <- 1      # |log2FC| >= 1 (≈2x)
cpm_min <- 0      # average expression filter

call_de <- function(tt) {
  tt %>%
    mutate(FDR = p.adjust(PValue, "BH"),
           call = case_when(
             FDR < alpha & logCPM > cpm_min & logFC >=  lfc_min ~ "Up",
             FDR < alpha & logCPM > cpm_min & logFC <= -lfc_min ~ "Down",
             TRUE ~ "NotSig"))
}

tbl_LB <- call_de(topTags(qlfLB, n=Inf)$table %>% tibble::rownames_to_column("gene"))
tbl_PB <- call_de(topTags(qlfPB, n=Inf)$table %>% tibble::rownames_to_column("gene"))
tbl_PL <- call_de(topTags(qlfPL, n=Inf)$table %>% tibble::rownames_to_column("gene"))

# Brain > both (Brain higher than Liver AND Pancreas)
brain_up    <- intersect(tbl_LB$gene[tbl_LB$call=="Down"],  # Brain > Liver
                         tbl_PB$gene[tbl_PB$call=="Down"])  # Brain > Pancreas
# Liver > both
liver_up    <- intersect(tbl_LB$gene[tbl_LB$call=="Up"],    # Liver > Brain
                         tbl_PL$gene[tbl_PL$call=="Down"])  # Liver > Pancreas
# Pancreas > both
pancreas_up <- intersect(tbl_PB$gene[tbl_PB$call=="Up"],    # Pancreas > Brain
                         tbl_PL$gene[tbl_PL$call=="Up"])    # Pancreas > Liver


# Universe = all genes that were tested and pass the expression filter
universe_symbols <- tbl_LB$gene  # all tested genes (same order as others)
# Map to Entrez
sym2ent <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(universe_symbols),
                                 keytype = "SYMBOL", columns = "ENTREZID") %>%
  dplyr::filter(!is.na(ENTREZID)) %>% distinct(SYMBOL, .keep_all = TRUE)
map_to_entrez <- function(symbols) {
  as.character(na.omit(sym2ent$ENTREZID[match(symbols, sym2ent$SYMBOL)]))
}
universe_entrez <- map_to_entrez(universe_symbols)

genes_list <- list(
  Brain    = map_to_entrez(brain_up),
  Liver    = map_to_entrez(liver_up),
  Pancreas = map_to_entrez(pancreas_up)
)

do_enrich <- function(entrez_vec, universe_entrez) {
  list(
    GO_BP = enrichGO(gene         = entrez_vec,
                     universe     = universe_entrez,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",
                     pAdjustMethod= "BH",
                     qvalueCutoff = 0.05,
                     readable     = TRUE),
    KEGG  = enrichKEGG(gene       = entrez_vec,
                       universe   = universe_entrez,
                       organism   = "hsa",
                       pAdjustMethod= "BH",
                       qvalueCutoff = 0.05) |>
      setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    REACT = enrichPathway(gene    = entrez_vec,
                          universe= universe_entrez,
                          organism= "human",
                          pAdjustMethod= "BH",
                          qvalueCutoff = 0.05,
                          readable = TRUE)
  )
}

enrich_res <- lapply(genes_list, do_enrich, universe_entrez = universe_entrez)

dir.create(file.path(outdir, "enrichment"), showWarnings = FALSE, recursive = TRUE)

simplify_go <- function(go_obj) {
  if (is.null(go_obj) || nrow(go_obj@result)==0) return(go_obj)
  simplify(go_obj, cutoff = 0.6, by = "p.adjust", select_fun = min)
}

plot_and_save <- function(obj, tissue, db) {
  if (is.null(obj) || nrow(obj@result)==0) return(invisible(NULL))
  # table
  fn_csv <- file.path(outdir, "enrichment",
                      sprintf("%s_%s_enrichment.csv", tissue, db))
  readr::write_csv(as.data.frame(obj@result), fn_csv)
  
  # top-15 dotplot
  p <- dotplot(obj, showCategory = 15, font.size = 11) +
    ggtitle(sprintf("%s — %s (top 15)", tissue, db)) +
    theme_minimal(base_size = 12)
  ggsave(file.path(outdir, "enrichment",
                   sprintf("%s_%s_dotplot.png", tissue, db)),
         p, width = 8.5, height = 6, dpi = 300, bg = "white")
}

for (tissue in names(enrich_res)) {
  go_s <- simplify_go(enrich_res[[tissue]]$GO_BP)
  enrich_res[[tissue]]$GO_BP_simplified <- go_s
  plot_and_save(go_s,          tissue, "GO_BP")
  plot_and_save(enrich_res[[tissue]]$KEGG,  tissue, "KEGG")
  plot_and_save(enrich_res[[tissue]]$REACT, tissue, "Reactome")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)   # tableGrob
  library(cowplot)     # ggdraw + draw_label
})

# STYLE
bg_col <- rgb(230,218,209, maxColorValue = 255)  # your RGB
header_fill <- "#0C5A6B"       # teal for the header strip (change if you want)
header_txt  <- "white"

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = bg_col, colour = NA),
    panel.background = element_rect(fill = bg_col, colour = NA)
  )

# PANEL MAKER
make_enrich_panel <- function(
    tissue, cp_obj, show_n = 6,
    alpha = 0.01, lfc_min = 1, cpm_min = 0,
    up_n = NA_integer_, universe_n = NA_integer_
){
  stopifnot(!is.null(cp_obj))
  df <- as.data.frame(cp_obj@result)
  if (nrow(df) == 0) {
    return(ggdraw() + draw_label(sprintf("%s: no significant terms", tissue),
                                 fontface = "bold", size = 12))
  }
  
  # pick columns & shorten text
  top <- df[order(df$p.adjust), ][seq_len(min(show_n, nrow(df))), ]
  top$`GO Term`    <- top$ID
  top$Description  <- top$Description
  top$Category     <- if ("ONTOLOGY" %in% names(top)) top$ONTOLOGY else "BP"
  top$`Adjusted p-value` <- formatC(top$p.adjust, format = "e", digits = 2)
  tbl <- top[, c("GO Term","Description","Category","Adjusted p-value")]
  
  # tableGrob
  tg <- tableGrob(
    tbl, rows = NULL,
    theme = ttheme_minimal(
      core = list(
        fg_params = list(fontface = c("plain"), fontsize = 10, x = 0.02),
        bg_params = list(fill = c("grey95","white"), alpha = 1)
      ),
      colhead = list(
        fg_params = list(fontface = "bold", col = "white"),
        bg_params = list(fill = header_fill, alpha = 1)
      )
    )
  )
  
  # header banner text with parameters
  hdr <- sprintf(
    "%s • Selection: FDR<%g  |  |log2FC|≥%g  |  mean logCPM>%g  |  Up genes=%s  |  Universe=%s",
    tissue, alpha, lfc_min, cpm_min,
    ifelse(is.na(up_n), "–", up_n),
    ifelse(is.na(universe_n), "–", universe_n)
  )
  
  # assemble: header strip + table
  header_g <- ggdraw() +
    draw_rect(fill = header_fill, color = header_fill) +
    draw_label(hdr, colour = header_txt, size = 12, fontface = "bold",
               x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5)
  
  plot_grid(
    header_g, ggdraw(tg),
    ncol = 1, rel_heights = c(0.18, 1)
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)   # tableGrob(), ttheme_minimal
  library(grid)        # unit()
  library(cowplot)     # plot_grid(), ggdraw()
  library(stringr)     # str_wrap()
})

# One function that builds the panel + header (no free vars)
make_enrich_panel <- function(
    tissue,
    cp_obj,                 # clusterProfiler object (e.g., GO_BP_simplified)
    show_n      = 6,
    alpha       = 0.01,
    lfc_min     = 1,
    cpm_min     = 0,
    up_n        = NA_integer_,
    universe_n  = NA_integer_,
    header_fill = "#0C5A6B",
    header_txt  = "white",
    bg_col      = rgb(255,255,255, maxColorValue = 255),
    wrap_width  = 60         # characters to wrap the Description
) {
  # empty
  if (is.null(cp_obj) || nrow(cp_obj@result) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x=.5, y=.5,
                      label = sprintf("%s: no significant terms", tissue),
                      size = 5, fontface = "bold"))
  }
  
  df  <- as.data.frame(cp_obj@result)
  top <- df[order(df$p.adjust), , drop = FALSE]
  top <- top[seq_len(min(show_n, nrow(top))), , drop = FALSE]
  
  tab <- data.frame(
    `GO Term`          = top$ID,
    Description        = str_wrap(top$Description, width = wrap_width),
    Category           = if ("ONTOLOGY" %in% names(top)) top$ONTOLOGY else "BP",
    `Adjusted p-value` = formatC(top$p.adjust, format = "e", digits = 2),
    check.names = FALSE
  )
  
  # table theme: left-aligned + tighter padding
  tt <- ttheme_minimal(
    core = list(
      fg_params = list(hjust = 0, x = 0.02, fontsize = 11),
      bg_params = list(fill = rep(c("grey95","white"), length.out = nrow(tab))),
      padding  = unit(c(5, 8), "pt")
    ),
    colhead = list(
      fg_params = list(hjust = 0, x = 0.02, fontface = "bold", col = "white", fontsize = 12),
      bg_params = list(fill = header_fill),
      padding  = unit(c(6, 8), "pt")
    )
  )
  
  tg <- tableGrob(tab, rows = NULL, theme = tt)
  
  # widen description col so values aren’t clipped
  #            GO ID  Description  Category  Adj p
  tg$widths <- unit(c(0.20,  0.56,   0.10,    0.14), "npc")
  
  # header text
  hdr <- sprintf(
    "%s • Selection: FDR<%g  |  |log2FC|≥%g  |  mean logCPM>%g  |  Up genes=%s  |  Universe=%s",
    tissue, alpha, lfc_min, cpm_min,
    ifelse(is.na(up_n), "–", up_n),
    ifelse(is.na(universe_n), "–", universe_n)
  )
  
  header_g <- ggplot() +
    geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=1),
              fill = header_fill, colour = header_fill) +
    annotate("text", x=.5, y=.5, label = hdr,
             colour = header_txt, size = 4, fontface = "bold") +
    theme_void()
  
  panel <- cowplot::plot_grid(
    header_g,
    cowplot::ggdraw(tg),
    ncol = 1,
    rel_heights = c(0.10, 0.90)   # shorter header
  )
  
  panel + theme(plot.background = element_rect(fill = bg_col, colour = NA))
}  


p_brain <- make_enrich_panel(
  tissue     = "Brain",
  cp_obj     = enrich_res$Brain$GO_BP_simplified,
  show_n     = 5,
  alpha      = alpha,
  lfc_min    = lfc_min,
  cpm_min    = cpm_min,
  up_n       = length(brain_up),
  universe_n = length(universe_entrez)
)

ggsave(file.path(outdir, "enrichment_brain_GO_panel.png"),
       p_brain, width = 12, height = 5.5, dpi = 300,
       bg = rgb(255,255,255, maxColorValue = 255))

p_liver <- make_enrich_panel(
  tissue     = "Liver",
  cp_obj     = enrich_res$Liver$GO_BP_simplified,
  show_n     = 5,
  alpha      = alpha,
  lfc_min    = lfc_min,
  cpm_min    = cpm_min,
  up_n       = length(liver_up),
  universe_n = length(universe_entrez)
)

ggsave(file.path(outdir, "enrichment_liver_GO_panel.png"),
       p_liver, width = 12, height = 5.5, dpi = 300,
       bg = rgb(255,255,255, maxColorValue = 255))

p_pancreas <- make_enrich_panel(
  tissue     = "Pancreas",
  cp_obj     = enrich_res$Pancreas$GO_BP_simplified,
  show_n     = 5,
  alpha      = alpha,
  lfc_min    = lfc_min,
  cpm_min    = cpm_min,
  up_n       = length(pancreas_up),
  universe_n = length(universe_entrez)
)

ggsave(file.path(outdir, "enrichment_pancreas_GO_panel.png"),
       p_pancreas, width = 12, height = 5.5, dpi = 300,
       bg = rgb(255,255,255, maxColorValue = 255))

# run GO for ALL ontologies
ego_all <- enrichGO(
  gene          = genes_list$Brain,   # Entrez IDs for the tissue
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",               # <— BP, CC, MF together
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# simplify separately within each ontology (avoids cross-ontology merging)
by_ont <- split(ego_all, ego_all@result$ONTOLOGY)
by_ont_simpl <- lapply(by_ont, function(x)
  simplify(x, cutoff = 0.6, by = "p.adjust", select_fun = min))

# make three panels (one for each)
p_brain_BP <- make_enrich_panel("Brain — GO:BP", by_ont_simpl$BP,
                                show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(brain_up), universe_n=length(universe_entrez))
p_brain_CC <- make_enrich_panel("Brain — GO:CC", by_ont_simpl$CC, ...)
p_brain_MF <- make_enrich_panel("Brain — GO:MF", by_ont_simpl$MF, ...)

# Combine top terms across ontologies for a single panel:
df_all <- dplyr::bind_rows(lapply(by_ont_simpl, as.data.frame))
df_all <- df_all[order(df_all$p.adjust), ]        # rank across BP/CC/MF

# run GO for ALL, simplify by ontology
go_all_simplified <- function(entrez_vec, universe_entrez) {
  ego_all <- enrichGO(
    gene          = entrez_vec,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",       # BP + CC + MF
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego_all) || nrow(ego_all@result) == 0) return(list(BP=NULL, CC=NULL, MF=NULL))
  split_obj <- split(ego_all, ego_all@result$ONTOLOGY)
  lapply(split_obj, function(x) simplify(x, cutoff = 0.6, by = "p.adjust", select_fun = min))
}

dir.create(file.path(outdir, "enrichment_GO_panels"), showWarnings = FALSE, recursive = TRUE)

# run per tissue 
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# run one ontology, return an enrichResult (or NULL if empty)
.run_go <- function(genes, universe, ont) {
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,          # "BP", "CC", or "MF"
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  # simplify expects an enrichResult, not a data.frame
  simplify(ego, cutoff = 0.6, by = "p.adjust", select_fun = min)
}

# wrapper: run all three ontologies
run_go3 <- function(genes, universe) {
  list(
    BP = .run_go(genes, universe, "BP"),
    CC = .run_go(genes, universe, "CC"),
    MF = .run_go(genes, universe, "MF")
  )
}


go_brain    <- run_go3(genes_list$Brain,    universe_entrez)
go_liver    <- run_go3(genes_list$Liver,    universe_entrez)
go_pancreas <- run_go3(genes_list$Pancreas, universe_entrez)

p_brain_BP <- make_enrich_panel("Brain — GO:BP", go_brain$BP,
                                show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Brain), universe_n=length(universe_entrez))

# make panels (one per ontology) 
p_brain_BP <- make_enrich_panel("Brain — GO:BP", go_brain$BP,
                                show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Brain), universe_n=length(universe_entrez))
p_brain_CC <- make_enrich_panel("Brain — GO:CC", go_brain$CC,
                                show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Brain), universe_n=length(universe_entrez))
p_brain_MF <- make_enrich_panel("Brain — GO:MF", go_brain$MF,
                                show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Brain), universe_n=length(universe_entrez))

p_liver_BP <- make_enrich_panel("Liver — GO:BP", go_liver$BP,   show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Liver), universe_n=length(universe_entrez))
p_liver_CC <- make_enrich_panel("Liver — GO:CC", go_liver$CC,   show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Liver), universe_n=length(universe_entrez))
p_liver_MF <- make_enrich_panel("Liver — GO:MF", go_liver$MF,   show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Liver), universe_n=length(universe_entrez))

p_panc_BP  <- make_enrich_panel("Pancreas — GO:BP", go_pancreas$BP, show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Pancreas), universe_n=length(universe_entrez))
p_panc_CC  <- make_enrich_panel("Pancreas — GO:CC", go_pancreas$CC, show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Pancreas), universe_n=length(universe_entrez))
p_panc_MF  <- make_enrich_panel("Pancreas — GO:MF", go_pancreas$MF, show_n=5, alpha=alpha, lfc_min=lfc_min, cpm_min=cpm_min,
                                up_n=length(genes_list$Pancreas), universe_n=length(universe_entrez))

# save each
save_panel <- function(p, fname) {
  ggsave(file.path(outdir, "enrichment_GO_panels", paste0(fname, ".png")),
         p, width = 12, height = 5.5, dpi = 300,
         bg = rgb(230,218,209, maxColorValue = 255))
}

save_panel(p_brain_BP, "brain_GO_BP"); save_panel(p_brain_CC, "brain_GO_CC"); save_panel(p_brain_MF, "brain_GO_MF")
save_panel(p_liver_BP, "liver_GO_BP"); save_panel(p_liver_CC, "liver_GO_CC"); save_panel(p_liver_MF, "liver_GO_MF")
save_panel(p_panc_BP,  "pancreas_GO_BP"); save_panel(p_panc_CC, "pancreas_GO_CC"); save_panel(p_panc_MF, "pancreas_GO_MF")

utils::packageVersion("clusterProfiler")
exists("simplify", where = asNamespace("clusterProfiler"), inherits = FALSE)
getFromNamespace("simplify", "clusterProfiler")  # should return a function, not error

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler","DOSE","enrichplot","org.Hs.eg.db"))

library(clusterProfiler); library(org.Hs.eg.db)

ego_bp <- enrichGO(gene=genes_list$Brain, universe=universe_entrez,
                   OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                   ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
# IMPORTANT: ego_bp is an enrichResult; DO NOT coerce to data.frame before this
ego_bp_s <- clusterProfiler::simplify(ego_bp, cutoff=0.6, by="p.adjust", select_fun=min)

bg_col      <- rgb(255,255,255, maxColorValue = 255)
header_fill <- "#0C5A6B"; header_txt <- "white"

add_header <- function(p, title_text) {
  header_g <- ggplot() +
    geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=1),
              fill=header_fill, colour=header_fill) +
    annotate("text", x=.5, y=.5, label=title_text,
             colour=header_txt, size=4, fontface="bold") +
    theme_void()
  cowplot::plot_grid(header_g, p + theme(plot.margin=margin(6,6,6,6)),
                     ncol=1, rel_heights=c(0.10, 0.90))
}

save_two <- function(p, fname, w=10, h=6, dpi=300) {
  dir.create(file.path(outdir, "go_bp_plots"), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(outdir, "go_bp_plots", paste0(fname, ".png")),
         p, width=w, height=h, dpi=dpi, bg=bg_col)
  ggsave(file.path(outdir, "go_bp_plots", paste0(fname, ".pdf")),
         p, width=w, height=h, device="pdf", bg=bg_col)
}

# Compose the header line once
hdr <- sprintf("Brain — GO:BP • Selection: FDR<%g | |log2FC|≥%g | mean logCPM>%g | Up genes=%d | Universe=%d",
               alpha, lfc_min, cpm_min, length(genes_list$Brain), length(universe_entrez))

# DOTPLOT (top 15) 
stopifnot(!is.null(ego_bp_s), nrow(ego_bp_s@result) > 0)
p_dot <- dotplot(ego_bp_s, showCategory=15, font.size=11) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = bg_col, colour = NA)
  )
p_dot_h <- add_header(p_dot, hdr)
save_two(p_dot_h, "brain_GO_BP_dotplot", w=11, h=6)

# BARPLOT (top 15) 
p_bar <- barplot(ego_bp_s, showCategory=15) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = bg_col, colour = NA)
  )
p_bar_h <- add_header(p_bar, hdr)
save_two(p_bar_h, "brain_GO_BP_barplot", w=11, h=6)

# ENRICHMENT MAP (term–term graph) 
# needs pairwise term similarity; compute on the simplified object
ego_bp_sim <- tryCatch(enrichplot::pairwise_termsim(ego_bp_s), error=function(e) NULL)
if (!is.null(ego_bp_sim) && nrow(ego_bp_sim@result) > 1) {
  p_emap <- emapplot(ego_bp_sim, showCategory=30, cex_label_category=0.9, layout="kk") +
    theme(plot.background = element_rect(fill = bg_col, colour = NA))
  p_emap_h <- add_header(p_emap, hdr)
  save_two(p_emap_h, "brain_GO_BP_emap", w=11, h=8)
}

# CNET PLOT (term–gene network) 
# choose up to 10 categories for clarity
if (!is.null(ego_bp_sim) && nrow(ego_bp_sim@result) > 0) {
  p_cnet <- cnetplot(ego_bp_sim, showCategory=10, circular=FALSE, node_label="category") +
    theme(plot.background = element_rect(fill = bg_col, colour = NA))
  p_cnet_h <- add_header(p_cnet, hdr)
  save_two(p_cnet_h, "brain_GO_BP_cnet", w=11, h=8)
}

message("Saved GO:BP plots to: ", normalizePath(file.path(outdir, "go_bp_plots")))


# libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

# helpers: GO per ontology + redundancy reduction 
run_go_one <- function(genes, universe, ont=c("BP","CC","MF")) {
  ont <- match.arg(ont)
  er <- clusterProfiler::enrichGO(
    gene = genes, universe = universe,
    OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
    ont = ont, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE
  )
  if (is.null(er) || nrow(er@result)==0) return(NULL)
  er
}

# if clusterProfiler::simplify exists, use it; else do a light redundancy trim
simplify_safe <- function(er, cutoff = 0.6) {
  if (is.null(er) || nrow(er@result) < 2) return(er)
  if (exists("simplify", where = asNamespace("clusterProfiler"), inherits = FALSE)) {
    return(clusterProfiler::simplify(er, cutoff=cutoff, by="p.adjust", select_fun=min))
  }
  er2 <- enrichplot::pairwise_termsim(er)
  sim <- er2@termsim
  if (is.null(sim) || nrow(sim) < 2) return(er)
  d  <- as.dist(1 - sim); hc <- hclust(d, "average")
  cl <- cutree(hc, h = 1 - cutoff)
  res <- er2@result; res$.__cl__ <- cl[match(res$ID, rownames(sim))]
  keep <- res %>% group_by(.__cl__) %>% slice_min(p.adjust, n=1, with_ties = FALSE) %>% pull(ID)
  er_subset <- er; er_subset@result <- res[match(keep, res$ID), setdiff(names(res), ".__cl__"), drop=FALSE]
  er_subset
}

run_go3 <- function(genes, universe) {
  list(
    BP = simplify_safe(run_go_one(genes, universe, "BP")),
    CC = simplify_safe(run_go_one(genes, universe, "CC")),
    MF = simplify_safe(run_go_one(genes, universe, "MF"))
  )
}

# build small, mixed table (BP+CC+MF) per tissue
# take top_n terms overall, but enforce up to max_per_cat per category for variety
top_terms_df <- function(go_list, top_n = 6, max_per_cat = 3, wrap = 42) {
  dfs <- lapply(names(go_list), function(k){
    er <- go_list[[k]]; if (is.null(er) || nrow(er@result)==0) return(NULL)
    out <- as.data.frame(er@result)
    out$ONTOLOGY <- k
    out
  })
  df <- bind_rows(dfs)
  if (nrow(df)==0) return(tibble())
  df <- df %>% arrange(p.adjust) %>%
    group_by(ONTOLOGY) %>% slice_head(n = max_per_cat) %>% ungroup() %>%
    slice_head(n = top_n)
  tibble(
    `GO Term`          = df$ID,
    Description        = stringr::str_wrap(df$Description, width = wrap),
    Category           = recode(df$ONTOLOGY, BP="Biological", CC="Cellular", MF="Molecular"),
    `Adjusted p-value` = formatC(df$p.adjust, format="e", digits=2)
  )
}

#  table grob 
bg_col      <- rgb(230,218,209, maxColorValue = 255)
header_fill <- "#9C5A7a"; header_txt <- "white"

table_panel_from_df <- function(tissue, tab, alpha, lfc_min, cpm_min, up_n, universe_n) {
  if (nrow(tab)==0) {
    body <- ggplot() + theme_void() + annotate("text", x=.5,y=.5,
                                               label = sprintf("%s: no significant GO terms", tissue), size=5, fontface="bold")
    header <- ggplot() + theme_void()
    return(cowplot::plot_grid(header, body, ncol=1, rel_heights=c(0.0,1)))
  }
  
  tt <- ttheme_minimal(
    core = list(
      fg_params = list(hjust=0, x=0.02, fontsize=11),
      bg_params = list(fill = rep(c("grey95","white"), length.out = nrow(tab))),
      padding  = unit(c(5,8), "pt")
    ),
    colhead = list(
      fg_params = list(hjust=0, x=0.02, fontface="bold", col="white", fontsize=12),
      bg_params = list(fill = header_fill),
      padding  = unit(c(6,8), "pt")
    )
  )
  tg <- tableGrob(tab, rows = NULL, theme = tt)
  tg$widths <- unit(c(0.22, 0.52, 0.10, 0.16), "npc")
  
  hdr <- sprintf(
    "%s • Selection: FDR<%g  |  |log2FC|≥%g  |  mean logCPM>%g  |  Up genes=%s  |  Universe=%s",
    tissue, alpha, lfc_min, cpm_min, up_n, universe_n
  )
  header <- ggplot() +
    geom_rect(aes(xmin=0,xmax=1,ymin=0,ymax=1), fill=header_fill, colour=header_fill) +
    annotate("text", x=.5, y=.5, label=hdr, colour=header_txt, size=4, fontface="bold") +
    theme_void()
  
  panel <- cowplot::plot_grid(header, cowplot::ggdraw(tg), ncol=1, rel_heights=c(0.10,0.90))
  panel + theme(plot.background = element_rect(fill = bg_col, colour = NA))
}

#  RUN enrichment per tissue 
go_brain    <- run_go3(genes_list$Brain,    universe_entrez)
go_liver    <- run_go3(genes_list$Liver,    universe_entrez)
go_pancreas <- run_go3(genes_list$Pancreas, universe_entrez)

tab_brain    <- top_terms_df(go_brain,    top_n=6, max_per_cat=3, wrap=40)
tab_liver    <- top_terms_df(go_liver,    top_n=6, max_per_cat=3, wrap=40)
tab_pancreas <- top_terms_df(go_pancreas, top_n=6, max_per_cat=3, wrap=40)

p_brain_tbl <- table_panel_from_df("Brain",    tab_brain,    alpha, lfc_min, cpm_min,
                                   length(genes_list$Brain),    length(universe_entrez))
p_liver_tbl <- table_panel_from_df("Liver",    tab_liver,    alpha, lfc_min, cpm_min,
                                   length(genes_list$Liver),    length(universe_entrez))
p_panc_tbl  <- table_panel_from_df("Pancreas", tab_pancreas, alpha, lfc_min, cpm_min,
                                   length(genes_list$Pancreas), length(universe_entrez))

#SAVE slides like your examples 
dir.create(file.path(outdir, "enrichment_GO_tables"), showWarnings = FALSE, recursive = TRUE)

# rain & Liver side-by-side
ggsave(file.path(outdir, "enrichment_GO_tables", "GO_tables_brain.png"),
       p_brain_tbl, width = 16, height = 9, dpi = 300, bg = bg_col)

# Pancreas only
ggsave(file.path(outdir, "enrichment_GO_tables", "GO_tables_pancreas.png"),
       p_panc_tbl, width = 12, height = 7, dpi = 300, bg = bg_col)

ggsave(file.path(outdir, "enrichment_GO_tables", "GO_tables_liver.png"),
       p_liver_tbl, width = 12, height = 7, dpi = 300, bg = bg_col)

# GO Enrichment slides

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)   # change to your organism if needed
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  library(grid)
})

if (!exists("outdir")) outdir <- getwd()
if (!exists("alpha"))   alpha   <- 0.01
if (!exists("lfc_min")) lfc_min <- 1
if (!exists("cpm_min")) cpm_min <- 0

dir.create(file.path(outdir, "enrichment_GO_tables"), showWarnings = FALSE, recursive = TRUE)

bg_col      <- rgb(230,218,209, maxColorValue = 255)
header_fill <- "#9C5A7A"
header_txt  <- "white"

# Guess ID type
.guess_type <- function(x) {
  x <- unique(na.omit(as.character(x)))
  if (length(x) == 0) return("ENTREZID")
  all_digits <- suppressWarnings(all(!is.na(as.integer(x))))
  if (all_digits) "ENTREZID" else "SYMBOL"
}

# Ensure ENTREZ IDs (for both target gene set and universe)
to_entrez <- function(ids, fromType = NULL) {
  ids <- unique(na.omit(as.character(ids)))
  if (is.null(fromType)) fromType <- .guess_type(ids)
  if (fromType == "ENTREZID") return(ids)
  suppressMessages({
    map <- bitr(ids, fromType = fromType, toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  })
  unique(map$ENTREZID)
}

# Map any variant ("bp", "biological process", etc.) to the 3 GO names
.norm_onto <- function(x){
  x <- tolower(as.character(x))
  ifelse(grepl("^bp|biolog", x), "Biological Process",
         ifelse(grepl("^cc|cellul", x), "Cellular Component",
                ifelse(grepl("^mf|molecul", x), "Molecular Function", NA_character_)))
}

# Balanced picker across BP/CC/MF with category normalization
top_terms_df <- function(go_list, top_n = 6, max_per_cat = 3, wrap = 42) {
  as_res_df <- function(er){
    if (is.null(er)) return(NULL)
    if (inherits(er, "enrichResult")) as.data.frame(er@result) else as.data.frame(er)
  }
  dfs <- lapply(go_list, as_res_df)
  df  <- dplyr::bind_rows(dfs)
  if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
  
  if (!"ONTOLOGY" %in% names(df)) {
    reps <- unlist(lapply(dfs, nrow), use.names = FALSE)
    df$ONTOLOGY <- rep(names(go_list), reps)
  }
  df$Category <- .norm_onto(df$ONTOLOGY)
  df <- df %>% filter(!is.na(Category)) %>% arrange(p.adjust)
  
  per_cat <- df %>% group_by(Category) %>% slice_head(n = max_per_cat) %>% ungroup()
  
  K <- n_distinct(per_cat$Category)
  base_each <- max(1, floor(top_n / max(1, K)))
  
  even_part <- per_cat %>% group_by(Category) %>% slice_head(n = base_each) %>% ungroup()
  remainder <- anti_join(per_cat, even_part, by = c("ID","Category"))
  
  pick <- even_part
  if (nrow(pick) < top_n && nrow(remainder) > 0) {
    fill <- remainder %>% arrange(p.adjust) %>% slice_head(n = top_n - nrow(pick))
    pick <- bind_rows(pick, fill)
  }
  pick <- pick %>%
    arrange(factor(Category, levels = c("Biological Process","Cellular Component","Molecular Function")),
            p.adjust)
  
  tibble(
    `GO ID`            = pick$ID,
    Description        = str_wrap(pick$Description, width = wrap),
    Category           = pick$Category,
    `Adjusted p-value` = formatC(pick$p.adjust, format = "e", digits = 2)
  )
}

# Pretty table grob with header
table_panel_from_df <- function(tissue, tab, alpha, lfc_min, cpm_min, up_n, universe_n) {
  if (nrow(tab) == 0) {
    body <- ggplot() + theme_void() +
      annotate("text", x=.5, y=.5,
               label=sprintf("%s: no significant GO terms", tissue),
               size=5, fontface="bold")
    return(body + theme(plot.background = element_rect(fill = bg_col, colour = NA)))
  }
  
  tt <- ttheme_minimal(
    core = list(
      fg_params = list(hjust=0, x=0.02, fontsize=11),
      bg_params = list(fill = rep(c("grey95","white"), length.out = nrow(tab))),
      padding  = unit(c(5,8), "pt")
    ),
    colhead = list(
      fg_params = list(hjust=0, x=0.02, fontface="bold", col="white", fontsize=12),
      bg_params = list(fill = header_fill),
      padding  = unit(c(6,8), "pt")
    )
  )
  tg <- tableGrob(tab, rows = NULL, theme = tt)
  tg$widths <- unit(c(0.20, 0.52, 0.12, 0.16), "npc")
  
  hdr <- sprintf(
    "%s • Selection: FDR<%g  |  |log2FC|≥%g  |  mean logCPM>%g  |  Up genes=%s  |  Universe=%s",
    tissue, alpha, lfc_min, cpm_min, up_n, universe_n
  )
  header <- ggplot() +
    geom_rect(aes(xmin=0,xmax=1,ymin=0,ymax=1), fill=header_fill, colour=header_fill) +
    annotate("text", x=.5, y=.5, label=hdr, colour=header_txt, size=4, fontface="bold") +
    theme_void()
  
  panel <- cowplot::plot_grid(header, cowplot::ggdraw(tg), ncol=1, rel_heights=c(0.10,0.90))
  panel + theme(plot.background = element_rect(fill = bg_col, colour = NA))
}

# Run GO over BP/CC/MF for a gene set (ENTREZ IDs)
run_go3 <- function(entrez_ids, universe = NULL, p_cut = 0.05, q_cut = 0.05) {
  ontos <- c(BP="BP", CC="CC", MF="MF")
  out <- lapply(ontos, function(onto) {
    if (length(entrez_ids) == 0) return(NULL)
    enrichGO(
      gene          = entrez_ids,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = onto,
      universe      = universe,
      pvalueCutoff  = p_cut,
      qvalueCutoff  = q_cut,
      pAdjustMethod = "BH",
      readable      = TRUE
    )
  })
  out
}

# Load inputs if missing 
if (!exists("genes_list")) {
  message("genes_list not found in env — trying to read gene lists from outdir")
  f_b <- file.path(outdir, "brain_specific_up.txt")
  f_l <- file.path(outdir, "liver_specific_up.txt")
  f_p <- file.path(outdir, "pancreas_specific_up.txt")
  genes_list <- list(
    Brain    = if (file.exists(f_b)) readLines(f_b) else character(0),
    Liver    = if (file.exists(f_l)) readLines(f_l) else character(0),
    Pancreas = if (file.exists(f_p)) readLines(f_p) else character(0)
  )
}

# Convert input gene IDs to ENTREZ (auto-detect SYMBOL vs ENTREZ)
genes_list_entrez <- lapply(genes_list, function(v) to_entrez(v, fromType = NULL))

# Universe handling
if (exists("universe_entrez")) {
  universe_e <- unique(universe_entrez)
} else {
  # Fallback: use union of the three sets (conservative, but keeps script runnable)
  universe_e <- unique(unlist(genes_list_entrez, use.names = FALSE))
  message("No universe_entrez provided — using union of input genes as background.")
}

# Run enrichment per tissue
go_brain    <- run_go3(genes_list_entrez$Brain,    universe_e, p_cut = alpha, q_cut = alpha)
go_liver    <- run_go3(genes_list_entrez$Liver,    universe_e, p_cut = alpha, q_cut = alpha)
go_pancreas <- run_go3(genes_list_entrez$Pancreas, universe_e, p_cut = alpha, q_cut = alpha)

# Build top-term tables (balanced BP/CC/MF)
tab_brain    <- top_terms_df(go_brain,    top_n = 6, max_per_cat = 3, wrap = 40)
tab_liver    <- top_terms_df(go_liver,    top_n = 6, max_per_cat = 3, wrap = 40)
tab_pancreas <- top_terms_df(go_pancreas, top_n = 6, max_per_cat = 3, wrap = 40)

# Compose panels
p_brain_tbl <- table_panel_from_df("Brain",    tab_brain,    alpha, lfc_min, cpm_min,
                                   length(genes_list_entrez$Brain),    length(universe_e))
p_liver_tbl <- table_panel_from_df("Liver",    tab_liver,    alpha, lfc_min, cpm_min,
                                   length(genes_list_entrez$Liver),    length(universe_e))
p_panc_tbl  <- table_panel_from_df("Pancreas", tab_pancreas, alpha, lfc_min, cpm_min,
                                   length(genes_list_entrez$Pancreas), length(universe_e))

# Save slides
slide1_path <- file.path(outdir, "enrichment_GO_tables", "GO_tables_brain.png")
slide2_path <- file.path(outdir, "enrichment_GO_tables", "GO_tables_pancreas.png")
slide3_path <- file.path(outdir, "enrichment_GO_tables", "GO_tables_liver.png")


slide1 <- cowplot::plot_grid(p_brain_tbl, p_liver_tbl, ncol = 2, rel_widths = c(1,1))
ggsave(slide1_path, p_brain_tbl, width = 12, height = 7, dpi = 300, bg = bg_col)

ggsave(slide2_path, p_panc_tbl, width = 12, height = 7, dpi = 300, bg = bg_col)
ggsave(slide3_path, p_liver_tbl, width = 12, height = 7, dpi = 300, bg = bg_col)

message("Saved:\n  ", normalizePath(slide1_path), "\n  ", normalizePath(slide2_path))

