library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
objs <- load("/Users/elena.lippolis/Library/CloudStorage/OneDrive-Htechnopole/Attachments/Exams/T/SingleCell/SRA703206_SRS3296613.sparse.RData")
counts <- get(objs[1])
n_cells_before <- ncol(pbmc)
n_cells_before

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = counts, project = "mydata", min.cells = 3, min.features = 200)
pbmc
head(colnames(pbmc), 3)
grep("^MT-",rownames(pbmc),value = TRUE)
library(ggplot2)

# mitochondrial %
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
grep("^RP[LS]",rownames(pbmc),value = TRUE)
pbmc[["percent.rbp"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
head(pbmc@meta.data, 5)

p <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"),
  ncol = 4,
  cols = c("#4C78A8", "#F58518", "#54A24B"),
  pt.size = 0
)

# save as PNG (high resolution)
ggsave("./qc_violin_v2.png", plot = p, width = 10, height = 4, dpi = 300)

# you can also save as PDF
ggsave("./qc_violin_v2.pdf", plot = p, width = 10, height = 4)


library(patchwork)
# simple palette for the two plots
col1 <- "#4C78A8"  # cool blue
col2 <- "#4C78A8"  # warm orange

# a clean minimal theme
theme_scatter <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0),
    panel.background = element_rect(fill = "white", color = NA),   # white plotting area
    plot.background  = element_rect(fill = "white", color = NA)    # white around panel
  )

# build plots with nicer aesthetics
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.6) +
  labs(title = "nCount_RNA vs percent.mt", x = "nCount_RNA", y = "percent.mt") +
  theme_scatter +
  scale_color_manual(values = col1)

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.6) +
  labs(title = "nCount_RNA vs nFeature_RNA", x = "nCount_RNA", y = "nFeature_RNA") +
  theme_scatter +
  scale_color_manual(values = col2)

# side-by-side
plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.rbp", pt.size = 0.6) +
  labs(title = "nCount_RNA vs percent.rbp", x = "nCount_RNA", y = "percent.rbp") +
  theme_scatter +
  scale_color_manual(values = col2)


# save combined and individual plots
ggsave("./scatter_qc_plot3.png", plot3, width = 10, height = 4, dpi = 300)
ggsave("./scatter_qc_plot1.png", plot1, width = 5, height = 4, dpi = 300)
ggsave("./scatter_qc_plot2.png", plot2, width = 5, height = 4, dpi = 300)
plot3
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#thresholds: number of genes
pbmc
# drop the "-ENSG...(.version)" suffix, keep SYMBOL
old <- rownames(pbmc)
new <- sub("-ENSG\\d+(\\.\\d+)?$", "", old)

# sanity check
head(cbind(old, new), 5)
sum(duplicated(new))  # how many symbols would collide?

rownames(pbmc) <- make.unique(new)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc@assays

apply(pbmc[["RNA"]]$data,1,mean) -> gene.expression

sort(gene.expression, decreasing = TRUE) -> gene.expression

head(gene.expression, n=50)
VlnPlot(pbmc, features = c("B2M","MALAT1"))

cc.genes.updated.2019
CellCycleScoring(pbmc, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> pbmc

pbmc[[]]

#the default method -vst- computes (or better, estimates) the mean-variance relationship of each gene, and chooses the 2000 genes with hte highest variance. 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc@assays$RNA

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 50, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

ElbowPlot(pbmc, ndims=50)

pc.touse <- (pbmc$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse


pbmc <- FindNeighbors(pbmc, dims = 1:21)
pbmc <- FindClusters(pbmc, resolution = 0.5)


head(Idents(pbmc), 5)

DimPlot(pbmc, reduction = "pca")
DimPlot(pbmc,reduction="pca", dims=c(4,9))


pbmc <- RunTSNE(pbmc, dims=1:21)
DimPlot(pbmc, reduction = "tsne")

pbmc <- RunUMAP(pbmc, dims = 1:15)
#if you cannot install UMAP, t_SNE is anyway ok for your project!
DimPlot(pbmc, reduction = "umap")

VlnPlot(pbmc,features="nCount_RNA", pt.size = 0)
VlnPlot(pbmc,features="percent.mt", pt.size = 0)
VlnPlot(pbmc,features="percent.rbp", pt.size = 0)
# make sure identities are clusters
Idents(pbmc) <- "seurat_clusters"

# counts
table(Idents(pbmc))

# as a data frame with percents
library(dplyr)
cluster_sizes <- pbmc@meta.data %>%
  count(seurat_clusters, name = "n") %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  arrange(as.numeric(as.character(seurat_clusters)))  # tidy order
cluster_sizes
# -------
Idents(pbmc) <- "seurat_clusters"

pbmc_21 <- RunTSNE(pbmc, dims = 1:21)
levs21 <- levels(pbmc_21)
# count with levels preserved (no drops)
cts21  <- table(factor(Idents(pbmc_21), levels = levs21))
names(labs21) <- levs21
labs21 <- paste0(levs21, " (n=", as.integer(cts21), ")")
pbmc_21 <- RenameIdents(pbmc_21, labs21)
levels(pbmc_21) <- paste0(levels(pbmc_21), " (n=", cts21[levels(pbmc_21)], ")")

DimPlot(pbmc_21, reduction = "tsne") +
  ggtitle("t-SNE (dims 1:21)")

# t-SNE with 10 PCs
pbmc_10 <- RunTSNE(pbmc, dims = 1:15)

levs10 <- levels(pbmc_10)
cts10  <- table(factor(Idents(pbmc_10), levels = levs10))
labs10 <- paste0(levs10, " (n=", as.integer(cts10), ")")
names(labs10) <- levs10
pbmc_10 <- RenameIdents(pbmc_10, labs10)

DimPlot(pbmc_10, reduction = "tsne", repel = TRUE) +
  ggtitle("t-SNE (dims 1:15)")

library(Seurat)
library(ggplot2)
library(patchwork)

Idents(pbmc) <- "seurat_clusters"

feats <- c("nCount_RNA","nFeature_RNA","percent.mt","percent.rbp","S.Score","G2M.Score")

# make separate violins (so we can theme them), then combine
vlist <- VlnPlot(pbmc, features = feats, pt.size = 0, combine = FALSE)
vlist <- lapply(vlist, function(p) {
  p + theme_minimal(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
})

v_all <- wrap_plots(vlist, ncol = 3)
v_all

DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "seurat_clusters"

markers_1vAll <- FindAllMarkers(
  pbmc,
  only.pos        = TRUE,     # up in the cluster vs all others
  test.use        = "wilcox", # or "MAST" for UMI
  min.pct         = 0.25,
  logfc.threshold = 0.25
) %>%
  dplyr::filter(p_val_adj < 0.05)

# Top 10 per cluster (adjust n as you like)
top10 <- markers_1vAll %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

# Seurat’s default DoHeatmap palette is magenta–yellow-ish, similar to the vignette
hm <- DoHeatmap(pbmc, features = top10$gene, group.by = "seurat_clusters") + NoLegend()

#ggsave("/.../one_vs_all_heatmap_top10.png", hm, width = 10, height = 8, dpi = 300, bg = "white")
hm

marker_counts <- markers_1vAll %>%
  dplyr::count(cluster, name = "n_markers") %>%
  dplyr::arrange(n_markers)
marker_counts

# Jaccard overlap on top markers
topN <- 30
top_by_cl <- markers_1vAll %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = topN, with_ties = FALSE)

clu <- sort(unique(top_by_cl$cluster))
jac <- matrix(NA, length(clu), length(clu), dimnames = list(clu, clu))
for (i in seq_along(clu)) for (j in seq_along(clu)) {
  gi <- top_by_cl$gene[top_by_cl$cluster==clu[i]]
  gj <- top_by_cl$gene[top_by_cl$cluster==clu[j]]
  jac[i,j] <- length(intersect(gi,gj))/length(union(gi,gj))
}
jac  # high values (e.g., >0.3) -> do pairwise comparison

DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "seurat_clusters"

marker_panels <- list(
  Endothelial = c("PECAM1","CDH5","KDR","VWF","ENG","THBD","PLVAP","ACKR1", "LYVE1"),
  Myofibroblasts = c("MYH11", "ACTG2", "ACTA2"),
  Stromal = c("APOE", "CCL8", "FABP5", "ADAMDEC1", "DCN", "SLIT2", "SOX6", "F3/CD142", "WNT5A", "WNT5B", "BMP2", "BMP5", "FRZB", "POSTN", "HSD17B2", "SMAD7", "LMOD1", "TAGLN", "CAV1"),
  Fibroblast  = c("COL1A1","COL1A2","DCN","LUM","PDGFRA","COL3A1"),
  SmoothMuscle= c("ACTA2","MYH11","TAGLN","CNN1","RGS5"),
  Tcell       = c("CD3D","CD3E","TRAC","IL7R","CCR7","GZMB"),
  Bcell       = c("MS4A1","CD79A","CD79B","BANK1","CD74"),
  Plasma      = c("SDC1","XBP1","MZB1","JCHAIN","IGHG1","IGKC"),
  Myeloid     = c("LYZ","S100A8","S100A9","LGALS3","FCGR3A","LST1","MNDA"),
  Neutrophil  = c("S100A8","S100A9","FCGR3B","CXCR2","MPO"),
  Epithelial  = c("EPCAM","KRT8","KRT18","KRT19"),
  Adipocyte   = c("PLIN1","ADIPOQ","LIPE","FABP4"),
  Pericytes = c("RGS5", "AOC3"),
  Glial = c("S100B")
)
library(dplyr)
scores <- lapply(marker_panels, function(genes){
  genes <- intersect(genes, rownames(pbmc))
  AddModuleScore(pbmc, features = list(genes), name = paste0(genes[1],"_tmp"))[[paste0(genes[1],"_tmp1")]]
})
scores <- as.data.frame(scores)                                # cells × panels
colnames(scores) <- names(marker_panels)
scores$cluster <- Idents(pbmc)
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "seurat_clusters"

# keep only genes present
gs <- lapply(marker_panels, function(x) intersect(x, rownames(pbmc)))

# score (creates columns CTS_1 … CTS_k in meta.data)
pbmc <- AddModuleScore(pbmc, features = gs, name = "CTS_")

# rename score columns to panel names
score_cols <- grep("^CTS_\\d+$", colnames(pbmc@meta.data), value = TRUE)
names(score_cols) <- names(gs)
colnames(pbmc@meta.data)[match(score_cols, colnames(pbmc@meta.data))] <- paste0("CTS_", names(gs))

library(dplyr); library(tidyr)
scoredf <- pbmc@meta.data %>%
  select(seurat_clusters, starts_with("CTS_"))

by_cl <- scoredf %>%
  group_by(seurat_clusters) %>%
  summarise(across(starts_with("CTS_"), median, na.rm = TRUE), .groups = "drop")

assignments <- by_cl %>%
  pivot_longer(-seurat_clusters, names_to = "celltype", values_to = "score") %>%
  group_by(seurat_clusters) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

assignments
lab_map <- setNames(sub("^CTS_", "", assignments$celltype), assignments$seurat_clusters)
pbmc$celltype <- plyr::mapvalues(pbmc$seurat_clusters, from = names(lab_map), to = unname(lab_map))
DimPlot(pbmc, reduction = "tsne", group.by = "celltype", label = TRUE, repel = TRUE)

# pick a concise set of decisive markers to display per type
feat_show <- unique(unlist(lapply(marker_panels, head, 4)))
feat_show <- intersect(feat_show, rownames(pbmc))
feat_show

DotPlot(pbmc, features = feat_show, group.by = "celltype", dot.scale = 6) +
  RotatedAxis() + ggtitle("Canonical markers per assigned cell type")

DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "seurat_clusters"

# One-vs-all markers
markers <- FindAllMarkers(
  pbmc, only.pos = TRUE, min.pct = 0.25,
  logfc.threshold = 0.25, test.use = "wilcox"
)
markers <- dplyr::filter(markers, p_val_adj < 0.05)

# De-bias (drop mito/ribo/housekeeping)
is_mito  <- function(g) grepl("^MT-", g)
is_ribo  <- function(g) grepl("^RP[SL]", g)
is_house <- function(g) g %in% c("MALAT1","B2M","RPLP0","TMSB10","FTL","FTH1")

markers_clean <- markers |>
  dplyr::filter(!is_mito(gene), !is_ribo(gene), !is_house(gene)) |>
  dplyr::mutate(spec_score = avg_log2FC * pmax(pct.1 - pct.2, 0))

# Keep exactly 9 clusters (pick which ones!)
all_clust <- levels(pbmc)
#keep the 9 biggest clusters
keep9 <- names(sort(table(Idents(pbmc)), decreasing = TRUE))[1:10]


markers9 <- dplyr::filter(markers_clean, cluster %in% keep9)

# Top K “important features” per cluster
K <- 8  # change to 5/10 as you like
important_feats_tbl <- markers9 |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = spec_score, n = K, with_ties = FALSE) |>
  dplyr::arrange(cluster, dplyr::desc(spec_score)) |>
  dplyr::ungroup()

# Save a CSV
#write.csv(important_feats_tbl,
          "/path/to/important_features_topK_per_cluster.csv",
          row.names = FALSE)

# produce a compact list (nice to paste in a report)
feat_list <- split(important_feats_tbl$gene, important_feats_tbl$cluster)
feat_list


DotPlot(pbmc, features = feat_show, group.by = "celltype", dot.scale = 6) +
  RotatedAxis() + ggtitle("Canonical markers per assigned cell type")

DotPlot(pbmc, features = c("RGCC","MEG3","MCTP1",
                           "APOE","NRXN1","NOTCH3","ACKR1",
                           "RGS5","TYROBP", "TNFRSF17")
)

# FeaturePlots for key markers
FeaturePlot(pbmc, features = c("RGCC","MEG3","MCTP1",
                               "APOE","NRXN1","NOTCH3","ACKR1",
                               "RGS5","TYROBP", "TNFRSF17"), reduction="tsne",ncol=4)


library(ggplot2)
pbmc@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")


cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
head(cluster2.markers, n = 5)

cluster2_01.markers <- FindMarkers(pbmc, ident.1 = 2, ident.2 = c(0, 1), min.pct = 0.25)
head(cluster2_01.markers, n = 5)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)



FeaturePlot(pbmc, features = cfeatures_raw <- c("RGCC", "RGS17", "LUM","CADM2","ENPEP","STXBP6","COX4I2","CCL4L2","IGLL5")
)


pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster0.markers <- cluster0.markers[order(-cluster0.markers$avg_log2FC),]
head(cluster0.markers, n = 10)
# RGCC TM4SF18 CYYR1 PLVAP FABP5 FAM167B

cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25, test.use = "wilcox")
cluster1.markers <- cluster1.markers[order(-cluster1.markers$avg_log2FC),]
head(cluster1.markers, n = 10)
# HELZ2

cluster10.markers <- FindMarkers(pbmc, ident.1 = 2, ident.2 = 0, min.pct = 0.25, test.use = "wilcox")
cluster10.markers <- cluster10.markers[order(-cluster10.markers$avg_log2FC),]
head(cluster10.markers, n = 10)


cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
cluster2.markers <- cluster2.markers[order(-cluster2.markers$avg_log2FC),]
head(cluster2.markers, n = 10)
# HECW2 MCTP1

cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25, test.use = "wilcox")
cluster3.markers <- cluster3.markers[order(-cluster3.markers$avg_log2FC),]
head(cluster3.markers, n = 20)
 
cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25, test.use = "wilcox")
cluster4.markers <- cluster4.markers[order(-cluster4.markers$avg_log2FC),]
head(cluster4.markers, n = 10)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25, test.use = "wilcox")
cluster5.markers <- cluster5.markers[order(-cluster5.markers$avg_log2FC),]
head(cluster5.markers, n = 10)

cluster6.markers <- FindMarkers(pbmc, ident.1 = 6, min.pct = 0.25, test.use = "wilcox")
cluster6.markers <- cluster6.markers[order(-cluster6.markers$avg_log2FC),]
head(cluster6.markers, n = 10)

cluster7.markers <- FindMarkers(pbmc, ident.1 = 7, min.pct = 0.25, test.use = "wilcox")
cluster7.markers <- cluster7.markers[order(-cluster7.markers$avg_log2FC),]
head(cluster7.markers, n = 10)

cluster9.markers <- FindMarkers(pbmc, ident.1 = 9, min.pct = 0.25, test.use = "wilcox")
cluster9.markers <- cluster9.markers[order(-cluster9.markers$avg_log2FC),]
head(cluster9.markers, n = 10)

VlnPlot(pbmc, features = c("IGLL5", "MCTP1"))

DotPlot(pbmc, features = c("RGCC","MEG3","THBD",
                           "LUM","NRXN1","RGS5","ACKR1","NDUFA4L2","CD79A")
)


new.cluster.ids <- c("Endothelial cells", "Fibroblasts", "Adipocytes", "Macrophages", "Glia cells", "Smooth muscle cells",
                     "Neutrophils", "Fibroblasts", "Monocytes", "Plasma Cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

library(dplyr)

Idents(pbmc) <- "seurat_clusters"

qc_by_cluster <- pbmc@meta.data %>%
  mutate(cluster = Idents(pbmc)) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    med_nCount = median(nCount_RNA),
    med_nFeature = median(nFeature_RNA),
    med_mt = median(percent.mt, na.rm=TRUE),
    med_rbp = median(percent.rbp %||% NA, na.rm=TRUE)  # if you computed it
  ) %>%
  ungroup()

# Z-score each metric across clusters to spot outliers
z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
qc_by_cluster <- qc_by_cluster %>%
  mutate(
    z_nCount = z(med_nCount),
    z_mt     = z(med_mt),
    z_nFeat  = z(med_nFeature),
    z_rbp    = if (!all(is.na(med_rbp))) z(med_rbp) else NA_real_
  )

qc_by_cluster

library(stringr)

markers <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 0.25,
                          min.pct = 0.1, return.thresh = 0.05)

is_mito  <- function(g) str_detect(g, "^MT-")
is_ribo  <- function(g) str_detect(g, "^RP[SL]")   # RPS/RPL
is_house <- function(g) g %in% c("MALAT1","FTL","FTH1","B2M","RPLP0","TMSB10")

bias_summary <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30, with_ties = FALSE) %>%  # look at top 30
  summarise(
    n_top = n(),
    frac_mito = mean(is_mito(gene)),
    frac_ribo = mean(is_ribo(gene)),
    frac_house= mean(is_house(gene))
  ) %>%
  mutate(flag_marker_bias = (frac_mito > 0.2) | (frac_ribo > 0.3) | (frac_house > 0.3))

bias_summary

DefaultAssay(pbmc) <- "RNA"

# genes you want to show
features_raw <- c("PLVAP","MEG3","THBD","LGALS1","NRXN1","RGS5","ACKR1","NDUFA4L2","CD79A")


# map plain symbols to rownames (handles '.1' etc. from make.unique)
map_features <- function(obj, feats) {
  rn   <- rownames(obj)
  base <- sub("\\.\\d+$", "", rn)           # strip .1/.2 suffixes
  mapped <- vapply(feats, function(g) {
    hits <- rn[base == g]
    if (length(hits)) hits[1] else NA_character_
  }, character(1))
  data.frame(requested = feats, mapped = mapped, present = !is.na(mapped), stringsAsFactors = FALSE)
}

m <- map_features(pbmc, features_raw)
if (any(!m$present)) {
  message("Not found (after symbol mapping): ", paste(m$requested[!m$present], collapse = ", "))
}
features_mapped <- m$mapped[m$present]
# Keep original order of requested features
features_mapped <- features_mapped[match(m$requested[m$present], m$requested[m$present])]

# Make a readable DotPlot:
dp <- DotPlot(
  pbmc,
  features  = features_mapped,
  dot.scale = 6,                  # larger dots
  cols      = c("grey90", "#3B5BFF"),  # low->high color
  col.min   = 0,                  # clamp low end at zero
  col.max   = 2,                  # increase dynamic range; adjust 1.5–2.5 if needed
  scale.by  = "size",             # size = % expressed (legend "Percent Expressed")
  scale.min = 0,                  # don't rescale colors across features
  scale.max = 0                   # (keeps absolute color scale with col.min/col.max)
) + RotatedAxis() +
  labs(x = "Features", y = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(".c/dotplot_markers_fixed.png",
       dp, width = 10, height = 6, dpi = 300, bg = "white")
dp

DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "seurat_clusters"

# Find markers
mk <- FindAllMarkers(pbmc, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, test.use="wilcox") %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(
    is_mito  = grepl("^MT-", gene),
    is_ribo  = grepl("^RP[SL]", gene),
    is_house = gene %in% c("MALAT1","B2M","RPLP0","TMSB10","FTL","FTH1"),
    delta_pct = pmax(pct.1 - pct.2, 0),
    score = avg_log2FC * delta_pct
  ) %>%
  dplyr::filter(!is_mito, !is_ribo, !is_house)

# Pick top 3 per cluster by specificity score
top_markers <- mk %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = score, n = 3) %>% dplyr::ungroup()

# Map gene symbols to actual rownames (handles ".1" etc.)
map_features <- function(obj, feats) {
  rn <- rownames(obj); base <- sub("\\.\\d+$", "", rn)
  vapply(feats, function(g){ hits <- rn[base==g]; if(length(hits)) hits[1] else NA_character_ }, "")
}
feat <- unique(stats::na.omit(map_features(pbmc, top_markers$gene)))

# Plot
DotPlot(pbmc, features = feat, dot.scale = 6, cols = c("grey90", "#3B5BFF"),
        col.min = 0, col.max = 2, scale.by = "size") + RotatedAxis()

library(dplyr)
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "seurat_clusters"

markers_1vAll <- FindAllMarkers(
  pbmc,
  only.pos       = TRUE,     # only genes up in the cluster vs rest
  min.pct        = 0.25,
  logfc.threshold= 0.25,
  test.use       = "wilcox"  # you can switch to "MAST" for UMI data if desired
) %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

# Quick bias flags (mito/ribo/housekeeping)
is_mito  <- function(g) grepl("^MT-", g)
is_ribo  <- function(g) grepl("^RP[SL]", g)
is_house <- function(g) g %in% c("MALAT1","B2M","RPLP0","TMSB10","FTL","FTH1")

markers_1vAll <- markers_1vAll %>%
  mutate(
    mito  = is_mito(gene),
    ribo  = is_ribo(gene),
    house = is_house(gene)
  )

# Save full table
write.csv(markers_1vAll,
          "./markers_one_vs_all.csv",
          row.names = FALSE)

# Top 5 per cluster, de-biased (drop mito/ribo/house if you want the ‘clean’ list)
top5_clean <- markers_1vAll %>%
  filter(!mito, !ribo, !house, p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(top5_clean,
          "./top5_markers_clean_one_vs_all.csv",
          row.names = FALSE)

top5_clean

library(dplyr)
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "seurat_clusters"

markers_1vAll <- FindAllMarkers(
  pbmc,
  only.pos       = TRUE,     # only genes up in the cluster vs rest
  min.pct        = 0.25,
  logfc.threshold= 0.25,
  test.use       = "wilcox"  # you can switch to "MAST" for UMI data if desired
) %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

# Quick bias flags (mito/ribo/housekeeping)
is_mito  <- function(g) grepl("^MT-", g)
is_ribo  <- function(g) grepl("^RP[SL]", g)
is_house <- function(g) g %in% c("MALAT1","B2M","RPLP0","TMSB10","FTL","FTH1")

markers_1vAll <- markers_1vAll %>%
  mutate(
    mito  = is_mito(gene),
    ribo  = is_ribo(gene),
    house = is_house(gene)
  )

# Save full table
write.csv(markers_1vAll,
          "./markers_one_vs_all.csv",
          row.names = FALSE)

# Top 5 per cluster, de-biased (drop mito/ribo/house if you want the ‘clean’ list)
top5_clean <- markers_1vAll %>%
  filter(!mito, !ribo, !house, p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(top5_clean,
          "./top5_markers_clean_one_vs_all.csv",
          row.names = FALSE)

levs <- levels(pbmc)
pair_list <- combn(levs, 2, simplify = FALSE)

run_pair <- function(a, b, direction = c("A_up","B_up")) {
  direction <- match.arg(direction)
  if (direction == "A_up") {
    res <- FindMarkers(pbmc, ident.1 = a, ident.2 = b, min.pct = 0.25,
                       logfc.threshold = 0.25, test.use = "wilcox")
    res$contrast <- paste0(a, "_vs_", b)
    res$up_in    <- a
  } else {
    res <- FindMarkers(pbmc, ident.1 = b, ident.2 = a, min.pct = 0.25,
                       logfc.threshold = 0.25, test.use = "wilcox")
    res$contrast <- paste0(b, "_vs_", a)
    res$up_in    <- b
  }
  res$gene <- rownames(res)
  res
}

pairwise_markers <- do.call(
  rbind,
  lapply(pair_list, function(p) {
    a <- p[1]; b <- p[2]
    rbind(
      run_pair(a,b,"A_up"),
      run_pair(a,b,"B_up")
    )
  })
) %>%
  relocate(gene, contrast, up_in) %>%
  arrange(contrast, desc(avg_log2FC), p_val_adj)

# Add bias flags and save
pairwise_markers <- pairwise_markers %>%
  mutate(
    mito  = is_mito(gene),
    ribo  = is_ribo(gene),
    house = is_house(gene)
  )

write.csv(pairwise_markers,
          "./markers_pairwise_all.csv",
          row.names = FALSE)

pairwise_markers

# Heatmap of top 10 per cluster
top10 <- markers_1vAll %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# DotPlot for selected genes
DotPlot(pbmc, features = c("PLVAP","MEG3","THBD","NRXN1","CD79A")) + RotatedAxis()

p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("UMAP — Seurat clusters")
print(p_umap)

# -------
library(Seurat)
objs <- load("/SRA703206_SRS3296613.sparse.RData")
counts <- get(objs[1])
n_cells_before <- ncol(pbmc)
n_cells_before
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = counts, project = "mydata", min.cells = 3, min.features = 200)
pbmc
head(colnames(pbmc), 3)
grep("^MT-",rownames(pbmc),value = TRUE)
library(ggplot2)

# mitochondrial %
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
grep("^RP[LS]",rownames(pbmc),value = TRUE)
pbmc[["percent.rbp"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
head(pbmc@meta.data, 5)

DefaultAssay(pbmc) <- "RNA"

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                      scale.factor = 1e4, verbose = FALSE)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Be explicit about assay (and layer if you have it) to avoid the warning:
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc),
                  assay = "RNA", verbose = FALSE)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 50, verbose = FALSE)
pbmc

# barplots of top loading genes for PC1 and PC2
VizDimLoadings(pbmc, dims = 1, nfeatures = 30, balanced = TRUE)
VizDimLoadings(pbmc, dims = 2, nfeatures = 30, balanced = TRUE)
