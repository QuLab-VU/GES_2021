library(Seurat) # version 2.3.4 to create file hashed dataframes; 3.x to run analysis
library(Matrix)
library(AnnotationDbi)
library(DoubletFinder)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(NbClust)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VISION)
library(limma)
library(stringr)
library(data.table)
library("org.Hs.eg.db")
setwd('~/git/GES_2020/scRNAseq')

# Read in hashing data
matrix_dir = "umi_count/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
dat_HTO <- mat

# Read in scRNAseq data
matrix_dir_rna = "read_count/"
barcode.path_rna <- paste0(matrix_dir_rna, "barcodes.tsv.gz")
features.path_rna <- paste0(matrix_dir_rna, "features.tsv.gz")
matrix.path_rna <- paste0(matrix_dir_rna, "matrix.mtx.gz")
mat_rna <- readMM(file = matrix.path_rna)
feature.names_rna = read.delim(features.path_rna, header = FALSE, stringsAsFactors = FALSE)
barcode.names_rna = read.delim(barcode.path_rna, header = FALSE, stringsAsFactors = FALSE)
barcodes_raw = barcode.names_rna$V1
colnames(mat_rna) = sapply(barcodes_raw, function(x) strsplit(x, "-")[[1]][1], 
                           USE.NAMES = F)
rownames(mat_rna) = feature.names_rna$V1
dat_expr <- mat_rna

# Add rownames to RNA dataframe
d <- as.data.frame(as.matrix(dat_expr))
d$SYMBOL <- mapIds(org.Hs.eg.db,keys=rownames(dat_expr),
                   column="SYMBOL",keytype="ENSEMBL",multiVals="first")

rownames(dat_expr) <- make.names(ifelse(is.na(d$SYMBOL), 
                                        rownames(dat_expr), 
                                        d$SYMBOL),
                                 unique = TRUE)

# Keep barcodes shared between dataframes
joint_bcs <- intersect(colnames(dat_expr),colnames(dat_HTO))
dat_expr <- dat_expr[,joint_bcs]
dat_HTO <- as.matrix(dat_HTO[,joint_bcs])

# save(dat_expr, dat_HTO, file = "/Volumes/Flash/preSeurat_DD.RData")
# load('/Volumes/Flash/preSeurat_DD.RData')

# Setup Seurat object
dat_hashtag <- CreateSeuratObject(counts = dat_expr)

# Normalize RNA data with log normalization
dat_hashtag <- NormalizeData(dat_hashtag)
# Find and scale variable features
dat_hashtag <- FindVariableFeatures(dat_hashtag, selection.method = "mean.var.plot")
dat_hashtag <- ScaleData(dat_hashtag, features = VariableFeatures(dat_hashtag))
dat_hashtag <- CellCycleScoring(dat_hashtag, s.features = cc.genes.updated.2019$s.genes, 
                                g2m.features = cc.genes.updated.2019$g2m.genes, 
                                set.ident = FALSE)

# Add HTO data as a new assay independent from RNA
dat_hashtag[["HTO"]] <- CreateAssayObject(counts = dat_HTO)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
dat_hashtag <- NormalizeData(dat_hashtag, assay = "HTO", normalization.method = "CLR")

# Demultiplex HTO data using default settings
dat_hashtag <- HTODemux(dat_hashtag, assay = "HTO", positive.quantile = 0.99)

# # Save demultiplexed matrix as a CSV file (for Reviewer #1)
# data_to_write_out <- as.data.frame(as.matrix(dat_hashtag@assays$RNA@counts))
# fwrite(x = data_to_write_out, row.names = TRUE, col.names = TRUE, file = "PC9_scRNAseqCounts_HTOdemux.csv")

# Group cells based on the max HTO signal
Idents(dat_hashtag) <- "HTO_maxID"
RidgePlot(dat_hashtag, assay = "HTO", features = rownames(dat_hashtag[["HTO"]]),
          ncol = 3, cols = c("green", "seagreen", "darkorchid", "red", "blue", 
                             "gold", "brown", "deeppink")) + 
  ggsave("FIG_S15.pdf", width = 15, height = 15)


# Run dimensionality reduction on all cells
dat_hashtag <- RunPCA(dat_hashtag, features = VariableFeatures(dat_hashtag))
dat_hashtag <- FindNeighbors(dat_hashtag, reduction = "pca", dims = 1:10)
dat_hashtag <- FindClusters(dat_hashtag, resolution = 0.6, verbose = FALSE)
dat_hashtag <- RunTSNE(dat_hashtag, reduction = "pca", dims = 1:10)
dat_hashtag <- RunUMAP(dat_hashtag, dims = 1:10)
# Projecting singlet identities on UMAP visualization
DimPlot(dat_hashtag, reduction = "umap", group.by = "HTO_classification") 

# Classify cells by primary HTO
Idents(dat_hashtag) <- "HTO_classification.global"

# Identify predicted doublets using doubletFinder
## Not remove until later - just for plotting
sweep.res.list_dat_hashtag <- paramSweep_v3(dat_hashtag, PCs = 1:10)
sweep.stats_dat_hashtag <- summarizeSweep(sweep.res.list_dat_hashtag, GT = FALSE)
bcmvn_dat_hashtag <- find.pK(sweep.stats_dat_hashtag)

homotypic.prop_dat_hashtag <- modelHomotypic(annotations = dat_hashtag@meta.data$seurat_clusters)
nExp_poi_dat_hashtag <- round(0.05*nrow(dat_hashtag@meta.data))
nExp_poi.adj <- round(nExp_poi_dat_hashtag*(1-homotypic.prop_dat_hashtag))

dat_hashtag <- doubletFinder_v3(dat_hashtag, PCs = 1:10, pN = 0.25,
                                pK = 0.09, nExp = nExp_poi_dat_hashtag, 
                                reuse.pANN = FALSE)

found_dH <- match(dat_hashtag@meta.data$HTO_classification,
                  old_tags)
met_dH <- data.frame(ifelse(is.na(found_dH), dat_hashtag@meta.data$HTO_classification,
                            new_tags[found_dH]), row.names = rownames(dat_hashtag@meta.data))
names(met_dH) <- "Population"
dat_hashtag <- AddMetaData(dat_hashtag, met_dH, 'Population')
dat_SDP <- FetchData(dat_hashtag, vars = c("UMAP_1", "UMAP_2", "nCount_RNA", "nFeature_RNA",
                                           "HTO_classification.global",
                                           "DF.classifications_0.25_0.09_737",
                                           "Population"))

dat_SDP$HTO_type <- paste0("HTO", "-", dat_SDP$HTO_classification.global)
dat_SDP$DD_type <- paste0("DD", "-", dat_SDP$DF.classifications_0.25_0.09_737)
dat_SDP$Quality <- ifelse(dat_SDP$nCount_RNA > 15000 & dat_SDP$nFeature_RNA > 3000, 
                          "Good", "Poor")

dat_SDP_singlet <- subset(dat_SDP, Population %in% c("VU", "MGH", "BR1", "DS3",
                                                     "DS6", "DS7", "DS8", "DS9"))

# Plot doublets (both types) and poor quality cells with singlets
## HTO singlet/doublet
plt_HTO <- ggplot() + theme_bw() +
  geom_density_2d(data = dat_SDP_singlet,
                  aes(x=UMAP_1, y=UMAP_2,
                      color = Population), 
                  alpha = 0.8) +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid",
                                "seagreen", "gold", "green", "blue")) +
  geom_point(data = dat_SDP, shape = 21, alpha = 0.7,
             aes(x = UMAP_1, y = UMAP_2, fill = HTO_type),
             size = 0.8, stroke = 0.01) + 
  scale_fill_manual(values = c("grey10", "grey90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12))

plt_HTO + ggsave("FIG_S16A.svg", width = 4, height = 3)
plt_HTO_leg <- ggpubr::get_legend(plt_HTO)
as_ggplot(plt_HTO_leg) + ggsave("FIG_S16A_legend.svg",
                                width = 2.5, height = 4)

## DD singlet/doublet
plt_DD <- ggplot() + theme_bw() +
  geom_density_2d(data = dat_SDP_singlet,
                  aes(x=UMAP_1, y=UMAP_2,
                      color = Population), 
                  alpha = 0.8) +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid",
                                "seagreen", "gold", "green", "blue")) +
  geom_point(data = dat_SDP, shape = 21, alpha = 0.7,
             aes(x = UMAP_1, y = UMAP_2, fill = DD_type),
             size = 0.8, stroke = 0.01) + 
  scale_fill_manual(values = c("grey10", "grey90")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12))

plt_DD + ggsave("FIG_S16B.svg", width = 4, height = 3)
plt_DD_leg <- ggpubr::get_legend(plt_DD)
as_ggplot(plt_DD_leg) + ggsave("FIG_S16B_legend.svg",
                               width = 2.5, height = 4)

## Cutoff good/bad cells
plt_qual <- ggplot() + theme_bw() +
  geom_density_2d(data = dat_SDP_singlet,
                  aes(x=UMAP_1, y=UMAP_2,
                      color = Population), 
                  alpha = 0.8) +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid",
                                "seagreen", "gold", "green", "blue")) +
  geom_point(data = dat_SDP, shape = 21, alpha = 0.7,
             aes(x = UMAP_1, y = UMAP_2, fill = Quality),
             size = 0.8, stroke = 0.01) + 
  scale_fill_manual(values = c("grey90", "grey10")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12))

plt_qual + ggsave("FIG_S16C.svg", width = 4, height = 3)
plt_qual_leg <- ggpubr::get_legend(plt_qual)
as_ggplot(plt_qual_leg) + ggsave("FIG_S16C_legend.svg",
                                 width = 2.5, height = 4)

# Extract the singlets (predicted by Seurat)
dat_singlet <- subset(dat_hashtag, idents = "Singlet")
# Select the most variable features
dat_singlet <- FindVariableFeatures(dat_singlet, selection.method = "mean.var.plot")
VariableFeaturePlot(dat_singlet) + ggsave("FIG_S5B.svg", 
                                          width = 7, height = 4)
# Scaling RNA data, only scale the variable features here for efficiency
dat_singlet <- ScaleData(dat_singlet, features = VariableFeatures(dat_singlet))
# Run dimensionality reduction on singlets only
dat_singlet <- RunPCA(dat_singlet, features = VariableFeatures(dat_singlet))
dat_singlet <- FindNeighbors(dat_singlet, reduction = "pca", dims = 1:10)
dat_singlet <- FindClusters(dat_singlet, resolution = 0.6, verbose = FALSE)
dat_singlet <- RunTSNE(dat_singlet, reduction = "pca", dims = 1:10)
dat_singlet <- RunUMAP(dat_singlet, dims = 1:10)
# Projecting singlet identities on UMAP visualization
DimPlot(dat_singlet, reduction = "umap", group.by = "HTO_classification") 

# Identify predicted doublets - Doublet Finder
sweep.res.list_dat_singlet <- paramSweep_v3(dat_singlet, PCs = 1:10)
sweep.stats_dat_singlet <- summarizeSweep(sweep.res.list_dat_singlet, GT = FALSE)
bcmvn_dat_singlet <- find.pK(sweep.stats_dat_singlet)

homotypic.prop <- modelHomotypic(annotations = dat_singlet@meta.data$seurat_clusters)
nExp_poi <- round(0.05*nrow(dat_singlet@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

dat_singlet <- doubletFinder_v3(dat_singlet, PCs = 1:10, pN = 0.25,
                                pK = 0.09, nExp = nExp_poi, 
                                reuse.pANN = FALSE)

old_tags <- c("HTO-A-GTCAACTCTTTAGCG",
              "HTO-B-TGATGGCCTATTGGG",
              "HTO-C-TTCCGCCTCTCTTTG",
              "HTO-D-AGTAAGTTCAGCGTA",
              "HTO-E-AAGTATCGTTTCGCA",
              "HTO-F-GGTTGCCAGATGTCA",
              "HTO-G-TGTCTTTCCTGCCAG",
              "HTO-H-CTCCTCTGCAATTAC")
new_tags <- c("VU", "MGH", "BR1", "DS3",
              "DS6", "DS7", "DS8", "DS9")
found <- match(dat_singlet@meta.data$HTO_classification,
               old_tags)
met <- data.frame(ifelse(is.na(found), dat_singlet@meta.data$HTO_classification,
                         new_tags[found]), row.names = rownames(dat_singlet@meta.data))
names(met) <- "Population"
dat_singlet <- AddMetaData(dat_singlet, met, 'Population')

# Remove doublets (predicted and multiple cell hashtags)
dat_singlet_doubRem <- subset(dat_singlet, DF.classifications_0.25_0.09_490 == "Singlet" &
                                HTO_classification.global == "Singlet")

# # Identify cell quality
# test <- FetchData(dat_singlet_doubRem, 
#                   vars = c("nCount_RNA", "nFeature_RNA", "Population", "seurat_clusters"))
# 
# test %>%
#   group_by(seurat_clusters) %>%
#   dplyr::summarize(Mean = mean(nCount_RNA, na.rm=TRUE))

# Remove low quality cells
## Weak - used for inferCNV detection - maximize number of cells
dat_singlet_doubRem_poorRem <- subset(dat_singlet_doubRem,
                                      nCount_RNA > 10000 &
                                        nFeature_RNA > 2000)

# save(dat_singlet_doubRem_poorRem, file = "inferCNV/Seurat_v3_scRNAseq_forInferCNV.RData")

## Strong - used for plotting and clustering
dat_singlet_doubRem_poorRem_ex <- subset(dat_singlet_doubRem,
                                         nCount_RNA > 15000 &
                                           nFeature_RNA > 3000)

# save(dat_singlet_doubRem_poorRem_ex, dat_hashtag,
#      file = "Seurat_v3_scRNAseq_Paper.RData")
# load("Seurat_v3_scRNAseq_Paper.RData")

# Plots
# UMAP - CLV, CLV density + sublines dots
## UMAP - by Population
plt_all_UMAP <- DimPlot(dat_singlet_doubRem_poorRem_ex, reduction = "umap", group.by = "Population") +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid",
                                "seagreen", "gold", "green", "blue")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12), legend.text = element_text(size=12),
        axis.title = element_text(size=12), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## All samples - not in paper
plt_all_UMAP + ggsave("UMAP_PopColored.svg", width = 4, height = 3)
plt_all_UMAP_leg <- ggpubr::get_legend(plt_all_UMAP)
as_ggplot(plt_all_UMAP_leg) + ggsave("UMAP_PopColored_legend.svg",
                                     width = 2.5, height = 4)


## UMAP - by Population - only CLV
dat_CLV <- subset(dat_singlet_doubRem_poorRem_ex, Population %in% c("VU", "MGH", "BR1"))
CLV_plotDat <- FetchData(dat_CLV, vars = c("UMAP_1", "UMAP_2", "Population"))

plt_CLV <- ggplot(CLV_plotDat, 
                  aes(x = UMAP_1, y = UMAP_2,
                      color = Population)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12), legend.text = element_text(size=12),
        axis.title = element_text(size=12), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plt_CLV + ggsave("FIG_3F.svg", width = 4, height = 3)
plt_CLV_leg <- ggpubr::get_legend(plt_CLV)
as_ggplot(plt_CLV_leg) + ggsave("UMAP_PopColored_CLVonly_legend.svg",
                                width = 2.5, height = 4)

## UMAP - by Population - sublines with CLV underlay
dat_sublines <- subset(dat_singlet_doubRem_poorRem_ex, Population %in% c("DS3", "DS6", "DS7", "DS8", "DS9"))
sublines_plotDat <- FetchData(dat_sublines, vars = c("UMAP_1", "UMAP_2", "Population"))

plt_sublines_CLVoverlay <- ggplot() + theme_bw() +
  geom_density_2d(data = CLV_plotDat, 
                  aes(x = UMAP_1, y = UMAP_2,
                      color = Population), alpha = 0.8) +
  scale_color_manual(values = c("red", "green", "blue")) +
  geom_point(data = sublines_plotDat, shape = 21,
             aes(x = UMAP_1, y = UMAP_2, fill = Population),
             size = 0.8, stroke = 0.01) + 
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid",
                               "seagreen", "gold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12)) 

plt_sublines_CLVoverlay + ggsave("FIG_4F.svg", 
                                 width = 4, height = 3)
plt_sublines_CLVoverlay_leg <- ggpubr::get_legend(plt_sublines_CLVoverlay)
as_ggplot(plt_sublines_CLVoverlay_leg) + ggsave("UMAP_PopColored_sublines_CLVoverlay_legend.svg", 
                                                width = 2.5, height = 4)

## UMAP - by cell cycle phase
dat_all8 <- subset(dat_singlet_doubRem_poorRem_ex, Population %in% c("VU", "MGH", "BR1", "DS3", 
                                                                     "DS6", "DS7", "DS8", "DS9"))
all8_plotDat <- FetchData(dat_all8, vars = c("UMAP_1", "UMAP_2", "Population", "Phase"))

plt_UMAP_CCphase <- ggplot(data = all8_plotDat, 
                           aes(x = UMAP_1, y = UMAP_2)) + theme_bw() +
  geom_density_2d(aes(color = Population), alpha = 0.8) +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid",
                                "seagreen", "gold", "green", "blue")) +
  geom_point(shape = 21, aes(fill = Phase), size = 0.8, stroke = 0.01) + 
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12)) 

plt_UMAP_CCphase + ggsave("FIG_S11.svg", 
                          width = 4, height = 3)
plt_UMAP_CCphase_leg <- ggpubr::get_legend(plt_UMAP_CCphase)
as_ggplot(plt_UMAP_CCphase_leg) + ggsave("UMAP_CCphaseColored_all8_leg.svg", 
                                         width = 2.5, height = 4)

## T-SNE - by Population
dat_singlet_doubRem_poorRem_ex$Population <- factor(dat_singlet_doubRem_poorRem_ex$Population,
                                                    levels = c("BR1", "MGH", "VU", "DS3",
                                                               "DS6", "DS7", "DS8", "DS9"))

plt_tsne <- DimPlot(dat_singlet_doubRem_poorRem_ex, reduction = "tsne", group.by = "Population") +
  scale_color_manual(values = c("red", "green", "blue", "brown", "deeppink", 
                                "darkorchid", "seagreen", "gold"),
                     labels = c("PC9-BR1", "PC9-MGH", "PC9-VU", "DS3",
                                "DS6", "DS7", "DS8", "DS9")) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size=12), legend.text = element_text(size=12),
        axis.title = element_text(size=12), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plt_tsne + ggsave("TSNE_PopColored.svg", width = 4, height = 3)
plt_tsne_leg <- ggpubr::get_legend(plt_tsne)
as_ggplot(plt_tsne_leg) + ggsave("FIG_S6F.svg",
                                 width = 2.5, height = 4)

## PCA - by Population
plt_pca <- DimPlot(dat_singlet_doubRem_poorRem_ex, reduction = "pca", group.by = "Population") +
  scale_color_manual(values = c("red", "green", "blue", "brown", "deeppink", 
                                "darkorchid", "seagreen", "gold"),
                     labels = c("PC9-BR1", "PC9-MGH", "PC9-VU", "DS3",
                                "DS6", "DS7", "DS8", "DS9")) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size=12), legend.text = element_text(size=12),
        axis.title = element_text(size=12), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plt_pca + ggsave("FIG_S6E.svg", width = 4, height = 3)


# Identify optimal number of clusters
res_CLV <- NbClust(CLV_plotDat[,c(1:2)], distance = "euclidean", min.nc=2, max.nc=6, 
                   method = "ward.D2", index = "all")
res_sublines <- NbClust(sublines_plotDat[,c(1:2)], distance = "euclidean", min.nc=2, max.nc=6, 
                        method = "ward.D2", index = "all")

# Kmeans overlay
CLV_plotDat$Cluster <- res_CLV$Best.partition
sublines_plotDat$Cluster <- res_sublines$Best.partition

## Plots as UMAP
plt_CLV_cluster <- ggplot(data = CLV_plotDat, aes(x=UMAP_1, y=UMAP_2)) + theme_bw() +
  geom_density_2d(aes(color = Population), 
                  alpha = 0.8) +
  xlim(-10,16) + ylim(-12,14) +
  scale_color_manual(values = c("red", "green", "blue")) +
  geom_point(aes(fill = as.factor(Cluster)), shape = 21, alpha = 1,
             size = 0.8, stroke = 0.01) + 
  scale_fill_manual(values = c("grey80", "grey50", "grey20")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12))

plt_CLV_cluster + ggsave("FIG_S6A.svg", 
                         width = 4, height = 3)
plt_CLV_cluster_leg <- ggpubr::get_legend(plt_CLV_cluster)
as_ggplot(plt_CLV_cluster_leg) + ggsave("FIG_S6A_leg.svg", 
                                        width = 2.5, height = 4)

plt_sublines_cluster <- ggplot(data = sublines_plotDat, aes(x=UMAP_1, y=UMAP_2)) + theme_bw() +
  geom_density_2d(aes(color = Population), 
                  alpha = 0.8) +
  xlim(-10,16) + ylim(-12,14) +
  scale_color_manual(values = c("brown", "deeppink", "darkorchid",
                                "seagreen", "gold")) +
  geom_point(aes(fill = as.factor(Cluster)), shape = 21, alpha = 1,
             size = 0.8, stroke = 0.01) + 
  scale_fill_manual(values = c("grey80", "grey20")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12))

plt_sublines_cluster + ggsave("FIG_S6C.svg", 
                              width = 4, height = 3)
plt_sublines_cluster_leg <- ggpubr::get_legend(plt_sublines_cluster)
as_ggplot(plt_sublines_cluster_leg) + ggsave("FIG_S6C_leg.svg", 
                                             width = 2.5, height = 4)


## Plots as bar
ggplot(data = CLV_plotDat, aes(Cluster)) +
  geom_bar(aes(fill = Population), color = "black", position = "fill", size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("red", "green", "blue")) +
  xlab("Cluster") + ylab("Cluster Fraction") +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12)) +
  ggsave("FIG_S6B.pdf", width = 5, height = 2.5)

ggplot(data = sublines_plotDat, aes(Cluster)) +
  geom_bar(aes(fill = Population), color = "black", position = "fill", size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid",
                               "seagreen", "gold")) +
  xlab("Cluster") + ylab("Cluster Fraction") +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  scale_x_continuous(breaks = seq(0,2,1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12)) +
  ggsave("FIG_S6D.pdf", width = 5, height = 2.5)


# VISION Analysis - functional interpretation
setwd('~/git/GES_2020/scRNAseq/VISION_gmt/hallmark/')
counts <- dat_singlet_doubRem_poorRem_ex@assays$RNA@counts
meta <- dat_singlet_doubRem_poorRem_ex@meta.data
# Scale counts within a sample
n.umi <- colSums(counts)
scaled_counts <- t(t(counts) / n.umi) * median(n.umi)
## HALLMARK GENES
gmt_files <- list.files(path = '~/git/GES_2020/scRNAseq/VISION_gmt/hallmark/', pattern = "\\.gmt$")
vis <- Vision(data = counts, signatures = gmt_files, meta = meta)
vis <- analyze(vis)
visScores <- as.data.frame(getSignatureScores(vis))

## Create dataframe of gene expression and signatures
dat_singlet_cleaned_VISION <- dat_singlet_doubRem_poorRem_ex
dat_singlet_cleaned_VISION <- AddMetaData(object = dat_singlet_cleaned_VISION, 
                                          metadata = visScores)

## Subset dataframe for only UMAP axes, sample name, and hallmark gene signatures
hallmark_metrics <- paste0("HALLMARK_", removeExt(gmt_files))
all_plotDat <- FetchData(dat_singlet_cleaned_VISION, 
                         vars = c("UMAP_1", "UMAP_2", "Population",
                                  hallmark_metrics, "HALLMARK_KRAS_SIGNALING",
                                  "HALLMARK_UV_RESPONSE"))

## Plot distribution of all 48 hallmark sets for each population
all_dat_hallmarks <- all_plotDat[3:ncol(all_plotDat)]
all_dat_hallmarks_melt <- reshape2::melt(all_dat_hallmarks, id.vars = "Population")
all_dat_hallmarks_melt$variable <- str_remove(as.character(all_dat_hallmarks_melt$variable), "HALLMARK_")
all_dat_hallmarks_melt$Population <- factor(all_dat_hallmarks_melt$Population,
                                            levels = c("BR1", "MGH", "VU", "DS3",
                                                       "DS6", "DS7", "DS8", "DS9"))

plt_hallmarks <- ggplot(all_dat_hallmarks_melt, aes(x=value, color = Population, fill = Population)) +
  geom_density(alpha = 0.5) + facet_wrap(~variable, ncol = 6, scales = "free") + theme_bw() +
  scale_color_manual(values = c("red", "green", "blue", "brown", "deeppink", 
                                "darkorchid", "seagreen", "gold"),
                     labels = c("PC9-BR1", "PC9-MGH", "PC9-VU", "DS3",
                                "DS6", "DS7", "DS8", "DS9")) +
  scale_fill_manual(values = c("red", "green", "blue", "brown", "deeppink", 
                               "darkorchid", "seagreen", "gold"),
                    labels = c("PC9-BR1", "PC9-MGH", "PC9-VU", "DS3",
                               "DS6", "DS7", "DS8", "DS9")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", axis.text = element_text(size = 18), axis.title = element_text(size = 18)) +
  xlab("Signature Score") + ylab("Density") 

plt_hallmarks + ggsave("FIG_S8.svg", width = 18, height = 24)
plt_hallmarks_leg <- ggpubr::get_legend(plt_hallmarks)
as_ggplot(plt_hallmarks_leg) + ggsave("FIG_S8_legend.svg",
                                      width = 6, height = 4)

## Plot example UMAP signature overlays
plt_all_overlay <- ggplot() + theme_bw() +
  geom_density_2d(data = all_plotDat, 
                  aes(x = UMAP_1, y = UMAP_2,
                      color = Population), alpha = 0.8) +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid",
                                "seagreen", "gold", "green", "blue")) +
  geom_point(data = all_plotDat, shape = 21,
             aes(x = UMAP_1, y = UMAP_2, fill = HALLMARK_INTERFERON_ALPHA_RESPONSE),
             size = 0.8, stroke = 0.01) + 
  scale_fill_gradient(guide = FALSE, 
                      low = "white", 
                      high = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), legend.position = "none", 
        legend.text = element_text(size = 12), legend.title = element_text(size=12), 
        axis.title=element_text(size=12)) 

plt_all_overlay


# Cluster biomarker identification
## Find differentially expressed genes between different populations
### Cell line versions
dat_CLV <- dat_singlet_doubRem_poorRem_ex
Idents(dat_CLV) <- "Population"

tags_CLV <- c("VU", "MGH", "BR1")
dat_CLV <- SubsetData(dat_CLV,ident.use = tags_CLV)
marker_genes_CLV = list()
for (tag in seq(length(tags_CLV))) {
  marker_dat <- FindMarkers(dat_CLV,
                            ident.1 = tags_CLV[tag],
                            ident.2 = NULL, only.pos = FALSE)
  marker_dat$Population <- tags_CLV[tag]
  marker_genes_CLV[[tag]] <- marker_dat
}

VU_DEGs <- rownames(subset(marker_genes_CLV[[1]], p_val_adj < 0.05))
MGH_DEGs <- rownames(subset(marker_genes_CLV[[2]], p_val_adj < 0.05))
BR1_DEGs <- rownames(subset(marker_genes_CLV[[3]], p_val_adj < 0.05))

#### Create DEG list of all CLV
samples_CLV <- list("VU" = VU_DEGs, 
                    "MGH" = MGH_DEGs, 
                    "BR1" = BR1_DEGs)

### Sublines
dat_sublines <- dat_singlet_doubRem_poorRem_ex
Idents(dat_sublines) <- "Population"

tags_sublines <- c('DS3','DS6','DS7','DS8','DS9')
dat_sublines <- SubsetData(dat_sublines,ident.use = tags_sublines)
marker_genes_sublines = list()
for (tag in seq(length(tags_sublines))) {
  marker_dat <- FindMarkers(dat_sublines,
                            ident.1 = tags_sublines[tag],
                            ident.2 = NULL, only.pos = FALSE)
  marker_dat$Population <- tags_sublines[tag]
  marker_genes_sublines[[tag]] <- marker_dat
}

DS3_DEGs <- rownames(subset(marker_genes_sublines[[1]], p_val_adj < 0.05))
DS6_DEGs <- rownames(subset(marker_genes_sublines[[2]], p_val_adj < 0.05))
DS7_DEGs <- rownames(subset(marker_genes_sublines[[3]], p_val_adj < 0.05))
DS8_DEGs <- rownames(subset(marker_genes_sublines[[4]], p_val_adj < 0.05))
DS9_DEGs <- rownames(subset(marker_genes_sublines[[5]], p_val_adj < 0.05))

#### Create DEG list of all CLV
samples_sublines <- list("DS3" = DS3_DEGs,
                         "DS6" = DS6_DEGs,
                         "DS7" = DS7_DEGs,
                         "DS8" = DS8_DEGs,
                         "DS9" = DS9_DEGs)


# # Save all the differentially expressed gene lists for GO semantic similarity analysis
# save(VU_DEGs, MGH_DEGs, BR1_DEGs,
#      DS3_DEGs, DS6_DEGs, DS7_DEGs, DS8_DEGs, DS9_DEGs,
#      file = "~/git/GES_2020/scRNAseq/DEGs_byCohort_hg38.RData")
# load("~/git/GES_2020/scRNAseq/DEGs_byCohort_hg38.RData")

# # Create lists of impactful mutations from GES_2020/WES/WES.R
# ## Load data from that folder
# load('~/git/GES_2020/WES/variants_byCohort.RData')
# ## Create list of IMPACT mutations for CLV
# samples_muts_CLV <- list("VU" = unique((subset(test_s1_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol)),
#                          "MGH" = unique((subset(test_s2_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol)),
#                          "BR1" = unique((subset(test_s3_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol)))
# 
# ## Create list of IMPACT mutations for sublines
# samples_muts_sublines <- list("DS3" = unique(subset(test_s4_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
#                               "DS6" = unique(subset(test_s5_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
#                               "DS7" = unique(subset(test_s6_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
#                               "DS8" = unique(subset(test_s7_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
#                               "DS9" = unique(subset(test_s8_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol))

# # Save all the differentially expressed gene lists for GO semantic similarity analysis
# save(samples_CLV, samples_muts_CLV,
#      samples_sublines, samples_muts_sublines,
#      file = "~/git/GES_2020/GO/mutations_DEGs-hg38.RData")