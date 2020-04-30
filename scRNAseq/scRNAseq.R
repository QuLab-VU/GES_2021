# devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
library(Seurat) # Use version 2.3.4 for the code below to work
library(Matrix)
library(AnnotationDbi)
library("org.Hs.eg.db")

setwd('~/Documents/QuarantaLab/GES_2020/scRNAseq/')

# Read in hashtag oligonucleotide (HTO) library
## This matrix includes all 8 HTOs. This library comes
## from the CITEseq Count function (see Methods and HTO_count.txt file)
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

# Read in single-cell RNAseq count data
## This count matrix will include all 8 samples. Output from Cell Ranger.
##Annotation with hashing labels below.
matrix_dir_rna = "/Volumes/quaranta/Data/RNAseq/PC9_scRNAseq/outs/filtered_feature_bc_matrix/"
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

## Save d for separate ENSG id conversion
d <- as.data.frame(as.matrix(dat_expr))
save(d, dat_HTO, file = "GEdata_ensgid.RData")

# Add gene symbols as row names for single-cell count matrix
## Replace Ensembl gene IDs with HGNC gene symbols (if possible)
d$SYMBOL <- mapIds(org.Hs.eg.db,keys=rownames(dat_expr),
                   column="SYMBOL",keytype="ENSEMBL",multiVals="first")

rownames(dat_expr) <- make.names(ifelse(is.na(d$SYMBOL), 
                                        rownames(dat_expr), 
                                        d$SYMBOL),
                                 unique = TRUE)

# Keep only the same cell barcodes in the HTO dataframe
# between gene expression and HTO libraries
joint_bcs <- intersect(colnames(dat_expr),colnames(dat_HTO))
dat_expr <- dat_expr[,joint_bcs]
dat_HTO <- as.matrix(dat_HTO[,joint_bcs])

# Setup Seurat object
dat_hashtag <- CreateSeuratObject(dat_expr, project = "PC9_hashing")

# Normalize RNA data with log normalization
dat_hashtag <- NormalizeData(dat_hashtag,display.progress = FALSE)

# Identify and scale variable genes 
dat_hashtag <- FindVariableGenes(dat_hashtag,do.plot = T,
                                    display.progress = FALSE)

# length(dat_hashtag@var.genes) # Identify the number of variable genes
dat_hashtag <- ScaleData(dat_hashtag,
                         genes.use = dat_hashtag@var.genes,
                         display.progress = FALSE)

# Add hashtag count dataframe to assay data
dat_hashtag <- SetAssayData(dat_hashtag,
                            assay.type = "HTO",
                            slot = "raw.data",
                            new.data = dat_HTO)

# Apply special normalization method to HTO libraries
dat_hashtag <- NormalizeData(dat_hashtag, assay.type = "HTO",
                             normalization.method = "genesCLR",
                             display.progress = FALSE)

# Basically, add the HTO information to the gene expression data
dat_hashtag <- HTODemux(dat_hashtag,
                        assay.type = "HTO",
                        positive_quantile =  0.99,
                        print.output = FALSE)

# Some sanity checks for multiplets and hashtags
# print (table(dat_hashtag@meta.data$hto_classification_global))
# print (table(dat_hashtag@meta.data$hto_classification))

# Plot hashtag representation for each hashtag library
## Changed axis label colors and subplot titles in external softward
dat_hashtag <- SetAllIdent(dat_hashtag,id = "hash_maxID")
RidgePlot(dat_hashtag,
          features.plot = rownames(GetAssayData(dat_hashtag,assay.type = "HTO")),
          nCol = 3, cols = c("blue", "green", "red", "brown", 
                             "deeppink", "darkorchid", "seagreen", "gold")) + 
  ggsave("FIG_S13.pdf", width = 15, height = 15)


# Annotate gene expression dataset with most likely HTO classification
dat_hashtag <- SetAllIdent(dat_hashtag,"hto_classification_global")

# Do same normalization techniques on singlet-only dataset
## Extract the singlets
dat_singlet <- SubsetData(dat_hashtag,ident.use = "Singlet")

## Select the most variable genes - This is FIG. S6B (saved as 4x5 PDF)
dat_singlet <- FindVariableGenes(dat_singlet, do.plot = TRUE,
                                 display.progress = FALSE)

# Scaling gene expression data
dat_singlet <- ScaleData(dat_singlet, 
                         genes.use = dat_singlet@var.genes,
                         display.progress = FALSE)

# Run PCA
dat_singlet <- RunPCA(dat_singlet,
                      pc.genes = dat_singlet@var.genes,
                      pcs.print = 0)

# Run tSNE
dat_singlet <- RunTSNE(dat_singlet, reduction.use = "pca",
                       dims.use = 1:10)

# Run UMAP
dat_singlet <- RunUMAP(dat_singlet, reduction.use = "pca", dims.use = 1:10)

# Changing Annotations on Plots from hashtag ID to sample name
## Add sample names as metadata
dat_singlet_clusters <- dat_singlet
old_tags <- c('HTO_A-GTCAACTCTTTAGCG',
              'HTO_B-TGATGGCCTATTGGG',
              'HTO_C-TTCCGCCTCTCTTTG',
              'HTO_D-AGTAAGTTCAGCGTA',
              'HTO_E-AAGTATCGTTTCGCA',
              'HTO_F-GGTTGCCAGATGTCA',
              'HTO_G-TGTCTTTCCTGCCAG',
              'HTO_H-CTCCTCTGCAATTAC')
new_tags <- c('VU', 'MGH', 'BR1', 'DS3',
              'DS6', 'DS7', 'DS8', 'DS9')
found <- match(dat_singlet_clusters@meta.data$hto_classification,
               old_tags)
met <- data.frame(ifelse(is.na(found), dat_singlet_clusters@meta.data$hto_classification,
              new_tags[found]), row.names = rownames(dat_singlet_clusters@meta.data))
names(met) <- "Population"
dat_singlet_clusters <- AddMetaData(dat_singlet_clusters, met,
                                    'Population')

# Plot dimensionality reduction approaches
PCAPlot(dat_singlet_clusters, group.by = "Population") + 
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid", 
                                "seagreen", "gold", "green", "blue")) +
  ggsave("FIG_S7C.svg", width = 6, height = 4)

TSNEPlot(dat_singlet_clusters,group.by = "Population") +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid", 
                                "seagreen", "gold", "green", "blue")) +
  ggsave("FIG_S7B.svg", width = 6, height = 4)

DimPlot(dat_singlet_clusters, reduction.use = "umap", group.by = "Population") +
  scale_color_manual(values = c("red", "brown", "deeppink", "darkorchid", 
                                "seagreen", "gold", "green", "blue")) +
  ggsave("FIG_S7A.svg", width = 6, height = 4)

# Make sample ID the main cell identity
dat_singlet_clusters <- SetAllIdent(dat_singlet_clusters, id = "Population")

# Plot cell line versions (CLV) in space of all samples
## Subset only cell line versions
dat_singlet_CLV <- SubsetData(object = dat_singlet_clusters, 
                              ident.use=c("VU", "MGH", "BR1"))
DimPlot(dat_singlet_CLV, reduction.use = "umap", group.by = "Population") +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14), legend.title = element_text(size=14), 
        axis.title=element_text(size=14), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("FIG_3G_left.svg", width = 6, height = 4)

# Plot sublines (not including DS8) in space of all samples
## Subset only sublines (not DS8)
dat_singlet_Sublines_noDS8 <- SubsetData(object = dat_singlet_clusters, 
                              ident.use=c("DS3", "DS6", "DS7", "DS9"))
DimPlot(dat_singlet_Sublines_noDS8, reduction.use = "umap", group.by = "Population") +
  scale_color_manual(values = c("brown", "deeppink", "darkorchid", "gold")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14), legend.title = element_text(size=14), 
        axis.title=element_text(size=14), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("FIG_4G_left.svg", width = 6, height = 4)

# Plot sublines in space of all samples
## Subset only sublines - highlight DS8
dat_singlet_Sublines <- SubsetData(object = dat_singlet_clusters, 
                                   ident.use=c("DS3", "DS6", "DS7", "DS8", "DS9"))
tint <- data.frame(ifelse(dat_singlet_Sublines@meta.data$Population == "DS8", 1, 0.3), 
                   row.names = rownames(dat_singlet_Sublines@meta.data))
names(tint) <- "Tint"
dat_singlet_Sublines <- AddMetaData(dat_singlet_Sublines,
                                    tint, 'Tint')
DimPlot(dat_singlet_Sublines, reduction.use = "umap", group.by = "Population") +
  scale_color_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
        axis.text=element_text(size=14), legend.title = element_text(size=14), 
        axis.title=element_text(size=14), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("RNA_UMAP_PopColored_Sublines_in8.svg", width = 6, height = 4)

## Make non-DS8 cells less transparent
dat_sublines_test <- data.frame(dat_singlet_Sublines@dr$umap@cell.embeddings)
dat_sublines_test$Population <- dat_singlet_Sublines@meta.data$Population[match(row.names(dat_sublines_test),
                                                                                row.names(dat_singlet_Sublines@meta.data))]
dat_sublines_test$Tint <- ifelse(dat_sublines_test$Population == "DS8", 1, 0.15)
cols <- c("DS3" = "brown", "DS6" = "deeppink", "DS7" = "darkorchid", 
            "DS8" = "seagreen", "DS9" = "gold")
labs <- c("DS3", "DS6", "DS7", "DS8", "DS9")

ggplot(dat_sublines_test, aes(x = UMAP1, y = UMAP2, color = Population, alpha = Tint)) +
  geom_point() + scale_alpha_continuous(range = c(0.15,1)) +
  theme_bw() + 
  scale_color_manual(values = cols, labels = labs) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14)) +
  ggsave("FIG_6C.svg", width = 6, height = 4)


# Run centroid calculations in common UMAP space
## Subset only UMAP values
umap_vals <- as.data.frame(GetCellEmbeddings(object = dat_singlet_clusters, reduction.type = "umap"))
umap_vals_df <- dplyr::bind_rows(umap_vals)
colnames(umap_vals_df) <- c("UMAP1", "UMAP2")
pops_all <- subset(dat_singlet_clusters@meta.data, select = c("Population"))
umap_vals_df$Population <- unlist(pops_all)

## Calculate UMAP centroids for each sample
UMAP1_centroid <- tapply(umap_vals_df$UMAP1, umap_vals_df$Population, mean)
UMAP2_centroid <- tapply(umap_vals_df$UMAP2, umap_vals_df$Population, mean)

## Combine centroids into 2D dataframe
centroid <- data.frame(UMAP1 <- UMAP1_centroid, UMAP2 <- UMAP2_centroid)

## Calculate pairwise distances between centroids
distMat <- harrietr::melt_dist(as.matrix(dist(centroid)))
distMat$Measure<- with(distMat, paste0(iso1, "-", iso2)) # Add pairwise labels

## Subset pairwise distances into cohorts (CLV, Sublines)
distMat_CLV <- subset(distMat, Measure %in% c("VU-BR1", "VU-MGH", "MGH-BR1"))
distMat_Sublines <- subset(distMat, 
                           Measure %in% c("DS6-DS3", "DS7-DS3", "DS8-DS3",
                                          "DS9-DS3", "DS7-DS6", "DS8-DS6",
                                          "DS9-DS6", "DS8-DS7", "DS9-DS7",
                                          "DS9-DS8"))
distMat_Sublines_noDS8 <- subset(distMat, 
                                 Measure %in% c("DS6-DS3", "DS7-DS3",
                                                "DS9-DS3", "DS7-DS6",
                                                "DS9-DS6", "DS9-DS7"))

## Plot cohort comparisons
ggplot(distMat_CLV, aes(x=reorder(Measure, dist), y=dist)) + #, fill = Measure)) +
  geom_bar(stat='identity') + coord_flip() + theme_bw() + ylim(0,16) +
  theme(legend.text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
        legend.title = element_text(size=14,face="bold"), axis.title=element_text(size=14),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Comparison", y = "Distance between Centroids") + 
  ggsave("FIG_3H.pdf", width = 7.5, height = 3)

ggplot(distMat_Sublines_noDS8, aes(x=reorder(Measure, dist), y=dist)) + #, fill = Measure)) +
  geom_bar(stat='identity') + coord_flip() + theme_bw() + ylim(0,16) +
  theme(legend.text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
        legend.title = element_text(size=14,face="bold"), axis.title=element_text(size=14),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Comparison", y = "Distance between Centroids") + 
  ggsave("FIG_4H.pdf", width = 7.5, height = 3)

ggplot(distMat_Sublines, aes(x=reorder(Measure, dist), y=dist)) + #, fill = Measure)) +
  geom_bar(stat='identity') + coord_flip() + theme_bw() + ylim(0,16) +
  theme(legend.text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
        legend.title = element_text(size=14,face="bold"), axis.title=element_text(size=14),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "Comparison", y = "Distance between Centroids") + 
  ggsave("FIG_6D.pdf", width = 7.5, height = 3)


# Run clustering on cohorts in a common UMAP space
## Subset only CLV cells
set.seed(1111)
umap_vals_df_CLV <- subset(umap_vals_df, Population %in% c("VU", "MGH", "BR1"))

## Subset only subline cells (no DS8)
umap_vals_df_sublines_noDS8 <- subset(umap_vals_df, Population %in% c("DS3", "DS6", "DS7", "DS9"))
umap_vals_df_sublines_noDS8 <- subset(umap_vals_df_sublines_noDS8, UMAP2 < 5 & UMAP1 < 0)

## Subset only subline cells (including DS8)
umap_vals_df_sublines <- subset(umap_vals_df, Population %in% c("DS3", "DS6", "DS7", "DS8", "DS9"))
umap_vals_df_sublines <- subset(umap_vals_df_sublines, UMAP1 < 0)

## Clustering using k-means (3 for CLV, 3 for sublines - see main text)
CLV_Kmeans <- kmeans(umap_vals_df_CLV[, 1:2], centers = 3)
sublines_noDS8_Kmeans <- kmeans(umap_vals_df_sublines_noDS8[,1:2], centers = 3)
sublines_Kmeans <- kmeans(umap_vals_df_sublines[,1:2], centers = 4)

## Annotate dataframes with cluster number
umap_vals_df_CLV$Kmeans <- CLV_Kmeans$cluster
umap_vals_df_sublines_noDS8$Kmeans <- sublines_noDS8_Kmeans$cluster
umap_vals_df_sublines$Kmeans <- sublines_Kmeans$cluster

## Plot cluster representation by cohort sample
### Cluster order may be slightly different than in paper
ggplot(umap_vals_df_CLV, aes(Kmeans)) + 
  geom_bar(aes(fill = Population), color = "black", position = "fill", size = 0.1) + 
  theme_classic() + 
  scale_fill_manual(values = c("red", "green", "blue")) +
  xlab("Cluster") + ylab("Cluster Fraction") +
  scale_y_continuous(breaks=seq(0,1,.5)) +
  theme(legend.position = "none", axis.text.y=element_text(size=20),
        axis.text.x = element_blank(), axis.title = element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.line.y = element_line(size = 0.1), axis.ticks.y = element_line(size = 0.1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("FIG_3G_right.pdf", width = 10, height = 3.5)

ggplot(umap_vals_df_sublines_noDS8, aes(Kmeans)) + 
  geom_bar(aes(fill = Population), color = "black", position = "fill", size = 0.1) + 
  theme_classic() + 
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "gold")) +
  xlab("Cluster") + ylab("Cluster Fraction") +
  scale_y_continuous(breaks=seq(0,1,.5)) +
  theme(legend.position = "none", axis.text.y=element_text(size=20),
        axis.text.x = element_blank(), axis.title = element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.line.y = element_line(size = 0.1), axis.ticks.y = element_line(size = 0.1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("FIG_4G_right.pdf", width = 10, height = 3.5)

ggplot(umap_vals_df_sublines, aes(Kmeans)) + 
  geom_bar(aes(fill = Population), color = "black", position = "fill", size = 0.1) + 
  theme_classic() + 
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold")) +
  xlab("Cluster") + ylab("Cluster Fraction") +
  scale_y_continuous(breaks=seq(0,1,.5)) +
  theme(legend.position = "none", axis.text.y=element_text(size=16),
        axis.text.x = element_blank(), axis.title = element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.line.y = element_line(size = 0.1), axis.ticks.y = element_line(size = 0.1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("clusterOverlap_allSublines_in8.pdf", width = 9, height = 3)

## Sanity check UMAP plots by cluster ID (not in paper)
ggplot(umap_vals_df_CLV, aes(x = UMAP1, y = UMAP2, color = Kmeans)) +
  geom_point() + 
  theme_bw() 

ggplot(umap_vals_df_sublines_noDS8, aes(x = UMAP1, y = UMAP2, color = Kmeans)) +
  geom_point() + 
  theme_bw() 

ggplot(umap_vals_df_sublines, aes(x = UMAP1, y = UMAP2, color = Kmeans)) +
  geom_point() + 
  theme_bw() 


# Due diligence
## CC regression 
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
## Separate into G2/M and S phase markers
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

## Perform Seurat cell cycle scoring
dat_singlet_clusters <- CellCycleScoring(object = dat_singlet_clusters, 
                                         s.genes = s.genes, 
                                         g2m.genes = g2m.genes,
                                         set.ident = TRUE)

## Plot data in common UMAP space by cell cycle score
### There is some separation, but not between samples.
### Non-regressed data used
DimPlot(dat_singlet_clusters, reduction.use = "umap", group.by = "Phase") +
  ggsave("FIG_S6C.svg", width = 6, height = 4)

## Make cell identities the sample name again
dat_singlet_clusters <- SetAllIdent(object = dat_singlet_clusters, id = "Population")