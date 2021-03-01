library(DESeq2)
library(biomaRt)

setwd('~/Documents/QuarantaLab/GES_2020/scRNAseq/')
load('GEdata_ensgid.RData')
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# Include only cell barcodes in both datasets
dat_expr <- d
joint_bcs <- intersect(colnames(dat_expr),colnames(dat_HTO))
dat_expr <- dat_expr[,joint_bcs]
dat_HTO <- as.matrix(dat_HTO[,joint_bcs])

# Setup Seurat object
dat_hashtag <- CreateSeuratObject(dat_expr, project = "PC9_hashing")

# Normalize RNA data with log normalization
dat_hashtag <- NormalizeData(dat_hashtag,display.progress = FALSE)

# Find and scale variable genes
dat_hashtag <- FindVariableGenes(dat_hashtag,do.plot = T,
                                 display.progress = FALSE)
dat_hashtag <- ScaleData(dat_hashtag,
                         genes.use = dat_hashtag@var.genes,
                         display.progress = FALSE)

# Add hashtag count dataframe to assay data
dat_hashtag <- SetAssayData(dat_hashtag,
                            assay.type = "HTO",
                            slot = "raw.data",
                            new.data = dat_HTO)

# Use special normalization method to hashtag libraries
dat_hashtag <- NormalizeData(dat_hashtag,assay.type = "HTO",
                             normalization.method = "genesCLR",
                             display.progress = FALSE)

# Basically, add the HTO information to the GE matrix
dat_hashtag <- HTODemux(dat_hashtag,
                        assay.type = "HTO",
                        positive_quantile =  0.99,
                        print.output = FALSE)

dat_hashtag <- SetAllIdent(dat_hashtag,id = "hash_maxID")
dat_hashtag <- SetAllIdent(dat_hashtag,"hto_classification_global")

# Identification of genes that distinguish samples between cohorts
## Will be used later for genomic-to-transcriptomic (i.e., genetic-to epigenetic) connection
### Cell Line Versions (CLV) - plus same normalization as above
#### Have to do normalization/scaling again because it's a new transcriptomic space
dat_hashtag_CLV <- dat_hashtag
dat_hashtag_CLV <- SubsetData(dat_hashtag_CLV,ident.use = "Singlet")
dat_hashtag_CLV <- SetAllIdent(dat_hashtag_CLV,"hash_ID")
dat_hashtag_CLV <- SubsetData(dat_hashtag_CLV,
                              ident.use = c('HTO_A-GTCAACTCTTTAGCG',
                                            'HTO_B-TGATGGCCTATTGGG',
                                            'HTO_C-TTCCGCCTCTCTTTG'))

dat_hashtag_CLV <- FindVariableGenes(dat_hashtag_CLV,do.plot = F)
dat_hashtag_CLV <- ScaleData(dat_hashtag_CLV,
                             genes.use = dat_hashtag_CLV@var.genes,
                             display.progress = FALSE)

### Changing from HTO identifiers to sample names
old_tags_CLV <- c('HTO_A-GTCAACTCTTTAGCG',
                  'HTO_B-TGATGGCCTATTGGG',
                  'HTO_C-TTCCGCCTCTCTTTG')
new_tags_CLV <- c('VU', 'MGH', 'BR1')
found_CLV <- match(dat_hashtag_CLV@meta.data$hash_ID,
                   old_tags_CLV)

# Add names as metadata tag
met_CLV <- data.frame(ifelse(is.na(found_CLV), 
                             dat_hashtag_CLV@meta.data$hash_ID,
                             new_tags_CLV[found_CLV]), 
                      row.names = rownames(dat_hashtag_CLV@meta.data))
names(met_CLV) <- "Population"
dat_hashtag_CLV <- AddMetaData(dat_hashtag_CLV,
                               met_CLV, 'Population')

### Get DEGs for each cell line version compared to other versions
dat_hashtag_CLV <- SetAllIdent(object = dat_hashtag_CLV, id = "Population")

marker_genes_CLV = list()
for (tag in seq(length(new_tags_CLV))) {
  marker_dat <- FindMarkers(dat_hashtag_CLV,
                            ident.1 = new_tags_CLV[tag],
                            ident.2 = NULL, only.pos = FALSE)
  marker_dat$Population <- new_tags_CLV[tag]
  marker_genes_CLV[[tag]] <- marker_dat
}

### Convert sample differentially expressed genes (DEGs) from Ensembl
### gene IDs to HGNC symbols
VU_DEGs <- getBM(attributes='hgnc_symbol', 
                 filters = 'ensembl_gene_id', 
                 values = rownames(marker_genes_CLV[[1]]), 
                 mart = ensembl)

MGH_DEGs <- getBM(attributes='hgnc_symbol', 
                  filters = 'ensembl_gene_id', 
                  values = rownames(marker_genes_CLV[[2]]), 
                  mart = ensembl)

BR1_DEGs <- getBM(attributes='hgnc_symbol', 
                  filters = 'ensembl_gene_id', 
                  values = rownames(marker_genes_CLV[[3]]), 
                  mart = ensembl)

### Sublines - same as above
dat_hashtag_sublines <- dat_hashtag
dat_hashtag_sublines <- SubsetData(dat_hashtag_sublines,ident.use = "Singlet")
dat_hashtag_sublines <- SetAllIdent(dat_hashtag_sublines,"hash_ID")
dat_hashtag_sublines <- SubsetData(dat_hashtag_sublines,
                                   ident.use = c('HTO_D-AGTAAGTTCAGCGTA',
                                                 'HTO_E-AAGTATCGTTTCGCA',
                                                 'HTO_F-GGTTGCCAGATGTCA',
                                                 'HTO_G-TGTCTTTCCTGCCAG',
                                                 'HTO_H-CTCCTCTGCAATTAC'))
dat_hashtag_sublines <- FindVariableGenes(dat_hashtag_sublines,
                                          do.plot = F)
dat_hashtag_sublines <- ScaleData(dat_hashtag_sublines,
                                  genes.use = dat_hashtag_sublines@var.genes,
                                  display.progress = FALSE)

old_tags_sublines <- c('HTO_D-AGTAAGTTCAGCGTA',
                       'HTO_E-AAGTATCGTTTCGCA',
                       'HTO_F-GGTTGCCAGATGTCA',
                       'HTO_G-TGTCTTTCCTGCCAG',
                       'HTO_H-CTCCTCTGCAATTAC')
new_tags_sublines <- c('DS3', 'DS6', 'DS7', 'DS8', 'DS9')
found_sublines <- match(dat_hashtag_sublines@meta.data$hash_ID,
                        old_tags_sublines)
met_sublines <- data.frame(ifelse(is.na(found_sublines), dat_hashtag_sublines@meta.data$hash_ID,
                                  new_tags_sublines[found_sublines]), row.names = rownames(dat_hashtag_sublines@meta.data))
names(met_sublines) <- "Population"
dat_hashtag_sublines <- AddMetaData(dat_hashtag_sublines,
                                    met_sublines,
                                    'Population')
dat_hashtag_sublines <- SetAllIdent(object = dat_hashtag_sublines, 
                                    id = "Population")

marker_genes_sublines = list()
for (tag in seq(length(new_tags_sublines))) {
  marker_dat <- FindMarkers(dat_hashtag_sublines,
                            ident.1 = new_tags_sublines[tag],
                            ident.2 = NULL, only.pos = FALSE)
  marker_dat$Population <- new_tags_sublines[tag]
  marker_genes_sublines[[tag]] <- marker_dat
}

DS3_DEGs <- getBM(attributes='hgnc_symbol', 
                  filters = 'ensembl_gene_id', 
                  values = rownames(marker_genes_sublines[[1]]), 
                  mart = ensembl)

DS6_DEGs <- getBM(attributes='hgnc_symbol', 
                  filters = 'ensembl_gene_id', 
                  values = rownames(marker_genes_sublines[[2]]), 
                  mart = ensembl)

DS7_DEGs <- getBM(attributes='hgnc_symbol', 
                  filters = 'ensembl_gene_id', 
                  values = rownames(marker_genes_sublines[[3]]), 
                  mart = ensembl)

DS8_DEGs <- getBM(attributes='hgnc_symbol', 
                  filters = 'ensembl_gene_id', 
                  values = rownames(marker_genes_sublines[[4]]), 
                  mart = ensembl)

DS9_DEGs <- getBM(attributes='hgnc_symbol', 
                  filters = 'ensembl_gene_id', 
                  values = rownames(marker_genes_sublines[[5]]), 
                  mart = ensembl)

### Sublines without DS8 - same as above
dat_hashtag_sublines_noDS8 <- dat_hashtag
dat_hashtag_sublines_noDS8 <- SubsetData(dat_hashtag_sublines_noDS8,ident.use = "Singlet")
dat_hashtag_sublines_noDS8 <- SetAllIdent(dat_hashtag_sublines_noDS8,"hash_ID")
dat_hashtag_sublines_noDS8 <- SubsetData(dat_hashtag_sublines_noDS8,
                                         ident.use = c('HTO_D-AGTAAGTTCAGCGTA',
                                                       'HTO_E-AAGTATCGTTTCGCA',
                                                       'HTO_F-GGTTGCCAGATGTCA',
                                                       'HTO_H-CTCCTCTGCAATTAC'))

dat_hashtag_sublines_noDS8 <- FindVariableGenes(dat_hashtag_sublines_noDS8,do.plot = F)
dat_hashtag_sublines_noDS8 <- ScaleData(dat_hashtag_sublines_noDS8,
                                        genes.use = dat_hashtag_sublines_noDS8@var.genes,
                                        display.progress = FALSE)

old_tags_sublines_noDS8 <- c('HTO_D-AGTAAGTTCAGCGTA',
                             'HTO_E-AAGTATCGTTTCGCA',
                             'HTO_F-GGTTGCCAGATGTCA',
                             'HTO_H-CTCCTCTGCAATTAC')
new_tags_sublines_noDS8 <- c('DS3', 'DS6', 'DS7', 'DS9')
found_sublines_noDS8 <- match(dat_hashtag_sublines_noDS8@meta.data$hash_ID,
                              old_tags_sublines_noDS8)
met_sublines_noDS8 <- data.frame(ifelse(is.na(found_sublines_noDS8), dat_hashtag_sublines_noDS8@meta.data$hash_ID,
                                        new_tags_sublines_noDS8[found_sublines_noDS8]), row.names = rownames(dat_hashtag_sublines_noDS8@meta.data))
names(met_sublines_noDS8) <- "Population"
dat_hashtag_sublines_noDS8 <- AddMetaData(dat_hashtag_sublines_noDS8,
                                          met_sublines_noDS8,
                                          'Population')

dat_hashtag_sublines_noDS8 <- SetAllIdent(object = dat_hashtag_sublines_noDS8, id = "Population")
marker_genes_sublines_noDS8 = list()
for (tag in seq(length(new_tags_sublines_noDS8))) {
  marker_dat <- FindMarkers(dat_hashtag_sublines_noDS8,
                            ident.1 = new_tags_sublines_noDS8[tag],
                            ident.2 = NULL, only.pos = FALSE)
  marker_dat$Population <- new_tags_sublines_noDS8[tag]
  marker_genes_sublines_noDS8[[tag]] <- marker_dat
}

DS3_DEGs_noDS8 <- getBM(attributes='hgnc_symbol', 
                        filters = 'ensembl_gene_id', 
                        values = rownames(marker_genes_sublines_noDS8[[1]]), 
                        mart = ensembl)

DS6_DEGs_noDS8 <- getBM(attributes='hgnc_symbol', 
                        filters = 'ensembl_gene_id', 
                        values = rownames(marker_genes_sublines_noDS8[[2]]), 
                        mart = ensembl)

DS7_DEGs_noDS8 <- getBM(attributes='hgnc_symbol', 
                        filters = 'ensembl_gene_id', 
                        values = rownames(marker_genes_sublines_noDS8[[3]]), 
                        mart = ensembl)

DS9_DEGs_noDS8 <- getBM(attributes='hgnc_symbol', 
                        filters = 'ensembl_gene_id', 
                        values = rownames(marker_genes_sublines_noDS8[[4]]), 
                        mart = ensembl)


# Save all the differentially expressed gene lists for GO semantic similarity analysis
save(VU_DEGs, MGH_DEGs, BR1_DEGs,
     DS3_DEGs_noDS8, DS6_DEGs_noDS8, DS7_DEGs_noDS8, DS9_DEGs_noDS8,
     DS3_DEGs, DS6_DEGs, DS7_DEGs, DS8_DEGs, DS9_DEGs,
     file = "../GO/DEGs_byCohort.RData")