library(BiocManager)
library(infercnv)
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd('~/git/GES_2020/scRNAseq/inferCNV')
# load("Seurat_v3_scRNAseq_forInferCNV.RData")

# # Code used to create inferCNV input files
# counts_matrix = GetAssayData(dat_singlet_doubRem_poorRem, 
#                              slot="counts")
# rownames(counts_matrix) = make.names(rownames(counts_matrix), unique=TRUE)
# 
# annot_file <- data.frame(rownames(dat_singlet_doubRem_poorRem@meta.data),
#                          dat_singlet_doubRem_poorRem@meta.data$Population)
# 
# write.table(x = annot_file, file = '~/Documents/QuarantaLab/PC9.10x.annotation.txt',
#           col.names = F, row.names = F, quote = F, sep = "\t")
# 
# a <- read.table("Homo_sapiens.GRCh38.98_geneOrder.txt")
# a$V5 <- mapIds(org.Hs.eg.db, keys = a$V1,
#                column = "SYMBOL", keytype = "ENSEMBL",
#                multiVals = "first")
# a$V6 <- make.names(ifelse(is.na(a$V5),
#                           a$V1,
#                           a$V5),
#                    unique = T)
# b <- a[,c("V6","V2","V3","V4")]
# 
# write.table(x = b, file = "Homo_sapiens.GRCh38.98_geneOrder_symbolAnnot.txt",
#             quote = F, sep = "\t", row.names = F, col.names = F)
# 


## Cell line versions with no background
# # Making input files
# annot_file_CLV <- subset(annot_file, dat_singlet_doubRem_poorRem@meta.data$Population %in% c("VU", "MGH", "BR1"))
# write.table(x = annot_file_CLV, file = '~/Documents/QuarantaLab/PC9.CLV.10x.annotation.txt',
#             col.names = F, row.names = F, quote = F, sep = "\t")
# counts_matrix_CLV <- counts_matrix[,colnames(counts_matrix) %in% 
#                                      unique(annot_file_CLV$rownames.dat_singlet_doubRem_poorRem.meta.data.)]
# saveRDS(round(counts_matrix_CLV, digits=3), "PC9.CLV.10x.counts.matrix.rds")

# Create inferCNV object (with no reference group)
infercnv_obj_CLV = infercnv::CreateInfercnvObject(raw_counts_matrix='PC9.CLV.10x.counts.matrix.rds',
                                                  annotations_file='PC9.CLV.10x.annotation.txt',
                                                  delim="\t",
                                                  gene_order_file='Homo_sapiens.GRCh38.98_geneOrder_symbolAnnot.txt',
                                                  ref_group_names = NULL)

# Run inferCNV 
infercnv_obj_CLV = infercnv::run(infercnv_obj_CLV,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir="PC9_cnvs_CLV",  # dir is auto-created for storing outputs
                                 cluster_by_groups=T,   # cluster
                                 denoise=T,
                                 HMM=T
)

## DS lines with VU as background
# # Making input files
# annot_file_VUDS <- subset(annot_file, dat_singlet_doubRem_poorRem@meta.data$Population %in% 
#                             c("VU", "DS3", "DS6", "DS7", "DS8", "DS9"))
# write.table(x = annot_file_VUDS, file = '~/Documents/QuarantaLab/PC9.VUDS.10x.annotation.txt',
#             col.names = F, row.names = F, quote = F, sep = "\t")
# counts_matrix_VUDS <- counts_matrix[,colnames(counts_matrix) %in% 
#                                       unique(annot_file_VUDS$rownames.dat_singlet_doubRem_poorRem.meta.data.)]
# saveRDS(round(counts_matrix_VUDS, digits=3), "PC9.VUDS.10x.counts.matrix.rds")

# Create inferCNV object (with VU as reference group)
infercnv_obj_VUDS = infercnv::CreateInfercnvObject(raw_counts_matrix='PC9.VUDS.10x.counts.matrix.rds',
                                              annotations_file='PC9.VUDS.10x.annotation.txt',
                                              delim="\t",
                                              gene_order_file='Homo_sapiens.GRCh38.98_geneOrder_symbolAnnot.txt',
                                              ref_group_names = "VU")

# Run inferCNV 
infercnv_obj_VUDS = infercnv::run(infercnv_obj_VUDS,
                             cutoff=0.1,  # 0.1 for 10x-genomics
                             out_dir="PC9_cnvs_VUrefGroup_DS",
                             cluster_by_groups=T,   # cluster by population
                             denoise=T,
                             HMM=T
)


## Whole dataset together
# # Create input counts matrix for whole dataset
# # saveRDS(round(counts_matrix, digits=3), "PC9.10x.counts.matrix.rds")
# # 
# # write.table(round(counts_matrix, digits=3), file='PC9.10x.counts.txt.matrix',
# #             quote=F, sep="\t")
# 
# # Create inferCNV object
# infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix='PC9.10x.counts.matrix.rds',
#                                               annotations_file='PC9.10x.annotation.txt',
#                                               delim="\t",
#                                               gene_order_file='Homo_sapiens.GRCh38.98_geneOrder_symbolAnnot.txt',
#                                               ref_group_names = NULL)
# 
# # Run inferCNV pipeline
# infercnv_obj = infercnv::run(infercnv_obj,
#                              cutoff=0.1,  # used 0.1 for 10x-genomics (in tool manual)
#                              out_dir="PC9_cnvs_noRefGroups",  # dir is auto-created for storing outputs
#                              cluster_by_groups=T,   # cluster by population
#                              denoise=T,
#                              HMM=T
# )


# # For plotting in inferCNV (if interested)
# library(tibble)
# library(tsvio)
# library(NGCHM)
# library(infercnvNGCHM)
# 
# ngchm(infercnv_obj = infercnv_obj,
#       out_dir = "PC9_cnvs_noRefGroups",
#       path_to_shaidyMapGen = '~/Downloads/ShaidyMapGen.jar')

