setwd("~/git/GES_2020/RNAseq/")

library(ggplot2)
library(DESeq2)
library(biomaRt)

# Load counts from featureCounts output - see RNAseq_processing.txt for how they are created
## See output text files in supplementary information files
VU_counts <- read.csv(file = "2959-CH-1-AAGACCGT-GTCGATTG_S10_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")
MGH_counts <- read.csv(file = "2959-CH-2-TTGCGAGA-TATGGCAC_S11_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")
BR1_counts <- read.csv(file = "2959-CH-3-GCAATTCC-CTCGAACA_S12_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")
DS3_counts <- read.csv(file = "2959-CH-4-GAATCCGT-CAACTCCA_S13_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")
DS6_counts <- read.csv(file = "2959-CH-5-CCGCTTAA-GTCATCGT_S14_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")
DS7_counts <- read.csv(file = "2959-CH-6-TACCTGCA-GGACATCA_S15_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")
DS8_counts <- read.csv(file = "2959-CH-7-GTCGATTG-AAGACCGT_S16_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")
DS9_counts <- read.csv(file = "2959-CH-8-TATGGCAC-TTGCGAGA_S17_R1_001.fastq.gz_featurecounts.txt",
                         sep = "", header = T, skip = 1, row.names = "Geneid")

# Compile into common dataframe
PC9_counts <- cbind(VU_counts[6], MGH_counts[6],
                    BR1_counts[6], DS3_counts[6],
                    DS6_counts[6], DS7_counts[6],
                    DS8_counts[6], DS9_counts[6])

# Rename column headers
colnames(PC9_counts) <- c("VU", "MGH", "BR1", "DS3",
                          "DS6", "DS7", "DS8", "DS9")

# Change Ensembl gene ID with gene symbol
PC9_counts$ensembl_gene_id <- rownames(PC9_counts)
ensembl <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
genes <- rownames(PC9_counts)
G_list <- getBM(attributes= c("ensembl_gene_id","hgnc_symbol"),
                filters= "ensembl_gene_id",
                values=genes,
                mart=mart)

# Reconfigure dataframe so that gene symbols are row names and counts
GE_data <- merge(PC9_counts, G_list, by = "ensembl_gene_id")
d <- GE_data[, -1]
d <- d[c(9, seq(1:8))]
rownames(d) <- make.names(d$hgnc_symbol, unique = T)
d <- d[, 2:9]
countdata <- d

# Save data as loadable RData file
# save(countdata, file = "PC9_countData.RData")

# DESeq2 RLD normalized (only one that works for single-replicate in DESeq2)
coldata <- data.frame(sample = colnames(countdata), row.names = colnames(countdata))
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ sample)

## Minimum count filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rld <- rlog(dds, blind=TRUE)

# Run principal component analysis (PCA)
pca_data <- prcomp(t(assay(rld)))
pca_data_perc <- round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
df_pca_data <- data.frame(PC1 = pca_data$x[,1], 
                          PC2 = pca_data$x[,2], 
                          Line = colnames(assay(rld)))

# Plot normalized data
ggplot(df_pca_data, aes(PC1, PC2, fill = Line))+
  geom_point(shape = 21, size=4, color = "black") + theme_bw() +
  labs(x="PC1", y="PC2") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_fill_manual(values = c("red", "brown", "deeppink",
                               "darkorchid", "seagreen",
                               "gold", "green", "blue")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    legend.position = "right", axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14)) +
  ggsave("FIG_S7A.pdf", width = 5, height = 4)

# Hierarchical clustering of PCA transcriptomic phenotypes
dists <- dist(t(assay(rld)))
hc <- hclust(dists, method = "ward.D2")
pdf("FIG_S7B.pdf")
plot(hc, xlab = "Population", main = "Cluster Tree by Transcriptome - PC9")
dev.off()