# BiocManager::install("GOSemSim")
library(stringr)
library(enrichR)
library(GOSemSim)
library(reshape2)
library(ggplot2)

# Load unique variants and differentially expressed genes
# for each cohort
setwd('~/Documents/QuarantaLab/GES_2020/GO/')
load('variants_byCohort.RData')
load('DEGs_byCohort.RData')

# Use these enrichr GO databases
dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018")

# Cell Line Versions (CLV)
# Create DEG list of all CLV
samples_CLV <- list("VU" = VU_DEGs$hgnc_symbol, 
                    "MGH" = MGH_DEGs$hgnc_symbol, 
                    "BR1" = BR1_DEGs$hgnc_symbol)

# Create annotated DEG GO dataframe
GO_df_DEG_CLV <- list()
i = 1
for (s in samples_CLV) {
  # print(names(samples)[i])
  enriched <- enrichr(s, dbs)
  bp <- enriched[["GO_Biological_Process_2018"]]
  mf <- enriched[["GO_Molecular_Function_2018"]]
  bp$Sample = names(samples_CLV)[i]
  mf$Sample = names(samples_CLV)[i]
  bp$GOtype <- "bp"
  mf$GOtype <- "mf"
  GO_samp <- rbind(bp, mf)
  GO_df_DEG_CLV[[names(samples_CLV)[i]]] <- GO_samp
  i = i + 1
}

## Add columns to be subset on later
allGO_data_DEG_CLV = do.call(rbind, GO_df_DEG_CLV)
allGO_data_DEG_CLV$logp <- -log10(allGO_data_DEG_CLV$P.value)
allGO_data_DEG_CLV$logq <- -log10(allGO_data_DEG_CLV$Adjusted.P.value)
allGO_data_DEG_CLV$numGO <- str_extract(allGO_data_DEG_CLV$Term, "GO:[0-9]{1,}")

## Keep only significant terms
sig_GO_CLV <- subset(allGO_data_DEG_CLV, Adjusted.P.value < 0.05)

## Significant GO terms by CLV and GO type for scRNAseq DEGs
sig_GO_DF_VU <- subset(sig_GO_CLV, Sample == "VU")
sig_GO_DF_VU_BP <- subset(sig_GO_DF_VU, GOtype == "bp")
sig_GO_DF_VU_MF <- subset(sig_GO_DF_VU, GOtype == "mf")

sig_GO_DF_MGH <- subset(sig_GO_CLV, Sample == "MGH")
sig_GO_DF_MGH_BP <- subset(sig_GO_DF_MGH, GOtype == "bp")
sig_GO_DF_MGH_MF <- subset(sig_GO_DF_MGH, GOtype == "mf")

sig_GO_DF_BR1 <- subset(sig_GO_CLV, Sample == "BR1")
sig_GO_DF_BR1_BP <- subset(sig_GO_DF_BR1, GOtype == "bp")
sig_GO_DF_BR1_MF <- subset(sig_GO_DF_BR1, GOtype == "mf")

# Create list of IMPACT mutations for CLV
samples_muts_CLV <- list("VU" = unique((subset(test_s1_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol)),
                         "MGH" = unique((subset(test_s2_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol)),
                         "BR1" = unique((subset(test_s3_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol)))

## Create dataframe of GO terms associated with cell line version mutations
muts_GO_CLV <- list()
i = 1
for (s in samples_muts_CLV) {
  enriched <- enrichr(s, dbs)
  bp <- enriched[["GO_Biological_Process_2018"]]
  mf <- enriched[["GO_Molecular_Function_2018"]]
  bp$Sample = names(samples_muts_CLV)[i]
  mf$Sample = names(samples_muts_CLV)[i]
  bp$GOtype <- "bp"
  mf$GOtype <- "mf"
  GO_samp <- rbind(bp, mf)
  muts_GO_CLV[[names(samples_muts_CLV)[i]]] <- GO_samp
  i = i + 1
}

## Annotate dataframe with key terms (subset on below)
muts_GO_all_CLV = do.call(rbind, muts_GO_CLV)
muts_GO_all_CLV$logp <- -log10(muts_GO_all_CLV$P.value)
muts_GO_all_CLV$logq <- -log10(muts_GO_all_CLV$Adjusted.P.value)
muts_GO_all_CLV$numGO <- str_extract(muts_GO_all_CLV$Term, "GO:[0-9]{1,}")

## Keep only significant GO terms
muts_GO_sig_CLV <- subset(muts_GO_all_CLV, P.value < 0.05)

## Identify significant GO terms by Sample and GO type for mutations
muts_GO_sig_VU <- subset(muts_GO_sig_CLV, Sample == "VU")
muts_GO_sig_VU_BP <- subset(muts_GO_sig_VU, GOtype == "bp")
muts_GO_sig_VU_MF <- subset(muts_GO_sig_VU, GOtype == "mf")

muts_GO_sig_MGH <- subset(muts_GO_sig_CLV, Sample == "MGH")
muts_GO_sig_MGH_BP <- subset(muts_GO_sig_MGH, GOtype == "bp")
muts_GO_sig_MGH_MF <- subset(muts_GO_sig_MGH, GOtype == "mf")

muts_GO_sig_BR1 <- subset(muts_GO_sig_CLV, Sample == "BR1")
muts_GO_sig_BR1_BP <- subset(muts_GO_sig_BR1, GOtype == "bp")
muts_GO_sig_BR1_MF <- subset(muts_GO_sig_BR1, GOtype == "mf")


## GOSemSim comparison 
### Load databases
hsGO_bp <- godata('org.Hs.eg.db', ont="BP")
hsGO_mf <- godata('org.Hs.eg.db', ont="MF")

### Genomic --> transcriptomic connection 
#### Will be 6 comparisons - 3 groups, 2 ontologies each
#### VU
VU_BP <- mgoSim(muts_GO_sig_VU_BP$numGO, sig_GO_DF_VU_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
VU_MF <- mgoSim(muts_GO_sig_VU_MF$numGO, sig_GO_DF_VU_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

#### MGH
MGH_BP <- mgoSim(muts_GO_sig_MGH_BP$numGO, sig_GO_DF_MGH_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
MGH_MF <- mgoSim(muts_GO_sig_MGH_MF$numGO, sig_GO_DF_MGH_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

#### BR1
BR1_BP <- mgoSim(muts_GO_sig_BR1_BP$numGO, sig_GO_DF_BR1_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
BR1_MF <- mgoSim(muts_GO_sig_BR1_MF$numGO, sig_GO_DF_BR1_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

### Compile dataframe of similarities
semSimDF_CLV <- data.frame(BP = c(VU_BP, MGH_BP, BR1_BP),
                            MF = c(VU_MF, MGH_MF, BR1_MF))
semSimDF_CLV$id <- c("VU", "MGH", "BR1")

### Unique GO ontology type similarity plots 
#### BP
ssDF_CLV_BP <- melt(semSimDF_CLV, measure.vars = c("BP"), id.vars = "id")
ggplot(ssDF_CLV_BP, aes(x = id, y = value, fill = id)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Cell Line Version", y = "GO Semantic Similarity") +
  ylim(0,0.5) +
  scale_fill_manual(values = c("red", "green", "blue"),
                    name = "none") +
  ggtitle("Biological Process") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("FIG_3I_left.pdf", width = 4, height = 5)

#### MF
ssDF_CLV_MF <- melt(semSimDF_CLV, measure.vars = c("MF"), id.vars = "id")
ggplot(ssDF_CLV_MF, aes(x = id, y = value, fill = id)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Cell Line Version", y = "GO Semantic Similarity") +
  ylim(0,0.35) +
  scale_fill_manual(values = c("red", "green", "blue"),
                    name = "none") +
  ggtitle("Molecular Function") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("FIG_3I_right.pdf", width = 4, height = 5)


# Sublines 
## Create a dataframe of all subline differentially expressed genes
samples_sublines <- list("DS3" = DS3_DEGs$hgnc_symbol,
                         "DS6" = DS6_DEGs$hgnc_symbol,
                         "DS7" = DS7_DEGs$hgnc_symbol,
                         "DS8" = DS8_DEGs$hgnc_symbol,
                         "DS9" = DS9_DEGs$hgnc_symbol)

## Create GO term dataframe for each GO type and subline
GO_df_DEG_sublines <- list()
i = 1
for (s in samples_sublines) {
  enriched <- enrichr(s, dbs)
  bp <- enriched[["GO_Biological_Process_2018"]]
  mf <- enriched[["GO_Molecular_Function_2018"]]
  bp$Sample = names(samples_sublines)[i]
  mf$Sample = names(samples_sublines)[i]
  bp$GOtype <- "bp"
  mf$GOtype <- "mf"
  GO_samp <- rbind(bp, mf)
  GO_df_DEG_sublines[[names(samples_sublines)[i]]] <- GO_samp
  i = i + 1
}

## Compile into common subline dataframe
allGO_data_DEG_sublines = do.call(rbind, GO_df_DEG_sublines)
allGO_data_DEG_sublines$logp <- -log10(allGO_data_DEG_sublines$P.value)
allGO_data_DEG_sublines$logq <- -log10(allGO_data_DEG_sublines$Adjusted.P.value)
allGO_data_DEG_sublines$numGO <- str_extract(allGO_data_DEG_sublines$Term, "GO:[0-9]{1,}")

## Keep only significant GO terms
sig_GO_sublines <- subset(allGO_data_DEG_sublines, Adjusted.P.value < 0.05)

## Identify significant GO terms by subline GO type for scRNAseq DEGs
sig_GO_DF_DS3 <- subset(sig_GO_sublines, Sample == "DS3")
sig_GO_DF_DS3_BP <- subset(sig_GO_DF_DS3, GOtype == "bp")
sig_GO_DF_DS3_MF <- subset(sig_GO_DF_DS3, GOtype == "mf")

sig_GO_DF_DS6 <- subset(sig_GO_sublines, Sample == "DS6")
sig_GO_DF_DS6_BP <- subset(sig_GO_DF_DS6, GOtype == "bp")
sig_GO_DF_DS6_MF <- subset(sig_GO_DF_DS6, GOtype == "mf")

sig_GO_DF_DS7 <- subset(sig_GO_sublines, Sample == "DS7")
sig_GO_DF_DS7_BP <- subset(sig_GO_DF_DS7, GOtype == "bp")
sig_GO_DF_DS7_MF <- subset(sig_GO_DF_DS7, GOtype == "mf")

sig_GO_DF_DS8 <- subset(sig_GO_sublines, Sample == "DS8")
sig_GO_DF_DS8_BP <- subset(sig_GO_DF_DS8, GOtype == "bp")
sig_GO_DF_DS8_MF <- subset(sig_GO_DF_DS8, GOtype == "mf")

sig_GO_DF_DS9 <- subset(sig_GO_sublines, Sample == "DS9")
sig_GO_DF_DS9_BP <- subset(sig_GO_DF_DS9, GOtype == "bp")
sig_GO_DF_DS9_MF <- subset(sig_GO_DF_DS9, GOtype == "mf")

# Create list of IMPACT mutations for sublines
samples_muts_sublines <- list("DS3" = unique(subset(test_s4_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
                              "DS6" = unique(subset(test_s5_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
                              "DS7" = unique(subset(test_s6_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
                              "DS8" = unique(subset(test_s7_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol),
                              "DS9" = unique(subset(test_s8_dedup, Impact %in% c("LOW", "MODERATE", "HIGH"))$Symbol))

## Create dataframe of GO terms associated with subline mutations
muts_GO_sublines <- list()
i = 1
for (s in samples_muts_sublines) {
  enriched <- enrichr(s, dbs)
  bp <- enriched[["GO_Biological_Process_2018"]]
  mf <- enriched[["GO_Molecular_Function_2018"]]
  bp$Sample = names(samples_muts_sublines)[i]
  mf$Sample = names(samples_muts_sublines)[i]
  bp$GOtype <- "bp"
  mf$GOtype <- "mf"
  GO_samp <- rbind(bp, mf)
  muts_GO_sublines[[names(samples_muts_sublines)[i]]] <- GO_samp
  i = i + 1
}

## Annotate dataframe with key terms (subset on below)
muts_GO_all_sublines = do.call(rbind, muts_GO_sublines)
muts_GO_all_sublines$logp <- -log10(muts_GO_all_sublines$P.value)
muts_GO_all_sublines$logq <- -log10(muts_GO_all_sublines$Adjusted.P.value)
muts_GO_all_sublines$numGO <- str_extract(muts_GO_all_sublines$Term, "GO:[0-9]{1,}")

## Keep only significant GO terms
muts_GO_sig_sublines <- subset(muts_GO_all_sublines, P.value < 0.05)

## Identify significant GO terms by subline and GO type for IMPACT mutations
muts_GO_sig_DS3 <- subset(muts_GO_sig_sublines, Sample == "DS3")
muts_GO_sig_DS3_BP <- subset(muts_GO_sig_DS3, GOtype == "bp")
muts_GO_sig_DS3_MF <- subset(muts_GO_sig_DS3, GOtype == "mf")

muts_GO_sig_DS6 <- subset(muts_GO_sig_sublines, Sample == "DS6")
muts_GO_sig_DS6_BP <- subset(muts_GO_sig_DS6, GOtype == "bp")
muts_GO_sig_DS6_MF <- subset(muts_GO_sig_DS6, GOtype == "mf")

muts_GO_sig_DS7 <- subset(muts_GO_sig_sublines, Sample == "DS7")
muts_GO_sig_DS7_BP <- subset(muts_GO_sig_DS7, GOtype == "bp")
muts_GO_sig_DS7_MF <- subset(muts_GO_sig_DS7, GOtype == "mf")

muts_GO_sig_DS8 <- subset(muts_GO_sig_sublines, Sample == "DS8")
muts_GO_sig_DS8_BP <- subset(muts_GO_sig_DS8, GOtype == "bp")
muts_GO_sig_DS8_MF <- subset(muts_GO_sig_DS8, GOtype == "mf")

muts_GO_sig_DS9 <- subset(muts_GO_sig_sublines, Sample == "DS9")
muts_GO_sig_DS9_BP <- subset(muts_GO_sig_DS9, GOtype == "bp")
muts_GO_sig_DS9_MF <- subset(muts_GO_sig_DS9, GOtype == "mf")

## Genomic --> transcriptomic connection 
## Will be 10 comparisons - 5 groups, 2 ontologies each
### DS3
DS3_BP <- mgoSim(muts_GO_sig_DS3_BP$numGO, sig_GO_DF_DS3_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
DS3_MF <- mgoSim(muts_GO_sig_DS3_MF$numGO, sig_GO_DF_DS3_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

### DS6
DS6_BP <- mgoSim(muts_GO_sig_DS6_BP$numGO, sig_GO_DF_DS6_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
DS6_MF <- mgoSim(muts_GO_sig_DS6_MF$numGO, sig_GO_DF_DS6_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

### DS7
DS7_BP <- mgoSim(muts_GO_sig_DS7_BP$numGO, sig_GO_DF_DS7_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
DS7_MF <- mgoSim(muts_GO_sig_DS7_MF$numGO, sig_GO_DF_DS7_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

### DS8
DS8_BP <- mgoSim(muts_GO_sig_DS8_BP$numGO, sig_GO_DF_DS8_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
DS8_MF <- mgoSim(muts_GO_sig_DS8_MF$numGO, sig_GO_DF_DS8_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

### DS9
DS9_BP <- mgoSim(muts_GO_sig_DS9_BP$numGO, sig_GO_DF_DS9_BP$numGO, semData=hsGO_bp, measure="Wang", combine="BMA")
DS9_MF <- mgoSim(muts_GO_sig_DS9_MF$numGO, sig_GO_DF_DS9_MF$numGO, semData=hsGO_mf, measure="Wang", combine="BMA")

## Compile into common dataframe
semSimDF_sublines <- data.frame(BP = c(DS3_BP, DS6_BP, DS7_BP, DS8_BP, DS9_BP),
                                MF = c(DS3_MF, DS6_MF, DS7_MF, DS8_MF, DS9_MF))
semSimDF_sublines$id <- c('DS3', 'DS6', 'DS7', 'DS8', 'DS9')
ssDF_sublines <- melt(semSimDF_sublines, measure.vars = c("BP", "MF"), id.vars = "id")

## Plot unique GO ontology type
### BP
ssDF_sublines_BP <- melt(semSimDF_sublines, measure.vars = c("BP"), id.vars = "id")
ggplot(ssDF_sublines_BP, aes(x = id, y = value, fill = id)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Subline", y = "GO Semantic Similarity") +
  ylim(0,0.5) +
  geom_hline(aes(yintercept = ssDF_CLV_BP[1,3]), linetype = "dashed", color = "blue") +
  geom_hline(aes(yintercept = ssDF_CLV_BP[2,3]), linetype = "dashed", color = "green") +
  geom_hline(aes(yintercept = ssDF_CLV_BP[3,3]), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold"),
                    name = "none") +
  ggtitle("Biological Process") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("sublines_GtoTconnection_BP.pdf", width = 4, height = 5)

##### Remove DS8 for paper figure
ssDF_sublines_BP_removeDS8 <- subset(ssDF_sublines_BP, id != "DS8")
ggplot(ssDF_sublines_BP_removeDS8, aes(x = id, y = value, fill = id)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Subline", y = "GO Semantic Similarity") +
  ylim(0,0.5) +
  geom_hline(aes(yintercept = ssDF_CLV_BP[1,3]), linetype = "dashed", color = "blue") +
  geom_hline(aes(yintercept = ssDF_CLV_BP[2,3]), linetype = "dashed", color = "green") +
  geom_hline(aes(yintercept = ssDF_CLV_BP[3,3]), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "gold"),
                    name = "none") +
  ggtitle("Biological Process") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("FIG_4I_left.pdf", width = 4, height = 5)

#### Include DS8, but tint colors not DS8 for emphasis
ssDF_sublines_BP$Tint <- c(0.3, 0.3, 0.3, 1, 0.3)
ggplot(ssDF_sublines_BP, aes(x = factor(id, levels = c("DS3", "DS6", "DS7", "DS9", "DS8")), 
                             y = value, fill = id, alpha = Tint)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Subline", y = "GO Semantic Similarity") +
  ylim(0,0.5) + scale_alpha_continuous(range = c(0.3,1)) +
  geom_hline(aes(yintercept = ssDF_CLV_BP[1,3]), linetype = "dashed", color = "blue") +
  geom_hline(aes(yintercept = ssDF_CLV_BP[2,3]), linetype = "dashed", color = "green") +
  geom_hline(aes(yintercept = ssDF_CLV_BP[3,3]), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold"),
                    name = "none") +
  ggtitle("Biological Process") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("FIG_6E_top.pdf", width = 6, height = 3)

### MF 
#### DS6 is removed from analyses because only 1 differentially expressed GO term
#### from the DEGs --> skews the semantic similarity score
ssDF_sublines_MF <- melt(semSimDF_sublines, measure.vars = c("MF"), id.vars = "id")
ggplot(ssDF_sublines_MF, aes(x = id, y = value, fill = id)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Subline", y = "GO Semantic Similarity") +
  ylim(0,0.35) +
  geom_hline(aes(yintercept = ssDF_CLV_MF[1,3]), linetype = "dashed", color = "blue") +
  geom_hline(aes(yintercept = ssDF_CLV_MF[2,3]), linetype = "dashed", color = "green") +
  geom_hline(aes(yintercept = ssDF_CLV_MF[3,3]), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold"),
                    name = "none") +
  ggtitle("Molecular Function") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("sublines_GtoTconnection_MF.pdf", width = 4, height = 5)

### Remove DS8 for paper (and DS6 - see above)
ssDF_sublines_MF_removeDS8_6 <- subset(ssDF_sublines_MF, !id  %in% c("DS6", "DS8"))
ggplot(ssDF_sublines_MF_removeDS8_6, aes(x = id, y = value, fill = id)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Subline", y = "GO Semantic Similarity") +
  ylim(0,0.35) +
  geom_hline(aes(yintercept = ssDF_CLV_MF[1,3]), linetype = "dashed", color = "blue") +
  geom_hline(aes(yintercept = ssDF_CLV_MF[2,3]), linetype = "dashed", color = "green") +
  geom_hline(aes(yintercept = ssDF_CLV_MF[3,3]), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("brown", "darkorchid", "gold"),
                    name = "none") +
  ggtitle("Molecular Function") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("FIG_4I_right.pdf", width = 4, height = 5)

### Include DS8 
ssDF_sublines_MF_removeDS6 <- subset(ssDF_sublines_MF, !id  %in% c("DS6"))
ssDF_sublines_MF_removeDS6$Tint <- c(0.3, 0.3, 1, 0.3)
ggplot(ssDF_sublines_MF_removeDS6, aes(x = factor(id, levels = c("DS3", "DS7", "DS9", "DS8")),
                                       y = value, fill = id, alpha = Tint)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_classic() + labs(x = "Subline", y = "GO Semantic Similarity") +
  ylim(0,0.35) + scale_alpha_continuous(range = c(0.3,1)) +
  geom_hline(aes(yintercept = ssDF_CLV_MF[1,3]), linetype = "dashed", color = "blue") +
  geom_hline(aes(yintercept = ssDF_CLV_MF[2,3]), linetype = "dashed", color = "green") +
  geom_hline(aes(yintercept = ssDF_CLV_MF[3,3]), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("brown", "darkorchid", "seagreen", "gold"),
                    name = "none") +
  ggtitle("Molecular Function") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
    legend.title = element_blank(), axis.title=element_text(size=14)) +
  ggsave("FIG_6E_bottom.pdf", width = 6, height = 3)

### Create legend for CLV lines on subline plots
aa <- ggplot(ssDF_sublines_MF, aes(x = id, y = value, fill = id)) +
  theme_classic() + labs(x = "Subline", y = "GO Semantic Similarity") +
  ylim(0,0.35) +
  geom_hline(aes(yintercept = ssDF_CLV_MF[1,3], color = "PC9-VU"), linetype = "dashed", show.legend = T) +
  geom_hline(aes(yintercept = ssDF_CLV_MF[2,3], color = "PC9-MGH"), linetype = "dashed", show.legend = T) +
  geom_hline(aes(yintercept = ssDF_CLV_MF[3,3], color = "PC9-BR1"), linetype = "dashed", show.legend = T) +
  scale_color_manual(values = c("PC9-VU" = "blue", "PC9-MGH" = "green", "PC9-BR1" = "red")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
    legend.title = element_blank(), axis.title=element_text(size=12)) 

### Plot legend
library(ggpubr)
leg <- get_legend(aa)
as_ggplot(leg) + ggsave("legend_FIG_4I_6E.pdf", width = 3, height = 1)