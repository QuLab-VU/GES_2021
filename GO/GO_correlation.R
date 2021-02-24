# Load libraries
library(stringr)
library(enrichR)
library(GOSemSim)
library(reshape2)
library(ggplot2)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Load unique variants and differentially expressed genes
# for each cohort
setwd('~/git/GES_2020/GO/')
load('~/git/GES_2020/GO/mutations_DEGs-hg38.RData')
source("~/git/GES_2020/GO/SummarySE.R")

# Load databases
hsGO_bp <- godata('org.Hs.eg.db', ont="BP")
hsGO_mf <- godata('org.Hs.eg.db', ont="MF")
hsGO_cc <- godata('org.Hs.eg.db', ont="CC")

# Use these enrichr GO databases
dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018",
         "GO_Cellular_Component_2018")


# Function to find enriched gene results - experimental
enrichResults_exp <- function(geneList) {
  GO_list <- list()
  i = 1
  for (s in geneList) {
    exp_genes <- s
    enriched <- enrichr(exp_genes, dbs)
    bp <- enriched[["GO_Biological_Process_2018"]]
    mf <- enriched[["GO_Molecular_Function_2018"]]
    cc <- enriched[["GO_Cellular_Component_2018"]]
    bp$Sample = names(geneList)[i]
    mf$Sample = names(geneList)[i]
    cc$Sample = names(geneList)[i]
    bp$GOtype <- "bp"
    mf$GOtype <- "mf"
    cc$GOtype <- "cc"
    GO_samp <- rbind(bp, mf, cc)
    GO_list[[names(geneList)[i]]] <- GO_samp
    i = i + 1
  }
  all_GO_list = do.call(rbind, GO_list)
  all_GO_list$logp <- -log10(all_GO_list$P.value)
  all_GO_list$logq <- -log10(all_GO_list$Adjusted.P.value)
  all_GO_list$numGO <- str_extract(all_GO_list$Term, "GO:[0-9]{1,}")
  all_GO_list
}

## DEGs - Cell Line Versions
sig_GO_DF_exp_RNA_VU <- subset(enrichResults_exp(geneList = samples_CLV), 
                               Sample == "VU")
sig_GO_DF_exp_RNA_VU_BP <- subset(sig_GO_DF_exp_RNA_VU, GOtype == "bp")
sig_GO_DF_exp_RNA_VU_MF <- subset(sig_GO_DF_exp_RNA_VU, GOtype == "mf")
sig_GO_DF_exp_RNA_VU_CC <- subset(sig_GO_DF_exp_RNA_VU, GOtype == "cc")

sig_GO_DF_exp_RNA_MGH <- subset(enrichResults_exp(geneList = samples_CLV), 
                                Sample == "MGH")
sig_GO_DF_exp_RNA_MGH_BP <- subset(sig_GO_DF_exp_RNA_MGH, GOtype == "bp")
sig_GO_DF_exp_RNA_MGH_MF <- subset(sig_GO_DF_exp_RNA_MGH, GOtype == "mf")
sig_GO_DF_exp_RNA_MGH_CC <- subset(sig_GO_DF_exp_RNA_MGH, GOtype == "cc")

sig_GO_DF_exp_RNA_BR1 <- subset(enrichResults_exp(geneList = samples_CLV), 
                                Sample == "BR1")
sig_GO_DF_exp_RNA_BR1_BP <- subset(sig_GO_DF_exp_RNA_BR1, GOtype == "bp")
sig_GO_DF_exp_RNA_BR1_MF <- subset(sig_GO_DF_exp_RNA_BR1, GOtype == "mf")
sig_GO_DF_exp_RNA_BR1_CC <- subset(sig_GO_DF_exp_RNA_BR1, GOtype == "cc")

## Mutations - Cell Line Versions
sig_GO_DF_exp_WXS_VU <- subset(enrichResults_exp(geneList = samples_muts_CLV), 
                               Sample == "VU")
sig_GO_DF_exp_WXS_VU_BP <- subset(sig_GO_DF_exp_WXS_VU, GOtype == "bp")
sig_GO_DF_exp_WXS_VU_MF <- subset(sig_GO_DF_exp_WXS_VU, GOtype == "mf")
sig_GO_DF_exp_WXS_VU_CC <- subset(sig_GO_DF_exp_WXS_VU, GOtype == "cc")

sig_GO_DF_exp_WXS_MGH <- subset(enrichResults_exp(geneList = samples_muts_CLV), 
                                Sample == "MGH")
sig_GO_DF_exp_WXS_MGH_BP <- subset(sig_GO_DF_exp_WXS_MGH, GOtype == "bp")
sig_GO_DF_exp_WXS_MGH_MF <- subset(sig_GO_DF_exp_WXS_MGH, GOtype == "mf")
sig_GO_DF_exp_WXS_MGH_CC <- subset(sig_GO_DF_exp_WXS_MGH, GOtype == "cc")

sig_GO_DF_exp_WXS_BR1 <- subset(enrichResults_exp(geneList = samples_muts_CLV), 
                                Sample == "BR1")
sig_GO_DF_exp_WXS_BR1_BP <- subset(sig_GO_DF_exp_WXS_BR1, GOtype == "bp")
sig_GO_DF_exp_WXS_BR1_MF <- subset(sig_GO_DF_exp_WXS_BR1, GOtype == "mf")
sig_GO_DF_exp_WXS_BR1_CC <- subset(sig_GO_DF_exp_WXS_BR1, GOtype == "cc")

## DEGs - Sublines
sig_GO_DF_exp_RNA_DS3 <- subset(enrichResults_exp(geneList = samples_sublines), 
                                Sample == "DS3")
sig_GO_DF_exp_RNA_DS3_BP <- subset(sig_GO_DF_exp_RNA_DS3, GOtype == "bp")
sig_GO_DF_exp_RNA_DS3_MF <- subset(sig_GO_DF_exp_RNA_DS3, GOtype == "mf")
sig_GO_DF_exp_RNA_DS3_CC <- subset(sig_GO_DF_exp_RNA_DS3, GOtype == "cc")

sig_GO_DF_exp_RNA_DS6 <- subset(enrichResults_exp(geneList = samples_sublines), 
                                Sample == "DS6")
sig_GO_DF_exp_RNA_DS6_BP <- subset(sig_GO_DF_exp_RNA_DS6, GOtype == "bp")
sig_GO_DF_exp_RNA_DS6_MF <- subset(sig_GO_DF_exp_RNA_DS6, GOtype == "mf")
sig_GO_DF_exp_RNA_DS6_CC <- subset(sig_GO_DF_exp_RNA_DS6, GOtype == "cc")

sig_GO_DF_exp_RNA_DS7 <- subset(enrichResults_exp(geneList = samples_sublines), 
                                Sample == "DS7")
sig_GO_DF_exp_RNA_DS7_BP <- subset(sig_GO_DF_exp_RNA_DS7, GOtype == "bp")
sig_GO_DF_exp_RNA_DS7_MF <- subset(sig_GO_DF_exp_RNA_DS7, GOtype == "mf")
sig_GO_DF_exp_RNA_DS7_CC <- subset(sig_GO_DF_exp_RNA_DS7, GOtype == "cc")

sig_GO_DF_exp_RNA_DS8 <- subset(enrichResults_exp(geneList = samples_sublines), 
                                Sample == "DS8")
sig_GO_DF_exp_RNA_DS8_BP <- subset(sig_GO_DF_exp_RNA_DS8, GOtype == "bp")
sig_GO_DF_exp_RNA_DS8_MF <- subset(sig_GO_DF_exp_RNA_DS8, GOtype == "mf")
sig_GO_DF_exp_RNA_DS8_CC <- subset(sig_GO_DF_exp_RNA_DS8, GOtype == "cc")

sig_GO_DF_exp_RNA_DS9 <- subset(enrichResults_exp(geneList = samples_sublines), 
                                Sample == "DS9")
sig_GO_DF_exp_RNA_DS9_BP <- subset(sig_GO_DF_exp_RNA_DS9, GOtype == "bp")
sig_GO_DF_exp_RNA_DS9_MF <- subset(sig_GO_DF_exp_RNA_DS9, GOtype == "mf")
sig_GO_DF_exp_RNA_DS9_CC <- subset(sig_GO_DF_exp_RNA_DS9, GOtype == "cc")

## Mutations - Sublines
sig_GO_DF_exp_WXS_DS3 <- subset(enrichResults_exp(geneList = samples_muts_sublines), 
                                Sample == "DS3")
sig_GO_DF_exp_WXS_DS3_BP <- subset(sig_GO_DF_exp_WXS_DS3, GOtype == "bp")
sig_GO_DF_exp_WXS_DS3_MF <- subset(sig_GO_DF_exp_WXS_DS3, GOtype == "mf")
sig_GO_DF_exp_WXS_DS3_CC <- subset(sig_GO_DF_exp_WXS_DS3, GOtype == "cc")

sig_GO_DF_exp_WXS_DS6 <- subset(enrichResults_exp(geneList = samples_muts_sublines), 
                                Sample == "DS6")
sig_GO_DF_exp_WXS_DS6_BP <- subset(sig_GO_DF_exp_WXS_DS6, GOtype == "bp")
sig_GO_DF_exp_WXS_DS6_MF <- subset(sig_GO_DF_exp_WXS_DS6, GOtype == "mf")
sig_GO_DF_exp_WXS_DS6_CC <- subset(sig_GO_DF_exp_WXS_DS6, GOtype == "cc")

sig_GO_DF_exp_WXS_DS7 <- subset(enrichResults_exp(geneList = samples_muts_sublines), 
                                Sample == "DS7")
sig_GO_DF_exp_WXS_DS7_BP <- subset(sig_GO_DF_exp_WXS_DS7, GOtype == "bp")
sig_GO_DF_exp_WXS_DS7_MF <- subset(sig_GO_DF_exp_WXS_DS7, GOtype == "mf")
sig_GO_DF_exp_WXS_DS7_CC <- subset(sig_GO_DF_exp_WXS_DS7, GOtype == "cc")

sig_GO_DF_exp_WXS_DS8 <- subset(enrichResults_exp(geneList = samples_muts_sublines), 
                                Sample == "DS8")
sig_GO_DF_exp_WXS_DS8_BP <- subset(sig_GO_DF_exp_WXS_DS8, GOtype == "bp")
sig_GO_DF_exp_WXS_DS8_MF <- subset(sig_GO_DF_exp_WXS_DS8, GOtype == "mf")
sig_GO_DF_exp_WXS_DS8_CC <- subset(sig_GO_DF_exp_WXS_DS8, GOtype == "cc")

sig_GO_DF_exp_WXS_DS9 <- subset(enrichResults_exp(geneList = samples_muts_sublines), 
                                Sample == "DS9")
sig_GO_DF_exp_WXS_DS9_BP <- subset(sig_GO_DF_exp_WXS_DS9, GOtype == "bp")
sig_GO_DF_exp_WXS_DS9_MF <- subset(sig_GO_DF_exp_WXS_DS9, GOtype == "mf")
sig_GO_DF_exp_WXS_DS9_CC <- subset(sig_GO_DF_exp_WXS_DS9, GOtype == "cc")

# Build the plots
## BR1
WXS_BR1 <- subset(sig_GO_DF_exp_WXS_BR1, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_BR1 <- subset(sig_GO_DF_exp_RNA_BR1, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_BR1 <- merge(x=WXS_BR1, y=RNA_BR1, by="Term", all.x = T, all.y = T)
merged_BR1_complete <- merged_BR1[complete.cases(merged_BR1),]
merged_BR1_complete$sort <- merged_BR1_complete$logp.x + merged_BR1_complete$logp.y
merged_BR1_complete_sorted <- merged_BR1_complete[order(-merged_BR1_complete$sort),]

## MGH
WXS_MGH <- subset(sig_GO_DF_exp_WXS_MGH, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_MGH <- subset(sig_GO_DF_exp_RNA_MGH, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_MGH <- merge(x=WXS_MGH, y=RNA_MGH, by="Term", all.x = T, all.y = T)
merged_MGH_complete <- merged_MGH[complete.cases(merged_MGH),]
merged_MGH_complete$sort <- merged_MGH_complete$logp.x + merged_MGH_complete$logp.y
merged_MGH_complete_sorted <- merged_MGH_complete[order(-merged_MGH_complete$sort),]

## VU
WXS_VU <- subset(sig_GO_DF_exp_WXS_VU, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_VU <- subset(sig_GO_DF_exp_RNA_VU, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_VU <- merge(x=WXS_VU, y=RNA_VU, by="Term", all.x = T, all.y = T)
merged_VU_complete <- merged_VU[complete.cases(merged_VU),]
merged_VU_complete$sort <- merged_VU_complete$logp.x + merged_VU_complete$logp.y
merged_VU_complete_sorted <- merged_VU_complete[order(-merged_VU_complete$sort),]

## DS3
WXS_DS3 <- subset(sig_GO_DF_exp_WXS_DS3, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_DS3 <- subset(sig_GO_DF_exp_RNA_DS3, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_DS3 <- merge(x=WXS_DS3, y=RNA_DS3, by="Term", all.x = T, all.y = T)
merged_DS3_complete <- merged_DS3[complete.cases(merged_DS3),]
merged_DS3_complete$sort <- merged_DS3_complete$logp.x + merged_DS3_complete$logp.y
merged_DS3_complete_sorted <- merged_DS3_complete[order(-merged_DS3_complete$sort),]

## DS6
WXS_DS6 <- subset(sig_GO_DF_exp_WXS_DS6, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_DS6 <- subset(sig_GO_DF_exp_RNA_DS6, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_DS6 <- merge(x=WXS_DS6, y=RNA_DS6, by="Term", all.x = T, all.y = T)
merged_DS6_complete <- merged_DS6[complete.cases(merged_DS6),]
merged_DS6_complete$sort <- merged_DS6_complete$logp.x + merged_DS6_complete$logp.y
merged_DS6_complete_sorted <- merged_DS6_complete[order(-merged_DS6_complete$sort),]

## DS7
WXS_DS7 <- subset(sig_GO_DF_exp_WXS_DS7, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_DS7 <- subset(sig_GO_DF_exp_RNA_DS7, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_DS7 <- merge(x=WXS_DS7, y=RNA_DS7, by="Term", all.x = T, all.y = T)
merged_DS7_complete <- merged_DS7[complete.cases(merged_DS7),]
merged_DS7_complete$sort <- merged_DS7_complete$logp.x + merged_DS7_complete$logp.y
merged_DS7_complete_sorted <- merged_DS7_complete[order(-merged_DS7_complete$sort),]

## DS8
WXS_DS8 <- subset(sig_GO_DF_exp_WXS_DS8, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_DS8 <- subset(sig_GO_DF_exp_RNA_DS8, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_DS8 <- merge(x=WXS_DS8, y=RNA_DS8, by="Term", all.x = T, all.y = T)
merged_DS8_complete <- merged_DS8[complete.cases(merged_DS8),]
merged_DS8_complete$sort <- merged_DS8_complete$logp.x + merged_DS8_complete$logp.y
merged_DS8_complete_sorted <- merged_DS8_complete[order(-merged_DS8_complete$sort),]

## DS9
WXS_DS9 <- subset(sig_GO_DF_exp_WXS_DS9, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
RNA_DS9 <- subset(sig_GO_DF_exp_RNA_DS9, P.value < 0.05)[,c("numGO", "Term", "Overlap", "logp")]
merged_DS9 <- merge(x=WXS_DS9, y=RNA_DS9, by="Term", all.x = T, all.y = T)
merged_DS9_complete <- merged_DS9[complete.cases(merged_DS9),]
merged_DS9_complete$sort <- merged_DS9_complete$logp.x + merged_DS9_complete$logp.y
merged_DS9_complete_sorted <- merged_DS9_complete[order(-merged_DS9_complete$sort),]

# Joint Plot
## Merge dataframes
merged_BR1_complete_sorted$Population <- "PC9-BR1"
merged_MGH_complete_sorted$Population <- "PC9-MGH"
merged_VU_complete_sorted$Population <- "PC9-VU"
merged_DS3_complete_sorted$Population <- "DS3"
merged_DS6_complete_sorted$Population <- "DS6"
merged_DS7_complete_sorted$Population <- "DS7"
merged_DS8_complete_sorted$Population <- "DS8"
merged_DS9_complete_sorted$Population <- "DS9"

merged_all_complete <- rbind(merged_BR1_complete_sorted, merged_MGH_complete_sorted,
                             merged_VU_complete_sorted, merged_DS3_complete_sorted,
                             merged_DS6_complete_sorted, merged_DS7_complete_sorted,
                             merged_DS8_complete_sorted, merged_DS9_complete_sorted)

## Removal of outlier GO terms that may strongly influence correlations
### 2 for PC9-BR1, 1 for PC9-MGH
merged_all_complete_sub <- subset(merged_all_complete, logp.y < 10)
merged_all_complete_sub$Term <- str_replace(merged_all_complete_sub$Term, " \\s*\\([^\\)]+\\)", "")

## Plot Cell Line Versions together
ggscatter(subset(merged_all_complete_sub, Population %in% c("PC9-VU", "PC9-MGH", "PC9-BR1")), 
          x = "logp.x", y = "logp.y", color = "Population",
          add = "reg.line", size = 1) + #, label = "Term", repel = "True") +
  # geom_label_repel(merged_all_complete_sub, label = Term, size = 4) + 
  # xlim(1.3, 2.5) + ylim(0, 6) +
  stat_cor(method = "spearman", aes(label = ..r.label..), label.x = 3.5, label.y = 4.75, size = 7) +
  ggrepel::geom_label_repel(data = subset(subset(merged_all_complete_sub, 
                                                 Population %in% c("PC9-VU", "PC9-MGH", "PC9-BR1")),
                                          logp.x > 2 | logp.y > 2),
                            aes(label = Term), size = 5) +
  facet_wrap(~Population, ncol = 1) +
  theme_bw() +
  scale_color_manual(values = c("red", "green", "blue")) +
  scale_fill_manual(values = c("red", "green", "blue")) +
  labs(x=expression("-"~Log[10]~"(p)"~GO~Terms[Genomics]),
       y=expression("-"~Log[10]~"(p)"~GO~Terms[Transcriptomics])) +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25),
        strip.text.x = element_text(size=25), strip.text.y = element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  ggsave("FIG_3G.svg", width = 12, height = 12.5)

## Plot Sublines together
### Removed correlations for DS3/DS6/DS7 because they have <=3 (normally 8-10 terms needed)
ggscatter(subset(merged_all_complete_sub, Population %in% c("DS3", "DS6", "DS7", "DS8", "DS9")), 
          x = "logp.x", y = "logp.y", color = "Population",
          add = "reg.line", size = 1) + #, label = "Term", repel = "True") +
  # geom_label_repel(merged_all_complete_sub, label = Term, size = 4) + 
  # xlim(1.3, 4.5) + ylim(0, 6) +
  stat_cor(method = "spearman", aes(label = ..r.label..), label.x = 3, label.y = 2.75, size = 7) +
  # stat_regline_equation(label.x = 3, label.y = 5.1, size = 7) +
  ggrepel::geom_label_repel(data = subset(subset(merged_all_complete_sub, 
                                                 Population %in% c("DS3", "DS6", "DS7", "DS8", "DS9")),
                                          logp.x > 2 | logp.y > 2),
                            aes(label = Term), size = 5) +
  facet_wrap(~Population, ncol = 1) +
  theme_bw() +
  scale_color_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold")) +
  scale_fill_manual(values = c("brown", "deeppink", "darkorchid", "seagreen", "gold")) +
  labs(x=expression("-"~Log[10]~"(p)"~GO~Terms[Genomics]),
       y=expression("-"~Log[10]~"(p)"~GO~Terms[Transcriptomics])) +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25),
        strip.text.x = element_text(size=25), strip.text.y = element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  ggsave("FIG_4G.svg", width = 12, height = 14)

## Calculate sharing across samples - per reviewer clarification
crossprod(table(merged_all_complete_sub[,c("Term","Population")]))