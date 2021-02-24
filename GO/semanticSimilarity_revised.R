# Load libraries
library(stringr)
library(enrichR)
library(GOSemSim)
library(reshape2)
library(ggplot2)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)

# Load unique variants and differentially expressed genes
# for each cohort
setwd('~/git/GES_2020/GO')
load('mutations_DEGs-hg38.RData')
source("SummarySE.R")

# Load databases
hsGO_bp <- godata('org.Hs.eg.db', ont="BP")
hsGO_mf <- godata('org.Hs.eg.db', ont="MF")
hsGO_cc <- godata('org.Hs.eg.db', ont="CC")

# Pull all possible genes from hg38 - download from ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/
txdb <- makeTxDbFromGFF(file = "Homo_sapiens.GRCh38.98.gtf", format = "gtf")
genes <- genes(txdb)
genes_symbol <- unname(mapIds(org.Hs.eg.db,keys=genes$gene_id,
                              column="SYMBOL",keytype="ENSEMBL",
                              multiVals="first"))
genes_symbol <- genes_symbol[!is.na(genes_symbol)]

# Use these enrichr GO databases
dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018",
         "GO_Cellular_Component_2018")

# Functions
## GO CI Method
CI_method <- function(dat) {
  ls_all = list()
  for (vs in seq(100)) {
    ls_keep = list()
    for (r in seq(nrow(dat))) {
      # print(r)
      val <- 1-dat$P.value[r]
      rn <- runif(n=1, min=1e-12, max=.9999999999)
      if (rn < val) {
        ls_keep = append(ls_keep, dat$numGO[r])
      }
    }
    ls_all = append(ls_all, list(ls_keep))
  }
  ls_df <- do.call("cbind", ls_all)
  ls_df
}

## Perform semantic similarity on resampled p-values
GSS_out <- function(d1, d2, hs_type) {
  ls_cxn = list()
  for (i in seq(10)) {
    # for (i in seq(100)) {
    print(i)
    # cxn_sub <- mgoSim(unlist(CI_method(dat = d1)[,i]),
    #                   unlist(CI_method(dat = d2)[,i]),
    #                   semData = hs_type, measure = "Wang", combine = "BMA")
    cxn_sub <- goSim(unlist(CI_method(dat = d1)[,i]),
                     unlist(CI_method(dat = d2)[,i]),
                     semData = hs_type, measure = "Wang")
    ls_cxn[[i]] <- cxn_sub
  }
  ls_cxn
}

## Function to find enriched gene results - experimental
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

## Function to find enriched gene results - simulated (sampling random genes)
enrichResults_sim <- function(geneList) {
  GO_list <- list()
  i = 1
  for (s in geneList) {
    sim_genes <- sample(x = genes_symbol, size = length(s), replace = FALSE)
    enriched <- enrichr(sim_genes, dbs)
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

# Experimental
## Quantify simulated CLV
### DEGs
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

### Mutations
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

## Quantify expulated sublines
### DEGs
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


### Mutations
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


# Simulated
## Quantify simulated CLV
### DEGs
sig_GO_DF_sim_RNA_VU <- subset(enrichResults_sim(geneList = samples_CLV),
                               Sample == "VU")
sig_GO_DF_sim_RNA_VU_BP <- subset(sig_GO_DF_sim_RNA_VU, GOtype == "bp")
sig_GO_DF_sim_RNA_VU_MF <- subset(sig_GO_DF_sim_RNA_VU, GOtype == "mf")
sig_GO_DF_sim_RNA_VU_CC <- subset(sig_GO_DF_sim_RNA_VU, GOtype == "cc")

sig_GO_DF_sim_RNA_MGH <- subset(enrichResults_sim(geneList = samples_CLV),
                                Sample == "MGH")
sig_GO_DF_sim_RNA_MGH_BP <- subset(sig_GO_DF_sim_RNA_MGH, GOtype == "bp")
sig_GO_DF_sim_RNA_MGH_MF <- subset(sig_GO_DF_sim_RNA_MGH, GOtype == "mf")
sig_GO_DF_sim_RNA_MGH_CC <- subset(sig_GO_DF_sim_RNA_MGH, GOtype == "cc")

sig_GO_DF_sim_RNA_BR1 <- subset(enrichResults_sim(geneList = samples_CLV),
                                Sample == "BR1")
sig_GO_DF_sim_RNA_BR1_BP <- subset(sig_GO_DF_sim_RNA_BR1, GOtype == "bp")
sig_GO_DF_sim_RNA_BR1_MF <- subset(sig_GO_DF_sim_RNA_BR1, GOtype == "mf")
sig_GO_DF_sim_RNA_BR1_CC <- subset(sig_GO_DF_sim_RNA_BR1, GOtype == "cc")

### Mutations
sig_GO_DF_sim_WXS_VU <- subset(enrichResults_sim(geneList = samples_muts_CLV),
                               Sample == "VU")
sig_GO_DF_sim_WXS_VU_BP <- subset(sig_GO_DF_sim_WXS_VU, GOtype == "bp")
sig_GO_DF_sim_WXS_VU_MF <- subset(sig_GO_DF_sim_WXS_VU, GOtype == "mf")
sig_GO_DF_sim_WXS_VU_CC <- subset(sig_GO_DF_sim_WXS_VU, GOtype == "cc")

sig_GO_DF_sim_WXS_MGH <- subset(enrichResults_sim(geneList = samples_muts_CLV),
                                Sample == "MGH")
sig_GO_DF_sim_WXS_MGH_BP <- subset(sig_GO_DF_sim_WXS_MGH, GOtype == "bp")
sig_GO_DF_sim_WXS_MGH_MF <- subset(sig_GO_DF_sim_WXS_MGH, GOtype == "mf")
sig_GO_DF_sim_WXS_MGH_CC <- subset(sig_GO_DF_sim_WXS_MGH, GOtype == "cc")

sig_GO_DF_sim_WXS_BR1 <- subset(enrichResults_sim(geneList = samples_muts_CLV),
                                Sample == "BR1")
sig_GO_DF_sim_WXS_BR1_BP <- subset(sig_GO_DF_sim_WXS_BR1, GOtype == "bp")
sig_GO_DF_sim_WXS_BR1_MF <- subset(sig_GO_DF_sim_WXS_BR1, GOtype == "mf")
sig_GO_DF_sim_WXS_BR1_CC <- subset(sig_GO_DF_sim_WXS_BR1, GOtype == "cc")

## Quantify simulated sublines
### DEGs
sig_GO_DF_sim_RNA_DS3 <- subset(enrichResults_sim(geneList = samples_sublines),
                                Sample == "DS3")
sig_GO_DF_sim_RNA_DS3_BP <- subset(sig_GO_DF_sim_RNA_DS3, GOtype == "bp")
sig_GO_DF_sim_RNA_DS3_MF <- subset(sig_GO_DF_sim_RNA_DS3, GOtype == "mf")
sig_GO_DF_sim_RNA_DS3_CC <- subset(sig_GO_DF_sim_RNA_DS3, GOtype == "cc")

sig_GO_DF_sim_RNA_DS6 <- subset(enrichResults_sim(geneList = samples_sublines),
                                Sample == "DS6")
sig_GO_DF_sim_RNA_DS6_BP <- subset(sig_GO_DF_sim_RNA_DS6, GOtype == "bp")
sig_GO_DF_sim_RNA_DS6_MF <- subset(sig_GO_DF_sim_RNA_DS6, GOtype == "mf")
sig_GO_DF_sim_RNA_DS6_CC <- subset(sig_GO_DF_sim_RNA_DS6, GOtype == "cc")

sig_GO_DF_sim_RNA_DS7 <- subset(enrichResults_sim(geneList = samples_sublines),
                                Sample == "DS7")
sig_GO_DF_sim_RNA_DS7_BP <- subset(sig_GO_DF_sim_RNA_DS7, GOtype == "bp")
sig_GO_DF_sim_RNA_DS7_MF <- subset(sig_GO_DF_sim_RNA_DS7, GOtype == "mf")
sig_GO_DF_sim_RNA_DS7_CC <- subset(sig_GO_DF_sim_RNA_DS7, GOtype == "cc")

sig_GO_DF_sim_RNA_DS8 <- subset(enrichResults_sim(geneList = samples_sublines),
                                Sample == "DS8")
sig_GO_DF_sim_RNA_DS8_BP <- subset(sig_GO_DF_sim_RNA_DS8, GOtype == "bp")
sig_GO_DF_sim_RNA_DS8_MF <- subset(sig_GO_DF_sim_RNA_DS8, GOtype == "mf")
sig_GO_DF_sim_RNA_DS8_CC <- subset(sig_GO_DF_sim_RNA_DS8, GOtype == "cc")

sig_GO_DF_sim_RNA_DS9 <- subset(enrichResults_sim(geneList = samples_sublines),
                                Sample == "DS9")
sig_GO_DF_sim_RNA_DS9_BP <- subset(sig_GO_DF_sim_RNA_DS9, GOtype == "bp")
sig_GO_DF_sim_RNA_DS9_MF <- subset(sig_GO_DF_sim_RNA_DS9, GOtype == "mf")
sig_GO_DF_sim_RNA_DS9_CC <- subset(sig_GO_DF_sim_RNA_DS9, GOtype == "cc")

### Mutations
sig_GO_DF_sim_WXS_DS3 <- subset(enrichResults_sim(geneList = samples_muts_sublines),
                                Sample == "DS3")
sig_GO_DF_sim_WXS_DS3_BP <- subset(sig_GO_DF_sim_WXS_DS3, GOtype == "bp")
sig_GO_DF_sim_WXS_DS3_MF <- subset(sig_GO_DF_sim_WXS_DS3, GOtype == "mf")
sig_GO_DF_sim_WXS_DS3_CC <- subset(sig_GO_DF_sim_WXS_DS3, GOtype == "cc")

sig_GO_DF_sim_WXS_DS6 <- subset(enrichResults_sim(geneList = samples_muts_sublines),
                                Sample == "DS6")
sig_GO_DF_sim_WXS_DS6_BP <- subset(sig_GO_DF_sim_WXS_DS6, GOtype == "bp")
sig_GO_DF_sim_WXS_DS6_MF <- subset(sig_GO_DF_sim_WXS_DS6, GOtype == "mf")
sig_GO_DF_sim_WXS_DS6_CC <- subset(sig_GO_DF_sim_WXS_DS6, GOtype == "cc")

sig_GO_DF_sim_WXS_DS7 <- subset(enrichResults_sim(geneList = samples_muts_sublines),
                                Sample == "DS7")
sig_GO_DF_sim_WXS_DS7_BP <- subset(sig_GO_DF_sim_WXS_DS7, GOtype == "bp")
sig_GO_DF_sim_WXS_DS7_MF <- subset(sig_GO_DF_sim_WXS_DS7, GOtype == "mf")
sig_GO_DF_sim_WXS_DS7_CC <- subset(sig_GO_DF_sim_WXS_DS7, GOtype == "cc")

sig_GO_DF_sim_WXS_DS8 <- subset(enrichResults_sim(geneList = samples_muts_sublines),
                                Sample == "DS8")
sig_GO_DF_sim_WXS_DS8_BP <- subset(sig_GO_DF_sim_WXS_DS8, GOtype == "bp")
sig_GO_DF_sim_WXS_DS8_MF <- subset(sig_GO_DF_sim_WXS_DS8, GOtype == "mf")
sig_GO_DF_sim_WXS_DS8_CC <- subset(sig_GO_DF_sim_WXS_DS8, GOtype == "cc")

sig_GO_DF_sim_WXS_DS9 <- subset(enrichResults_sim(geneList = samples_muts_sublines),
                                Sample == "DS9")
sig_GO_DF_sim_WXS_DS9_BP <- subset(sig_GO_DF_sim_WXS_DS9, GOtype == "bp")
sig_GO_DF_sim_WXS_DS9_MF <- subset(sig_GO_DF_sim_WXS_DS9, GOtype == "mf")
sig_GO_DF_sim_WXS_DS9_CC <- subset(sig_GO_DF_sim_WXS_DS9, GOtype == "cc")


# Calculate Semantic similarity
## CC
### Experimental
BR1_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_BR1_CC, d2 = sig_GO_DF_exp_RNA_BR1_CC, hs_type = hsGO_cc)
MGH_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_MGH_CC, d2 = sig_GO_DF_exp_RNA_MGH_CC, hs_type = hsGO_cc)
VU_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_VU_CC, d2 = sig_GO_DF_exp_RNA_VU_CC, hs_type = hsGO_cc)
DS3_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS3_CC, d2 = sig_GO_DF_exp_RNA_DS3_CC, hs_type = hsGO_cc)
DS6_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS6_CC, d2 = sig_GO_DF_exp_RNA_DS6_CC, hs_type = hsGO_cc)
DS7_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS7_CC, d2 = sig_GO_DF_exp_RNA_DS7_CC, hs_type = hsGO_cc)
DS8_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS8_CC, d2 = sig_GO_DF_exp_RNA_DS8_CC, hs_type = hsGO_cc)
DS9_CC_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS9_CC, d2 = sig_GO_DF_exp_RNA_DS9_CC, hs_type = hsGO_cc)

### Simulated
BR1_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_BR1_CC, d2 = sig_GO_DF_sim_RNA_BR1_CC, hs_type = hsGO_cc)
MGH_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_MGH_CC, d2 = sig_GO_DF_sim_RNA_MGH_CC, hs_type = hsGO_cc)
VU_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_VU_CC, d2 = sig_GO_DF_sim_RNA_VU_CC, hs_type = hsGO_cc)
DS3_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS3_CC, d2 = sig_GO_DF_sim_RNA_DS3_CC, hs_type = hsGO_cc)
DS6_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS6_CC, d2 = sig_GO_DF_sim_RNA_DS6_CC, hs_type = hsGO_cc)
DS7_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS7_CC, d2 = sig_GO_DF_sim_RNA_DS7_CC, hs_type = hsGO_cc)
DS8_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS8_CC, d2 = sig_GO_DF_sim_RNA_DS8_CC, hs_type = hsGO_cc)
DS9_CC_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS9_CC, d2 = sig_GO_DF_sim_RNA_DS9_CC, hs_type = hsGO_cc)

## MF
### Experimental
BR1_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_BR1_MF, d2 = sig_GO_DF_exp_RNA_BR1_MF, hs_type = hsGO_mf)
MGH_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_MGH_MF, d2 = sig_GO_DF_exp_RNA_MGH_MF, hs_type = hsGO_mf)
VU_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_VU_MF, d2 = sig_GO_DF_exp_RNA_VU_MF, hs_type = hsGO_mf)
DS3_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS3_MF, d2 = sig_GO_DF_exp_RNA_DS3_MF, hs_type = hsGO_mf)
DS6_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS6_MF, d2 = sig_GO_DF_exp_RNA_DS6_MF, hs_type = hsGO_mf)
DS7_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS7_MF, d2 = sig_GO_DF_exp_RNA_DS7_MF, hs_type = hsGO_mf)
DS8_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS8_MF, d2 = sig_GO_DF_exp_RNA_DS8_MF, hs_type = hsGO_mf)
DS9_MF_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS9_MF, d2 = sig_GO_DF_exp_RNA_DS9_MF, hs_type = hsGO_mf)

### Simulated
BR1_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_BR1_MF, d2 = sig_GO_DF_sim_RNA_BR1_MF, hs_type = hsGO_mf)
MGH_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_MGH_MF, d2 = sig_GO_DF_sim_RNA_MGH_MF, hs_type = hsGO_mf)
VU_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_VU_MF, d2 = sig_GO_DF_sim_RNA_VU_MF, hs_type = hsGO_mf)
DS3_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS3_MF, d2 = sig_GO_DF_sim_RNA_DS3_MF, hs_type = hsGO_mf)
DS6_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS6_MF, d2 = sig_GO_DF_sim_RNA_DS6_MF, hs_type = hsGO_mf)
DS7_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS7_MF, d2 = sig_GO_DF_sim_RNA_DS7_MF, hs_type = hsGO_mf)
DS8_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS8_MF, d2 = sig_GO_DF_sim_RNA_DS8_MF, hs_type = hsGO_mf)
DS9_MF_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS9_MF, d2 = sig_GO_DF_sim_RNA_DS9_MF, hs_type = hsGO_mf)

## BP
### Experimental
BR1_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_BR1_BP, d2 = sig_GO_DF_exp_RNA_BR1_BP, hs_type = hsGO_bp)
MGH_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_MGH_BP, d2 = sig_GO_DF_exp_RNA_MGH_BP, hs_type = hsGO_bp)
VU_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_VU_BP, d2 = sig_GO_DF_exp_RNA_VU_BP, hs_type = hsGO_bp)
DS3_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS3_BP, d2 = sig_GO_DF_exp_RNA_DS3_BP, hs_type = hsGO_bp)
DS6_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS6_BP, d2 = sig_GO_DF_exp_RNA_DS6_BP, hs_type = hsGO_bp)
DS7_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS7_BP, d2 = sig_GO_DF_exp_RNA_DS7_BP, hs_type = hsGO_bp)
DS8_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS8_BP, d2 = sig_GO_DF_exp_RNA_DS8_BP, hs_type = hsGO_bp)
DS9_BP_exp <- GSS_out(d1 = sig_GO_DF_exp_WXS_DS9_BP, d2 = sig_GO_DF_exp_RNA_DS9_BP, hs_type = hsGO_bp)

### Simulated
BR1_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_BR1_BP, d2 = sig_GO_DF_sim_RNA_BR1_BP, hs_type = hsGO_bp)
MGH_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_MGH_BP, d2 = sig_GO_DF_sim_RNA_MGH_BP, hs_type = hsGO_bp)
VU_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_VU_BP, d2 = sig_GO_DF_sim_RNA_VU_BP, hs_type = hsGO_bp)
DS3_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS3_BP, d2 = sig_GO_DF_sim_RNA_DS3_BP, hs_type = hsGO_bp)
DS6_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS6_BP, d2 = sig_GO_DF_sim_RNA_DS6_BP, hs_type = hsGO_bp)
DS7_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS7_BP, d2 = sig_GO_DF_sim_RNA_DS7_BP, hs_type = hsGO_bp)
DS8_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS8_BP, d2 = sig_GO_DF_sim_RNA_DS8_BP, hs_type = hsGO_bp)
DS9_BP_sim <- GSS_out(d1 = sig_GO_DF_sim_WXS_DS9_BP, d2 = sig_GO_DF_sim_RNA_DS9_BP, hs_type = hsGO_bp)

# # Save the data (takes hours to run)
# save(BR1_CC_exp, MGH_CC_exp, VU_CC_exp, DS3_CC_exp, DS6_CC_exp, DS7_CC_exp, DS8_CC_exp, DS9_CC_exp,
#      BR1_MF_exp, MGH_MF_exp, VU_MF_exp, DS3_MF_exp, DS6_MF_exp, DS7_MF_exp, DS8_MF_exp, DS9_MF_exp,
#      BR1_BP_exp, MGH_BP_exp, VU_BP_exp, DS3_BP_exp, DS6_BP_exp, DS7_BP_exp, DS8_BP_exp, DS9_BP_exp,
#      BR1_CC_sim, MGH_CC_sim, VU_CC_sim, DS3_CC_sim, DS6_CC_sim, DS7_CC_sim, DS8_CC_sim, DS9_CC_sim,
#      BR1_MF_sim, MGH_MF_sim, VU_MF_sim, DS3_MF_sim, DS6_MF_sim, DS7_MF_sim, DS8_MF_sim, DS9_MF_sim,
#      BR1_BP_sim, MGH_BP_sim, VU_BP_sim, DS3_BP_sim, DS6_BP_sim, DS7_BP_sim, DS8_BP_sim, DS9_BP_sim,
#      file = "semSim_comparison_hg38_mac_distributionsNoMax.RData")

# Plot simulated data
# all_dat_sim <- list(vu_bp = VU_BP_sim[[1]], mgh_bp = MGH_BP_sim[[1]], br1_bp = BR1_BP_sim[[1]], ds3_bp = DS3_BP_sim[[1]],
#                     ds6_bp = DS6_BP_sim[[1]], ds7_bp = DS7_BP_sim[[1]], ds8_bp = DS8_BP_sim[[1]], ds9_bp = DS9_BP_sim[[1]],
#                     vu_mf = VU_MF_sim[[1]], mgh_mf = MGH_MF_sim[[1]], br1_mf = BR1_MF_sim[[1]], ds3_mf = DS3_MF_sim[[1]],
#                     ds6_mf = DS6_MF_sim[[1]], ds7_mf = DS7_MF_sim[[1]], ds8_mf = DS8_MF_sim[[1]], ds9_mf = DS9_MF_sim[[1]],
#                     vu_cc = VU_CC_sim[[1]], mgh_cc = MGH_CC_sim[[1]], br1_cc = BR1_CC_sim[[1]], ds3_cc = DS3_CC_sim[[1]],
#                     ds6_cc = DS6_CC_sim[[1]], ds7_cc = DS7_CC_sim[[1]], ds8_cc = DS8_CC_sim[[1]], ds9_cc = DS9_CC_sim[[1]])

all_dat_sim <- list(vu_bp = sort(VU_BP_sim[[1]], decreasing = T)[1:1000], mgh_bp = sort(MGH_BP_sim[[1]], decreasing = T)[1:1000],
                    br1_bp = sort(BR1_BP_sim[[1]], decreasing = T)[1:1000], ds3_bp = sort(DS3_BP_sim[[1]], decreasing = T)[1:1000],
                    ds6_bp = sort(DS6_BP_sim[[1]], decreasing = T)[1:1000], ds7_bp = sort(DS7_BP_sim[[1]], decreasing = T)[1:1000],
                    ds8_bp = sort(DS8_BP_sim[[1]], decreasing = T)[1:1000], ds9_bp = sort(DS9_BP_sim[[1]], decreasing = T)[1:1000],
                    vu_mf = sort(VU_MF_sim[[1]], decreasing = T)[1:1000], mgh_mf = sort(MGH_MF_sim[[1]], decreasing = T)[1:1000],
                    br1_mf = sort(BR1_MF_sim[[1]], decreasing = T)[1:1000], ds3_mf = sort(DS3_MF_sim[[1]], decreasing = T)[1:1000],
                    ds6_mf = sort(DS6_MF_sim[[1]], decreasing = T)[1:1000], ds7_mf = sort(DS7_MF_sim[[1]], decreasing = T)[1:1000],
                    ds8_mf = sort(DS8_MF_sim[[1]], decreasing = T)[1:1000], ds9_mf = sort(DS9_MF_sim[[1]], decreasing = T)[1:1000],
                    vu_cc = sort(VU_CC_sim[[1]], decreasing = T)[1:1000], mgh_cc = sort(MGH_CC_sim[[1]], decreasing = T)[1:1000],
                    br1_cc = sort(BR1_CC_sim[[1]], decreasing = T)[1:1000], ds3_cc = sort(DS3_CC_sim[[1]], decreasing = T)[1:1000],
                    ds6_cc = sort(DS6_CC_sim[[1]], decreasing = T)[1:1000], ds7_cc = sort(DS7_CC_sim[[1]], decreasing = T)[1:1000],
                    ds8_cc = sort(DS8_CC_sim[[1]], decreasing = T)[1:1000], ds9_cc = sort(DS9_CC_sim[[1]], decreasing = T)[1:1000])

# Put in common dataframe
all_dat_df_sim <- as.data.frame(stack(all_dat_sim))
all_dat_df_sim <- within(all_dat_df_sim,
                         ind <-data.frame(do.call('rbind',
                                                  strsplit(as.character(toupper(ind)),
                                                           '_', fixed=TRUE))))
all_dat_df_sim <- do.call('data.frame', all_dat_df_sim)
names(all_dat_df_sim) <- c("value", "Population", "GOtype")

all_dat_df_sim$Population <- factor(all_dat_df_sim$Population,
                                    levels = c("VU", "MGH", "BR1", "DS3",
                                               "DS6", "DS7", "DS8", "DS9"))
all_dat_df_sim$Group <- "Simulated"

# Plot experimental data
# all_dat_exp <- list(vu_bp = VU_BP_exp[[1]], mgh_bp = MGH_BP_exp[[1]], br1_bp = BR1_BP_exp[[1]], ds3_bp = DS3_BP_exp[[1]],
#                     ds6_bp = DS6_BP_exp[[1]], ds7_bp = DS7_BP_exp[[1]], ds8_bp = DS8_BP_exp[[1]], ds9_bp = DS9_BP_exp[[1]],
#                     vu_mf = VU_MF_exp[[1]], mgh_mf = MGH_MF_exp[[1]], br1_mf = BR1_MF_exp[[1]], ds3_mf = DS3_MF_exp[[1]],
#                     ds6_mf = DS6_MF_exp[[1]], ds7_mf = DS7_MF_exp[[1]], ds8_mf = DS8_MF_exp[[1]], ds9_mf = DS9_MF_exp[[1]],
#                     vu_cc = VU_CC_exp[[1]], mgh_cc = MGH_CC_exp[[1]], br1_cc = BR1_CC_exp[[1]], ds3_cc = DS3_CC_exp[[1]],
#                     ds6_cc = DS6_CC_exp[[1]], ds7_cc = DS7_CC_exp[[1]], ds8_cc = DS8_CC_exp[[1]], ds9_cc = DS9_CC_exp[[1]])

all_dat_exp <- list(vu_bp = sort(VU_BP_exp[[1]], decreasing = T)[1:1000], mgh_bp = sort(MGH_BP_exp[[1]], decreasing = T)[1:1000],
                    br1_bp = sort(BR1_BP_exp[[1]], decreasing = T)[1:1000], ds3_bp = sort(DS3_BP_exp[[1]], decreasing = T)[1:1000],
                    ds6_bp = sort(DS6_BP_exp[[1]], decreasing = T)[1:1000], ds7_bp = sort(DS7_BP_exp[[1]], decreasing = T)[1:1000],
                    ds8_bp = sort(DS8_BP_exp[[1]], decreasing = T)[1:1000], ds9_bp = sort(DS9_BP_exp[[1]], decreasing = T)[1:1000],
                    vu_mf = sort(VU_MF_exp[[1]], decreasing = T)[1:1000], mgh_mf = sort(MGH_MF_exp[[1]], decreasing = T)[1:1000],
                    br1_mf = sort(BR1_MF_exp[[1]], decreasing = T)[1:1000], ds3_mf = sort(DS3_MF_exp[[1]], decreasing = T)[1:1000],
                    ds6_mf = sort(DS6_MF_exp[[1]], decreasing = T)[1:1000], ds7_mf = sort(DS7_MF_exp[[1]], decreasing = T)[1:1000],
                    ds8_mf = sort(DS8_MF_exp[[1]], decreasing = T)[1:1000], ds9_mf = sort(DS9_MF_exp[[1]], decreasing = T)[1:1000],
                    vu_cc = sort(VU_CC_exp[[1]], decreasing = T)[1:1000], mgh_cc = sort(MGH_CC_exp[[1]], decreasing = T)[1:1000],
                    br1_cc = sort(BR1_CC_exp[[1]], decreasing = T)[1:1000], ds3_cc = sort(DS3_CC_exp[[1]], decreasing = T)[1:1000],
                    ds6_cc = sort(DS6_CC_exp[[1]], decreasing = T)[1:1000], ds7_cc = sort(DS7_CC_exp[[1]], decreasing = T)[1:1000],
                    ds8_cc = sort(DS8_CC_exp[[1]], decreasing = T)[1:1000], ds9_cc = sort(DS9_CC_exp[[1]], decreasing = T)[1:1000])

# Put in common dataframe
all_dat_df_exp <- as.data.frame(stack(all_dat_exp))
all_dat_df_exp <- within(all_dat_df_exp,
                         ind <-data.frame(do.call('rbind',
                                                  strsplit(as.character(toupper(ind)),
                                                           '_', fixed=TRUE))))
all_dat_df_exp <- do.call('data.frame', all_dat_df_exp)
names(all_dat_df_exp) <- c("value", "Population", "GOtype")

all_dat_df_exp$Population <- factor(all_dat_df_exp$Population,
                                    levels = c("VU", "MGH", "BR1", "DS3",
                                               "DS6", "DS7", "DS8", "DS9"))
all_dat_df_exp$Group <- "Experimental"


# Put data together
all_dat_df_all <- rbind(all_dat_df_exp, all_dat_df_sim)
all_dat_df_all$GOtype <- factor(all_dat_df_all$GOtype, levels = c("BP", "MF", "CC"))

# # Plot as boxPlot
# ggplot(all_dat_df_all, aes(x=Population, y=value, fill = Population, linetype = Group)) +
#   geom_boxplot() + facet_grid(GOtype ~ ., scales = "free") + theme_bw() +
#   xlab("Population") + ylab("GO Semantic Similarity") +
#   scale_fill_manual(values = c("blue", "green", "red", "brown",
#                                "deeppink", "darkorchid", "seagreen", "gold")) +
#   theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=14),
#     legend.title = element_blank(), axis.title=element_text(size=14),
#     strip.text = element_text(size = 14), legend.position = "bottom") +
#   ggsave("GOSemSim_simulated_violinNoOutliers_fromPvalues_hg38_distributionsMax1000.pdf", width = 12, height = 15)


## Normalize data
all_dat_df_sum <- summarySE(all_dat_df_all,
                            measurevar = "value",
                            groupvars = c("GOtype", "Population", "Group"))
all_dat_df_sum_exp <- subset(all_dat_df_sum, Group == "Experimental")
all_dat_df_sum_sim <- subset(all_dat_df_sum, Group == "Simulated")
highVal <- all_dat_df_sum_sim$value + all_dat_df_sum_sim$sd
all_dat_df_sum_exp$compVal <- all_dat_df_sum_exp$value - highVal


## Plot as confidence interval relative to simulated baseline (+1SD)
ggplot(all_dat_df_sum_exp, aes(x=Population, y=compVal, color=Population, group=Population)) +
  geom_errorbar(aes(ymin=compVal-ci, ymax=compVal+ci), width=0.2, size = 1) +
  geom_point(size = 1, shape = 21, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(GOtype ~ ., scales = "free") + theme_bw() +
  xlab("Population") + ylab("GO Semantic Similarity (relative to baseline)") +
  scale_color_manual(values = c("blue", "green", "red", "brown",
                                "deeppink", "darkorchid", "seagreen", "gold")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5), axis.text=element_text(size=18),
    legend.title = element_blank(), axis.title=element_text(size=18),
    strip.text = element_text(size = 18), legend.position = "none") +
  ggsave("GOSemSim_medianCI_fromPvalues_hg38_distributionsMax1000.pdf", width = 12, height = 15)