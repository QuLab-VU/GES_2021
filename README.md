# GES_2020
Code (data analysis and model simulations) for GES paper (2020).

clonal Fractional Proliferation (cFP) folder includes data for each 
dataset (one csv = one 384-well plate, each well with single-colony growth
rates). cFP.R analysis script plots figures in text associated with drug-induced
proliferation (DIP) rate distributions.

DrugResponse folder includes drug-response count data, annotated by files in
the platemap folder. Two scripts are included, functions from the previously
published diprate R package (Harris et al., Nat Meth (2016)) and the script
to reproduce drug response figures in the paper.

Joint_functions includes plotting functions shared between the cFP and DrugResponse
folders.

Simulations folder includes the monoclonal and polyclonal growth model (M/PGM) 
simulations and data creation. Plotting happens in cFP folder. In addition, Simulations
also contains code to complete parameter scans of both models (also plotted in 
cFP folder script).

Whole exome sequencing (WES) folder includes summarized metrics from the data, an .xlm
file to identify the canonical ex19del mutation in all cell populations, and code to
pull data from separate sources and create figures in paper.

Single-cell RNA sequencing (scRNAseq) folder includes three data folders, two scripts, 
and some resource files. The read_count folder includes the gene expression library for
all 8 cell populations, which were sequenced together. The umi_count folder includes the 
hashtag oligonucleotide (HTO) library, created by CITE-seq Count function on the 'hashed'
data. The HTO_identification folder provides additional information on the HTO library. 
scRNAseq.R creates figures from the paper, and findDEGs.R identifies lists of differentially
expressed genes (DEGs) to be used in the GO folder. The cell cycle genes file is used for 
a cell cycle regression analysis in the scRNAseq.R script.

The gene ontology (GO) folder has a script that performs a GO, graph-structure based 
semantic similarity analysis on mutations and DEGs that define each cell population, in 
order to find a genetic-to-epigenetic connection between the datasets. Two RData files 
included in the folder are calculated separately in the WES and scRNAseq folders.

RNA sequencing (RNAseq) folder includes a script that runs a simple clustering and visualization
of bulk RNAseq data for all cell populations.


Data can be accessed in the following databases:
	SRA (#######)
	GEO (#######)
