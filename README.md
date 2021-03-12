## Data repository for &quot;An *in vitro* model of tumor heterogeneity resolves genetic, epigenetic, and stochastic sources of cell state variability,&quot; Hayford et al. (2021), *[journal vol : article]; [DOI: xxx](http://dx.doi.org/xxx)*

---

### ***Instructions for creating panels in all main and supplementary figures based on experimental and simulated data in this repository***

- #### <ins>MAIN FIGURES</ins>

	- #### <ins>FIGURE 1</ins>: *N/A*
	
	- #### <ins>FIGURE 2</ins>
		
		**Panels A and C**: In the [DrugResponse](https://github.com/QuLab-VU/GES_2021/tree/master/DrugResponse) directory, run ``DrugResponse.R``, which pulls data from the two `Parental-*.csv` files in the directory and the well conditions in the [DrugResponse/Platemaps](https://github.com/QuLab-VU/GES_2021/tree/master/DrugResponse/Platemaps) subdirectory.
		
		**Panels B and D**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run ``cFP.R``, which pulls data from the 10 `cFP_*.csv` files in the directory.
		
	- #### <ins>FIGURE 3</ins>
		
		**Panel A**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `mutations_byChromosome.csv`.
		
		**Panels B, C, and D**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from the `vep_*.txt` files in the directory and uses the database in the RData object in `RefCDS_human_GRCH38.p12.rda` to cross-reference variants. *NOTE: The `vep_*.txt` files must be manually unzipped before running `WES.R`.*
		
		**Panel E**: In the [scRNAseq/inferCNV](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/inferCNV) subdirectory, run `inferCNV.R`, which pulls a counts matrix from the RData object in `PC9.CLV.10x.counts.matrix.rds`, included in the directory. Necessary annotation and gene order files are also provided.
		
		**Panel F**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories. Scripts to de-multiplex hashed raw data and outputs are included in the [scRNAseq/HTO_identification](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/HTO_identification) subdirectory. A full matrix of de-multiplexed counts is included as `PC9_scRNAseqCounts_HTOdemux.csv.zip`.
		
		**Panel G**: In the [GO](https://github.com/QuLab-VU/GES_2021/tree/master/GO) directory, run `GO_correlation.R`, which pulls data from `mutations_DEGs-hg38.RData`, a file that compiles all IMPACT genetic mutations (from the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory) and differentially expressed genes (DEGs; from the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory).
		
	- #### <ins>FIGURE 4</ins>
		
		**Panel A**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) folder, run `WES.R`, which pulls data from `mutations_byChromosome.csv`.
		
		**Panels B, C, and D**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from the `vep_*.txt` files in the directory and uses the database in the RData object in `RefCDS_human_GRCH38.p12.rda` to cross-reference variants.
		
		**Panel E**: In the [scRNAseq/inferCNV](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/inferCNV) subdirectory, run `inferCNV.R`, which pulls a counts matrix from the RData object in `PC9.VUDS.10x.counts.matrix.rds` (created in `inferCNV.R`). Necessary annotation and gene order files are also provided.
		
		**Panel F**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories. Scripts to de-multiplex hashed raw data and outputs are provided in the [scRNAseq/HTO_identification](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/HTO_identification) subdirectory. A full matrix of de-multiplexed counts is included as `PC9_scRNAseqCounts_HTOdemux.csv.zip`.
		
		**Panel G**: In the [GO](https://github.com/QuLab-VU/GES_2021/tree/master/GO) folder, run `GO_correlation.R`, which pulls data from `mutations_DEGs-hg38.RData`, a file that compiles all IMPACT genetic mutations (from the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory) and differentially expressed genes (DEGs; from the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory).
		
	- #### <ins>FIGURE 5</ins>
		
		**Panels A and E**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory.
		
		**Panels B and F**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls simulated data from the `trajectories_*.csv` files in the directory. Model trajectories are representative examples of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory).
		
		**Panels C and G**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls simulated data from the `distributions_*.csv` files in the directory. Model distributions were calculated from example trajectories as part of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory). For each subline, the mean and confidence interval reported on the plot is calculated based on 100 bootstrapped p-values provided in one of the `ADbootstrap*.csv` files.
		
		**Panels D and H**: In the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory, run `plotParameterScan.R`, which pulls data from the `*_lowVal.csv` files in the directory.
		
	- #### <ins>FIGURE 6</ins>: *N/A*

- #### <ins>SUPPLEMENTARY FIGURES</ins>

	- #### <ins>SUPPLEMENTARY FIGURE S1</ins>
		
		**Panel A**: Screenshot of the EGFR gene from the [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/) (IGV) based on raw exome sequencing data (avaiable in the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA) at accession \#PRJNA632351). Image is stored as `PC9-EGFRgene_mutations_ex19delCommon.svg` in the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory.
		
		**Panel B**: ***N/A***
		
	- #### <ins>SUPPLEMENTARY FIGURE S2</ins>
		
		**Panels A, B, and C**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory. Data from overlays in panel *C* come from the `PopD_trajectories.RData` object.
		
	- #### <ins>SUPPLEMENTARY FIGURE S3</ins>
		
		**Panel A**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `number_mutations.csv` in the directory.
		
		**Panel B**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `samples_called_vars_named.vcf.gz` in the directory. Directions to download reference FASTA and GTF files are provided in `WES.R`.
		
		**Panel C**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `shared_variants_CLV.csv` in the directory.
		
		**Panel D**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `shared_variants_sublines.csv` in the directory.
		
		**Panel E**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `shared_variants_VUDSlines.csv` in the directory.
		
	- #### <ins>SUPPLEMENTARY FIGURE S4</ins>
		
		**Panels A and B**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from the `vep_*.txt` files in the directory and uses the database in the RData object in `RefCDS_human_GRCH38.p12.rda` to cross-reference variants.
		
	- #### <ins>SUPPLEMENTARY FIGURE S5</ins>
		
		**Panel A**: Screenshot of the summarized output from the [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) quality control analysis on the scRNA-seq library (available in the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) (GEO) data repository at accession \#GSE150084). Settings are shown in the image, which is stored as `CellRanger_PC9.svg` in the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory. 
		
		**Panel B**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		 
	- #### <ins>SUPPLEMENTARY FIGURE S6</ins>
		
		**Panels A and B**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories and subsets data by cell line versions.
		
		**Panels C and D**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories and subsets data by sublines.
		
		**Panels E and F**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		
	- #### <ins>SUPPLEMENTARY FIGURE S7</ins>
		
		**Panels A and B**: In the [RNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/RNAseq) directory, run `RNAseq.R`, which pulls from all 8 `*_featurecounts.txt` files in the directory. These files were created using the Bash script in `RNAseq_processing.txt`. *NOTE: The `*_featurecounts.txt` files must be manually unzipped before running `RNAseq.R`.*
		
	- #### <ins>SUPPLEMENTARY FIGURE S8</ins>
		
		In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories. Input hallmark gene signature (`.gmt`) files can be found in the [scRNAseq/VISION_gmt/hallmark](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/VISION_gmt/hallmark) subdirectory.
		
	- #### <ins>SUPPLEMENTARY FIGURE S9</ins>
		
		In the [GO](https://github.com/QuLab-VU/GES_2021/tree/master/GO) directory, run `semanticSimilarity.R`, which pulls data from `mutations_DEGs-hg38.RData`, a file that compiles all IMPACT genetic mutations (from the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory) and differentially expressed genes (DEGs; from the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory). Directions for downloading reference GTF file are provided in `semanticSimilarity.R`.
		
	- #### <ins>SUPPLEMENTARY FIGURE S10</ins>
		
		In the [scRNAseq/inferCNV](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/inferCNV) directory, run `inferCNV.R`, which pulls a counts matrix from the RData object in `PC9.VUDS.10x.counts.matrix.rds` (created in `inferCNV.R`). Necessary annotation and gene order files are also provided in the directory.
		
	- #### <ins>SUPPLEMENTARY FIGURE S11</ins>: *N/A*
	
	- #### <ins>SUPPLEMENTARY FIGURE S12</ins>
		
		**Panel A**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory.
		
		**Panel B**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory. Model trajectories are representative examples of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory).
		
		**Panel C**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls simulated data from the `distributions_*.csv` files in the directory. Model distributions were calculated from example trajectories as part of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory). For each subline, the mean and confidence interval reported on the plot is calculated based on 100 bootstrapped p-values provided in one of the `ADbootstrap*.csv` files.
		
		**Panel D**: In the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory, run `plotParameterScan.R`, which pulls from the `*_lowVal.csv` files in the directory.
		
	- #### <ins>SUPPLEMENTARY FIGURE S13</ins>
		
		**Panels A and B**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `samples_called_vars_named.vcf.gz` in the directory.
		
	- #### <ins>SUPPLEMENTARY FIGURE S14</ins>
		
		In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		
	- #### <ins>SUPPLEMENTARY FIGURE S15</ins>
		
		In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		
