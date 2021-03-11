###Data repository for &quot;An *in vitro* model of tumor heterogeneity resolves genetic, epigenetic, and stochastic sources of cell state variability,&quot; Hayford et al. (2021), *[journal vol : article]; [DOI: xxx](http://dx.doi.org/xxx)*

---

####\**Instructions for creating panels in all main and supplementary figures based on experimental and simulated data*

- ##### <u>MAIN FIGURES </u>

	- ##### FIGURE 1: *N/A*
	
	- ##### FIGURE 2
		
		**Panels A and C**: In the [DrugResponse](https://github.com/QuLab-VU/GES_2021/tree/master/DrugResponse) directory, run ``DrugResponse.R``, which pulls data from the two `Parental-*.csv` files in the directory and the well conditions in the [DrugResponse/Platemaps](https://github.com/QuLab-VU/GES_2021/tree/master/DrugResponse/Platemaps) subdirectory.
		
		**Panels B and D**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run ``cFP.R``, which pulls data from the 10 `cFP_*.csv` files in the directory.
		
	- ##### FIGURE 3
		
		**Panel A**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `mutations_byChromosome.csv`.
		
		**Panels B, C, and D**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from the `vep_*.txt` files in the directory and uses the database in the RData object in `RefCDS_human_GRCH38.p12.rda` to cross-reference variants.
		
		**Panel E**: In the [scRNAseq/inferCNV](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/inferCNV) subdirectory, run `inferCNV.R`, which pulls a counts matrix from the RData object in `PC9.CLV.10x.counts.matrix.rds`, included in the directory. Necessary annotation and gene order files are also provided.
		
		**Panel F**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories. Scripts to de-multiplex hashed raw data and outputs are included in the [scRNAseq/HTO_identification](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/HTO_identification) subdirectory. A full matrix of de-multiplexed counts is included as `PC9_scRNAseqCounts_HTOdemux.csv.zip`.
		
		**Panel G**: In the [GO](https://github.com/QuLab-VU/GES_2021/tree/master/GO) directory, run `GO_correlation.R`, which pulls data from `mutations_DEGs-hg38.RData`, a file that compiles all IMPACT genetic mutations (from the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory) and differentially expressed genes (DEGs; from the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory).
		
	- ##### FIGURE 4
		
		**Panel A**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) folder, run `WES.R`, which pulls data from `mutations_byChromosome.csv`.
		
		**Panels B, C, and D**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from the `vep_*.txt` files in the directory and uses the database in the RData object in `RefCDS_human_GRCH38.p12.rda` to cross-reference variants.
		
		**Panel E**: In the [scRNAseq/inferCNV](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/inferCNV) subdirectory, run `inferCNV.R`, which pulls a counts matrix from the RData object in `PC9.VUDS.10x.counts.matrix.rds` (created in `inferCNV.R`). Necessary annotation and gene order files are also provided.
		
		**Panel F**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories. Scripts to de-multiplex hashed raw data and outputs are provided in the [scRNAseq/HTO_identification](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/HTO_identification) subdirectory. A full matrix of de-multiplexed counts is included as `PC9_scRNAseqCounts_HTOdemux.csv.zip`.
		
		**Panel G**: In the [GO](https://github.com/QuLab-VU/GES_2021/tree/master/GO) folder, run `GO_correlation.R`, which pulls data from `mutations_DEGs-hg38.RData`, a file that compiles all IMPACT genetic mutations (from the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory) and differentially expressed genes (DEGs; from the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory).
		
	- ##### FIGURE 5
		
		**Panels A and E**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory.
		
		**Panels B and F**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls simulated data from the `trajectories_*.csv` files in the directory. Model trajectories are representative examples of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory).
		
		**Panels C and G**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls simulated data from the `distributions_*.csv` files in the directory. Model distributions were calculated from example trajectories as part of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory). Statistics come from the `ADbootstrap*.csv` files.
		
		**Panels D and H**: In the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory, run `plotParameterScan.R`, which pulls data from the `*_lowVal.csv` files in the directory.
		
	- ##### FIGURE 6: *N/A*

- ##### <u>SUPPLEMENTARY FIGURES</u>

	- ##### SUPPLEMENTARY FIGURE S1
		
		**Panel A**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, screenshot from **\*\[TOOL\]\***, stored as `PC9-EGFRgene_mutations_ex19delCommon.svg`.
		
		**Panel B**: ***N/A***
		
	- ##### SUPPLEMENTARY FIGURE S2
		
		**Panels A, B, and C**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory. Data from overlays in panel *C* come from the `PopD_trajectories.RData` object.
		
	- ##### SUPPLEMENTARY FIGURE S3
		
		**Panel A**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `number_mutations.csv` in the directory.
		.
		**Panel B**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `samples_called_vars_named.vcf.gz` in the directory. Directions to download reference FASTA and GTF files are provided in `WES.R`.
		
		**Panel C**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `shared_variants_CLV.csv` in the directory.
		
		**Panel D**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `shared_variants_sublines.csv` in the directory.
		
		**Panel E**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `shared_variants_VUDSlines.csv` in the directory.
		
	- ##### SUPPLEMENTARY FIGURE S4
		
		**Panels A and B**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from the `vep_*.txt` files in the directory and uses the database in the RData object in `RefCDS_human_GRCH38.p12.rda` to cross-reference variants.
		
	- ##### SUPPLEMENTARY FIGURE S5
		
		**Panel A**: Screenshot from the CellRanger output HTML
		
		**Panel B**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		
	- ##### SUPPLEMENTARY FIGURE S6
		
		**Panels A and B**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories and subsets data by cell line versions.
		
		**Panels C and D**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories and subsets data by sublines.
		
		**Panels E and F**: In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		
	- ##### SUPPLEMENTARY FIGURE S7
		
		**Panels A and B**: In the [RNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/RNAseq) directory, run `RNAseq.R`, which pulls from all 8 `*_featurecounts.txt` files in the directory. These files were created using the Bash script in `RNAseq_processing.txt`.
		
	- ##### SUPPLEMENTARY FIGURE S8
		
		In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories. Input hallmark gene signature (`.gmt`) files can be found in the [scRNAseq/VISION_gmt/hallmark](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/VISION_gmt/hallmark) subdirectory.
		
	- ##### SUPPLEMENTARY FIGURE S9
		
		In the [GO](https://github.com/QuLab-VU/GES_2021/tree/master/GO) directory, run `semanticSimilarity.R`, which pulls data from `mutations_DEGs-hg38.RData`, a file that compiles all IMPACT genetic mutations (from the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory) and differentially expressed genes (DEGs; from the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory). Directions for downloading reference GTF file are provided in `semanticSimilarity.R`.
		
	- ##### SUPPLEMENTARY FIGURE S10
		
		In the [scRNAseq/inferCNV](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/inferCNV) directory, run `inferCNV.R`, which pulls a counts matrix from the RData object in `PC9.VUDS.10x.counts.matrix.rds` (created in `inferCNV.R`). Necessary annotation and gene order files are also provided in the directory.
		
	- ##### SUPPLEMENTARY FIGURE S11: *N/A*
	
	- ##### SUPPLEMENTARY FIGURE S12
		
		**Panel A**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory.
		
		**Panel B**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls data from the `trajectories_*.csv` files in the directory. Model trajectories are representative examples of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory).
		
		**Panel C**: In the [cFP](https://github.com/QuLab-VU/GES_2021/tree/master/cFP) directory, run `cFP.R`, which pulls simulated data from the `distributions_*.csv` files in the directory. Model distributions were calculated from example trajectories as part of a larger simulation scan (`*.py` models in the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory). Statistics come from `ADbootstrap*.csv` files.
		
		**Panel D**: In the [Simulations](https://github.com/QuLab-VU/GES_2021/tree/master/Simulations) directory, run `plotParameterScan.R`, which pulls from the `*_lowVal.csv` files in the directory.
		
	- ##### SUPPLEMENTARY FIGURE S13
		
		**Panels A and B**: In the [WES](https://github.com/QuLab-VU/GES_2021/tree/master/WES) directory, run `WES.R`, which pulls data from `samples_called_vars_named.vcf.gz` in the directory.
		
	- ##### SUPPLEMENTARY FIGURE S14
		
		In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		
	- ##### SUPPLEMENTARY FIGURE S15
		
		In the [scRNAseq](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq) directory, run `scRNAseq.R`, which pulls from 10x Genomics reduced data in the [scRNAseq/read_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/read_count) and [scRNAseq/umi_count](https://github.com/QuLab-VU/GES_2021/tree/master/scRNAseq/umi_count) subdirectories.
		