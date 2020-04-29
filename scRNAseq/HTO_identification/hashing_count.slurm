#!/bin/bash
#SBATCH --mail-user=corey.hayford@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hashing_count.txt

cd /data/lola/hayforc/scRNAseq/PC9
module load Anaconda3/5.0.1
source activate RNAvelocityEnvReal
CITE-seq-Count -R1 2958-CH-1-hash-ATTACTCG-TCTTTCCC_S381_R1_001.fastq.gz -R2 2958-CH-1-hash-ATTACTCG-TCTTTCCC_S381_R2_001.fastq.gz -t tags.csv -cbf 1 -cbl 16 -umif 17 -umil 26 -cells 15400 -o /data/lola/hayforc/scRNAseq/PC9/
