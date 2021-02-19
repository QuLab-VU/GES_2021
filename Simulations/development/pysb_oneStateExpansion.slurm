#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --output=pysb_oneStateExpansion_ACCRE.out

module load Anaconda2/5.0.1
source activate RNAvelocityEnvReal

cd /data/lola/hayforc/ParamScan/DSlines_expanded
python pysb_oneStateExpansion.py
