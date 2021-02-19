#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem=16G
#SBATCH --output=PC9-DS1_param-scan_Expansion.out

module load Anaconda2/5.0.1
source activate RNAvelocityEnvReal

cd /data/lola/hayforc/ParamScan/DSlines_expanded
python PC9-DS1_parameterscan_Expansion.py
