#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem=16G
#SBATCH --output=PC9-DS8_parameterscan_twoStateExpansion_Divs11-allDIP.out

module load Anaconda2/5.0.1
source activate RNAvelocityEnvReal

cd /data/lola/hayforc/ParamScan/DS8_mod
python PC9-DS8_parameterscan_twoStateExpansion_Divs11-allDIP.py
