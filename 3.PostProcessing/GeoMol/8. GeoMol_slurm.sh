#!/bin/bash
#SBATCH -J RemoveJob
#SBATCH -p HUNTER,SWAN,CLUSTER
#SBATCH -N 3
#SBATCH -t 48:00:00

module unload anaconda/python3/2022.05v
module add schrodinger/2022-2
source ~/miniconda3/etc/profile.d/conda.sh
conda activate GeoMol

python GeoMol_MacroModel.py
