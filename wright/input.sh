#!/bin/bash

#SBATCH -p general
#SBATCH -J test_rmg
##SBATCH -n 24
#SBATCH --mem=20GB
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -o wright2009.slurm.log
#SBATCH -e wright2009.slurm.log


export PYTHONPATH=$PYTHONPATH:$HOME/Code/RMG-Py/

source activate rmg_env
python $HOME/Code/RMG-Py/rmg.py wright2009.py > wright2009.log
source deactivate
