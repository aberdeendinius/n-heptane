#!/bin/bash

#SBATCH -p interactive
#SBATCH -J test_rmg
##SBATCH -n 24
#SBATCH --mem=20GB
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -o zhang2013_2atm.slurm.log
#SBATCH -e zhang2013_2atm.slurm.log


export PYTHONPATH=$PYTHONPATH:$HOME/Code/RMG-Py/

source activate rmg_env
python $HOME/Code/RMG-Py/rmg.py zhang2013_2atm.py > zhang2013_2atm.log
source deactivate
