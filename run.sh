#!/bin/bash
#SBATCH -p gen
#SBATCH --constraint=broadwell
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --output=/mnt/home/atanelus/project/logs/%j.log
#SBATCH --mem=128000mb
pwd; hostname; date;

module load matlab;
cd ~/project
matlab -nodisplay -nosplash -nodesktop -r "run('SEM_run.m');exit;"

