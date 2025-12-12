#!/usr/bin/env bash
#SBATCH -o slurm-rscript_%j.out
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=1        # 
#SBATCH --mem-per-cpu=240G       # 

cd /jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster
Rscript $1
