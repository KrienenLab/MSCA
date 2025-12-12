#!/usr/bin/env bash
#SBATCH -o slurm-python_gpu_%j.out
#SBATCH --ntasks=1
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=12       # 
#SBATCH --mem=500G               # 
#SBATCH --gres=gpu:1             #

cd /jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster
conda activate scvi
python $1
