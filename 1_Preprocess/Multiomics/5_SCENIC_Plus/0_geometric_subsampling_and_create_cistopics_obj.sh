#!/bin/bash
#SBATCH -o ./job_outputs/slurm_0_geometric_subsampling_and_create_cistopics_obj_%j.out
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mem=500G       # 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=vn0027@princeton.edu

module load anacondapy/2024.02
conda activate /jukebox/scratch/vn0027/envs/scenicplus

cd /jukebox/krienen/victor/hmba/pycisTopic

python3 /jukebox/krienen/victor/marm_hmba_atac_50k/scripts/0_geometric_subsampling_and_create_cistopics_obj.py

