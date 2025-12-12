#!/usr/bin/env bash
#SBATCH -o slurm_snap_callpeaks_%j.out
#SBATCH --time=72:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=1        # 
#SBATCH --mem=600G               # 

cd /jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/slurm_logs_atac
source $(conda info --base)/etc/profile.d/conda.sh
conda activate scvi
python ../pipeline_atac/call_peaks.py
