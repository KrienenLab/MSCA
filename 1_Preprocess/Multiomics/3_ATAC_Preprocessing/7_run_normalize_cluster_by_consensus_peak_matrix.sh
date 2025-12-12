#!/usr/bin/env bash
#SBATCH -o slurm_snap_norm_%j.out
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=1        # 
#SBATCH --mem-per-cpu=500G       # 

cd /jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/slurm_logs_atac
source $(conda info --base)/etc/profile.d/conda.sh
conda activate crested
python ../pipeline_atac/normalize_cluster_by_consensus_peak_matrix.py
