#!/usr/bin/env bash
#SBATCH -o slurm_snap_fragments_%j.out
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=md6347@princeton.edu
PDIR=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster
sample_name=$1
echo "Running sample: $sample_name"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /jukebox/scratch/md6347/envs/scvi_daisy
python3 ${PDIR}/pipeline_atac/make_one_snap_fragment.py "$sample_name"

