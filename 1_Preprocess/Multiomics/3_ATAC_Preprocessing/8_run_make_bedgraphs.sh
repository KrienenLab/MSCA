#!/usr/bin/env bash
#SBATCH -o snap_bw_%j.out
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=1        # 
#SBATCH --mem=950G               # 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=md6347@princeton.edu

PDIR=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster
source $(conda info --base)/etc/profile.d/conda.sh
conda activate scvi
export RUST_BACKTRACE="1"
export RUST_LOG="bigtools=debug,snapatac2=info"
python ${PDIR}/pipeline_atac/make_bedgraphs.py
