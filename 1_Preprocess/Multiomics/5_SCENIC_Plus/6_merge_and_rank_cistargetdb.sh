#!/bin/bash
#SBATCH -o ./job_outputs/slurm_5_merge_and_rank_cistargetdb_%j.out
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mem-per-cpu=250G       # 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=vn0027@princeton.edu

module load anacondapy/2024.02
conda activate /jukebox/scratch/vn0027/envs/scenicplus

cd /jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db

export PATH="/jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db:$PATH"

base_path="/jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db"
script_path="${base_path}/create_cisTarget_databases"
output_name="marm_subcortex"
db_dir="${base_path}/dbs"

mkdir -p "${db_dir}"

# Step 1: Merge partial motifs_vs_regions scores
python "${script_path}/combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py" \
    -i "${output_name}" \
    -o "${db_dir}/"

# Step 2: Clean up intermediate part files (optional but recommended)
#rm -f "${output_name}"*part*

# Step 3: Convert scores to rankings
python "${script_path}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py" \
    -i "${db_dir}/${output_name}.motifs_vs_regions.scores.feather" \
    -s 555
    