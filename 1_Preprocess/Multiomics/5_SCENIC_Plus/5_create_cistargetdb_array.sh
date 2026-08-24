#!/bin/bash
#SBATCH --job-name=cistarget_array
#SBATCH --output=./job_outputs/slurm_4_cistarget_array_%A_%a.out
#SBATCH --array=1-10
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vn0027@princeton.edu

'''
## To prepare fasta from consensus regions:

cd /jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db
REGION_BED="/jukebox/krienen/victor/marm_hmba_atac_50k/data/marm_subcortex_filtered_peaks.bed"
GENOME_FASTA="/jukebox/krienen/genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/mCalJac1-pat-X_rs_mitos2/fasta/genome.fa"
CHROMSIZES="/jukebox/krienen/genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/mCalJac1-pat-X_rs_mitos2/star/chrNameLength.txt"
DATABASE_PREFIX="marm_subcortex_1kb_bg_with_mask"
SCRIPT_DIR="/jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        marm_subcortex.with_1kb_bg_padding.filtered.fa \
        1000 \
        yes

'''      

module load anacondapy/2024.02
conda activate /jukebox/scratch/vn0027/envs/scenicplus

cd /jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db

export PATH="/jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db:$PATH"

base_path=/jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db
script_path=${base_path}/create_cisTarget_databases
output_name=marm_subcortex
fasta_file="${base_path}/marm_subcortex.with_1kb_bg_padding.filtered.fa"
motif_db="${base_path}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
motif_list="${base_path}/motifs.txt"

# Use SLURM_ARRAY_TASK_ID as the current part
python ${script_path}/create_cistarget_motif_databases.py \
    -f ${fasta_file} \
    -M ${motif_db} \
    -m ${motif_list} \
    -p ${SLURM_ARRAY_TASK_ID} 10 \
    -o ${output_name} \
    -t 20
    