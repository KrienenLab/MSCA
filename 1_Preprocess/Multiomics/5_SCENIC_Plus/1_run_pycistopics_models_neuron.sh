#!/bin/bash
#SBATCH -o ./job_outputs/slurm_1_pycistopics_models_neuron_%j.out
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=10        
#SBATCH --mem=250G        
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=vn0027@princeton.edu

module load anacondapy/2024.02
conda activate /jukebox/scratch/vn0027/envs/scenicplus

cd /jukebox/krienen/victor/marm_hmba_atac_50k/pycisTopic/neuron
mkdir ./outs

# Set custom TMPDIR to avoid conflict in /tmp and reuse errors
export TMPDIR="/jukebox/scratch/vn0027/tmp/mallet_models_neuron"
mkdir -p "$TMPDIR"

pycistopic topic_modeling mallet \
    --input "/jukebox/krienen/victor/marm_hmba_atac_50k/data/cistopic_obj_neuron_geosketch_50kcells_marm_subcortex.pkl" \
    --output "/jukebox/krienen/victor/marm_hmba_atac_50k/pycisTopic/neuron/outs/models.pkl" \
    --topics 200 175 150 125 100 50  \
    --parallel 10 \
    --keep True \
    --mallet_path "/jukebox/scratch/vn0027/Mallet-202108/bin/mallet"
