#!/bin/bash
#SBATCH -o ./job_outputs/slurm_7_scplus_grn_inference_%j.out
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=12
#SBATCH --mem=500G
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=vn0027@princeton.edu

module load anacondapy/2024.02
conda activate /jukebox/scratch/vn0027/envs/scenicplus

cd /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/Snakemake

#scenicplus grn_inference AUCell \
 #   --eRegulon_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/eRegulon_direct.tsv \
  #  --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/ACC_GEX.h5mu \
   # --aucell_out_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/AUCell_direct.h5mu \
    #--n_cpu 12

scenicplus grn_inference AUCell \
    --eRegulon_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/eRegulons_extended.tsv \
    --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/ACC_GEX.h5mu \
    --aucell_out_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/AUCell_extended.h5mu \
    --n_cpu 12

scenicplus grn_inference create_scplus_mudata \
    --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/ACC_GEX.h5mu \
    --e_regulon_auc_direct_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/AUCell_direct.h5mu \
    --e_regulon_auc_extended_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/AUCell_extended.h5mu \
    --e_regulon_metadata_direct_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/eRegulon_direct.tsv \
    --e_regulon_metadata_extended_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/eRegulons_extended.tsv \
    --out_file /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_neuron/outs/scplus_mdata.h5mu