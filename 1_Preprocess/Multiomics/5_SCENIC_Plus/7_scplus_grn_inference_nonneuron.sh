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

cd /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/Snakemake

#scenicplus grn_inference region_to_gene \
 #   --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/ACC_GEX.h5mu \
  #  --search_space_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/search_space.tsv \
   # --out_region_to_gene_adjacencies /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/region_to_gene_adj.tsv \
    #--temp_dir "/jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/tmp" \
    #--n_cpu 12

#scenicplus grn_inference TF_to_gene \
 #   --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/ACC_GEX.h5mu \
  #  --tf_names /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/tf_names.txt \
   # --out_tf_to_gene_adjacencies /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/tf_to_gene_adj.tsv \
    #--n_cpu 12 \
    #--temp_dir "/jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/tmp" \
    #--seed 666

#grn_inference eGRN 
#grn_inference eGRN extended

#scenicplus grn_inference AUCell \
 #   --eRegulon_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/eRegulon_direct.tsv \
  #  --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/ACC_GEX.h5mu \
   # --aucell_out_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/AUCell_direct.h5mu \
    #--n_cpu 12

scenicplus grn_inference eGRN \
    --is_extended \
    --TF_to_gene_adj_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/tf_to_gene_adj.tsv \
    --region_to_gene_adj_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/region_to_gene_adj.tsv \
    --cistromes_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/cistromes_extended.h5ad \
    --ranking_db_fname /jukebox/krienen/victor/marm_hmba_atac_50k/cistarget_db/dbs/marm_subcortex.regions_vs_motifs.rankings.feather \
    --eRegulon_out_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/eRegulons_extended.tsv \
    --temp_dir "/jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/tmp" \
    --n_cpu 12

scenicplus grn_inference AUCell \
    --eRegulon_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/eRegulons_extended.tsv \
    --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/ACC_GEX.h5mu \
    --aucell_out_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/AUCell_extended.h5mu \
    --n_cpu 12

scenicplus grn_inference create_scplus_mudata \
    --multiome_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/ACC_GEX.h5mu \
    --e_regulon_auc_direct_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/AUCell_direct.h5mu \
    --e_regulon_auc_extended_mudata_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/AUCell_extended.h5mu \
    --e_regulon_metadata_direct_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/eRegulon_direct.tsv \
    --e_regulon_metadata_extended_fname /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/eRegulons_extended.tsv \
    --out_file /jukebox/krienen/victor/marm_hmba_atac_50k/scplus_nonneuron/outs/scplus_mdata.h5mu
