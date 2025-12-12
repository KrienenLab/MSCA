#!/usr/bin/env bash
#SBATCH -o cellr-a_%j.out
#SBATCH -p all
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=8        # 
#SBATCH --mem-per-cpu=4G         # 

# Make sure to run this script within the directory made for each sample, for instance: 
#   .../alignments/240109_cellranger-arc/THE_SAMPLE_NAME;
# 
# and run this script with the sample name as a positional argument, such as:
#   sbatch this_script.sh THE_SAMPLE_NAME

module purge
module load cellranger-arc/2.0.2

cellranger-arc count \
--id=$1 \
--reference=/jukebox/krienen/genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/mCalJac1-pat-X_rs_mitos2 \
--libraries=./libraries.csv \
--localcores=8 \
--localmem=30;

exit 0
