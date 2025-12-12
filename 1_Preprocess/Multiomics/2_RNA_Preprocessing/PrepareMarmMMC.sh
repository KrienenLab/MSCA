#!/usr/bin/env bash
#SBATCH -J "run-mmc"
#SBATCH -o /jukebox/krienen/daisy/logs/mmc_%j.log
#SBATCH -e /jukebox/krienen/daisy/logs/mmc_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sdan@princeton.edu


module purge;
module load anacondapy/2024.02;
conda activate /scratch/sdan/tools/cell_type_mapper;

## use the interim BG labels as the reference for future incoming data batches
h5ad_file="/scratch/sdan/reference_tmp.h5ad";
output_path='/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/analysis/work/cell_type_mapper_output';
scratch_path='/scratch/sdan/';
species='marmoset';
query_path="/scratch/sdan/query_tmp.h5ad"
save_prefix='HMBA_marmoset_0621'

# Evidently this AllenInstitue 'argschema' module doesn't sanitize command line
# inputs. Instead you have to pass in strings to be evaluated as python 
# expressions.

python -m cell_type_mapper.cli.precompute_stats_scrattch \
--h5ad_path "${h5ad_file}" \
--hierarchy '[ "Subcortex_Class_v4", "Subcortex_Group_v4", "Cluster_v4"]' \
--output_path="${output_path}/${species}_precomputed_stats.h5" \
--normalization raw \
--tmp_dir "${scratch_path}/temp";

python -m cell_type_mapper.cli.reference_markers \
--precomputed_path_list="['${output_path}/${species}_precomputed_stats.h5']" \
--output_dir="${output_path}";

python -m cell_type_mapper.cli.query_markers \
--reference_marker_path_list="['${output_path}/reference_markers.h5']" \
--output_path="${output_path}/${species}_marker.json"

python -m cell_type_mapper.cli.from_specified_markers \
    --query_path="${query_path}" \
    --extended_result_path="${output_path}/${save_prefix}_extended_result_hierachical.json" \
    --csv_result_path="${output_path}/${save_prefix}_result_table_hierachical.csv" \
    --precomputed_stats.path="${output_path}/${species}_precomputed_stats.h5" \
    --query_markers.serialized_lookup="${output_path}/${species}_marker.json" \
    --type_assignment.normalization="raw" 


exit 0;

# ------- #

# python -m cell_type_mapper.cli.precompute_stats_abc --help
# usage: precompute_stats_abc.py [-h] [--hierarchy HIERARCHY] [--output_path OUTPUT_PATH] [--h5ad_path_list H5AD_PATH_LIST] [--cell_metadata_path CELL_METADATA_PATH] [--cluster_annotation_path CLUSTER_ANNOTATION_PATH]
#                                [--cluster_membership_path CLUSTER_MEMBERSHIP_PATH] [--input_json INPUT_JSON] [--output_json OUTPUT_JSON] [--log_level LOG_LEVEL] [--clobber CLOBBER] [--normalization NORMALIZATION] [--n_processors N_PROCESSORS]
#                                [--tmp_dir TMP_DIR] [--split_by_dataset SPLIT_BY_DATASET]

# options:
#   -h, --help            show this help message and exit

# PrecomputedStatsABCSchema:
#   --hierarchy HIERARCHY
#                         List of term_set_labels in our cell types taxonomy ordered from most gross to most fine (default=None) (REQUIRED)
#   --output_path OUTPUT_PATH
#                         Path to the HDF5 file that will be written with the precomputed stats. The serialized taxonomy tree will also be saved here (default=None) (REQUIRED)
#   --h5ad_path_list H5AD_PATH_LIST
#                         List of paths to h5ad files that contain the cell-by-gene data for which we are precomputing statistics (default=None) (REQUIRED)
#   --cell_metadata_path CELL_METADATA_PATH
#                         Path to cell_metadata.csv; the file mapping cells to clusters in our cell types taxonomy. (default=None) (REQUIRED)
#   --cluster_annotation_path CLUSTER_ANNOTATION_PATH
#                         Path to cluster_annotation_term.csv; the file containing parent-child relationships within our cell types taxonomy (default=None) (REQUIRED)
#   --cluster_membership_path CLUSTER_MEMBERSHIP_PATH
#                         Path to cluster_to_cluster_annotation_membership.csv; the file containing the mapping between cluster labels and aliases in our cell types taxonomy (default=None) (REQUIRED)
#   --input_json INPUT_JSON
#                         file path of input json file
#   --output_json OUTPUT_JSON
#                         file path to output json file
#   --log_level LOG_LEVEL
#                         set the logging level of the module (default=ERROR)
#   --clobber CLOBBER     Set to True to allow the code to overwrite an existing file. (default=False)
#   --normalization NORMALIZATION
#                         Normalization of the h5ad files; must be either 'raw' or 'log2CPM' (default=raw)
#   --n_processors N_PROCESSORS
#                         Number of worker processes to spin up. (default=3)
#   --tmp_dir TMP_DIR     Directory where temprorary scratch files will be written out if n_processors > 1 (default=None)
#   --split_by_dataset SPLIT_BY_DATASET
#                         If true, split the dataset by the 'dataset_label' field in cell_metadata.csv, storing each dataset in a separate HDF5 file. Files will be named like output_path but with a secondary suffix added before .h5 specifying
#                         which dataset they contain. (default=False)


# cd /jukebox/krienen/marm_hmba_integration/240905_reprocess_and_recluster/data/ortho/;

# out_dir="/jukebox/krienen/marm_hmba_integration/240905_reprocess_and_recluster/data/mapmycell/macaque_BG";
# [ -e $out_dir ] || mkdir $out_dir

# SAMPLE="rna_merged_qc_pass";

# python -m cell_type_mapper.cli.from_specified_markers \
#     --query_path="${SAMPLE}_macaque_ortho.h5ad" \
#     --extended_result_path="${out_dir}/${SAMPLE}_extended_result_hierachical.json" \
#     --csv_result_path="${out_dir}/${SAMPLE}_result_table_hierachical.csv" \
#     --precomputed_stats.path="/jukebox/krienen/abc_atlas/taxonomies_082024/mapmycell/macaque_precomputed_stats.h5" \
#     --query_markers.serialized_lookup="/jukebox/krienen/abc_atlas/taxonomies_082024/mapmycell/macaque_marker.json" \
#     --type_assignment.normalization="raw" 
