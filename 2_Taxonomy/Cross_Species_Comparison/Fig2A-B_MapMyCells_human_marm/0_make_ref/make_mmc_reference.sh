#!/usr/bin/env bash
#SBATCH -J mmc_make_ref
#SBATCH -o mmc_make_ref_%j.out
#SBATCH -p all
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G

###############################################################################
# Build a MapMyCells (cell_type_mapper) reference.
#
# All MapMyCells reference files used in this project (precomputed stats,
# reference markers, query marker lookup) were made with this pipeline.
#
# Steps
#   1. precompute stats  : per-cluster expression statistics of the reference
#   2. reference_markers : marker genes for every pair of sibling taxonomy nodes
#   3. query_markers     : marker lookup (json) used when mapping query data
#
# Reference input, selected with INPUT_FORMAT
#   scrattch : one raw-count h5ad (cells x genes) whose obs columns hold the
#              taxonomy levels                               (cli: precompute_stats_scrattch)
#                H5AD_FILE   path to the h5ad
#                HIERARCHY   obs columns, coarsest to finest, comma separated
#   abc      : Allen Brain Cell Atlas file format            (cli: precompute_stats_abc)
#                H5AD_DIRS   directories searched for raw-count h5ads (space separated)
#                H5AD_PATTERN  file name pattern (default '*raw.h5ad')
#                CELL_METADATA, CLUSTER_ANNOTATION, CLUSTER_MEMBERSHIP
#                            cell_metadata.csv, cluster_annotation_term.csv,
#                            cluster_to_cluster_annotation_membership.csv
#                HIERARCHY   term_set labels, coarsest to finest, comma separated
#
# Always required: REF_NAME, OUT_DIR, INPUT_FORMAT, CONDA_ENV
# Optional (only passed to cell_type_mapper when set; its defaults apply otherwise):
#   NORMALIZATION ('raw' or 'log2CPM'; precompute), N_PROCESSORS, TMP_DIR,
#   CLOBBER=True (overwrite existing outputs; otherwise finished steps are skipped)
#
# Example
#   REF_NAME=my_ref OUT_DIR=/path/to/out CONDA_ENV=my_env INPUT_FORMAT=scrattch \
#   H5AD_FILE=/path/to/ref.h5ad HIERARCHY="class,subclass,cluster" \
#   sbatch make_mmc_reference.sh
#
# Outputs (in OUT_DIR)
#   ${REF_NAME}_precomputed_stats.h5, reference_markers.h5, ${REF_NAME}_marker_genes.json
#   (reference_markers.h5 has a fixed name, so use one OUT_DIR per reference)
###############################################################################

module purge
module load "${ANACONDA_MODULE:-anacondapy/2024.02}"
conda activate "${CONDA_ENV:?Set CONDA_ENV (environment with cell_type_mapper)}"

set -euo pipefail

REF_NAME="${REF_NAME:?Set REF_NAME}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR}"
INPUT_FORMAT="${INPUT_FORMAT:?Set INPUT_FORMAT to 'scrattch' or 'abc'}"
HIERARCHY="${HIERARCHY:?Set HIERARCHY}"

mkdir -p "${OUT_DIR}"
precomputed_stats="${OUT_DIR}/${REF_NAME}_precomputed_stats.h5"
reference_markers="${OUT_DIR}/reference_markers.h5"
query_markers="${OUT_DIR}/${REF_NAME}_marker_genes.json"

# cell_type_mapper takes list arguments as python expressions
to_py_list() { printf "['%s']" "$(paste -sd, - | sed "s/,/', '/g")"; }
hierarchy_list="$(tr ',' '\n' <<< "${HIERARCHY}" | to_py_list)"

# optional arguments, passed only when set
clobber_opt=();   [[ -n "${CLOBBER:-}"       ]] && clobber_opt=(--clobber="${CLOBBER}")
tmp_opt=();       [[ -n "${TMP_DIR:-}"       ]] && tmp_opt=(--tmp_dir="${TMP_DIR}")
nproc_opt=();     [[ -n "${N_PROCESSORS:-}"  ]] && nproc_opt=(--n_processors="${N_PROCESSORS}")
norm_opt=();      [[ -n "${NORMALIZATION:-}" ]] && norm_opt=(--normalization="${NORMALIZATION}")

# finished steps are skipped unless CLOBBER=True
run_step() { [[ ! -s "$1" || "${CLOBBER:-}" == "True" ]]; }

# 1. precompute cluster statistics
if run_step "${precomputed_stats}"; then
    case "${INPUT_FORMAT}" in
        scrattch)
            python -m cell_type_mapper.cli.precompute_stats_scrattch \
                --h5ad_path="${H5AD_FILE:?Set H5AD_FILE}" \
                --hierarchy="${hierarchy_list}" \
                --output_path="${precomputed_stats}" \
                ${norm_opt[@]+"${norm_opt[@]}"} ${clobber_opt[@]+"${clobber_opt[@]}"} ${tmp_opt[@]+"${tmp_opt[@]}"}
            ;;
        abc)
            # shellcheck disable=SC2086
            h5ad_list="$(find ${H5AD_DIRS:?Set H5AD_DIRS} -name "${H5AD_PATTERN:-*raw.h5ad}" | sort | to_py_list)"
            [[ "${h5ad_list}" != "['']" ]] || { echo "No h5ad files found in ${H5AD_DIRS}" >&2; exit 1; }
            python -m cell_type_mapper.cli.precompute_stats_abc \
                --hierarchy="${hierarchy_list}" \
                --output_path="${precomputed_stats}" \
                --h5ad_path_list="${h5ad_list}" \
                --cell_metadata_path="${CELL_METADATA:?Set CELL_METADATA}" \
                --cluster_annotation_path="${CLUSTER_ANNOTATION:?Set CLUSTER_ANNOTATION}" \
                --cluster_membership_path="${CLUSTER_MEMBERSHIP:?Set CLUSTER_MEMBERSHIP}" \
                ${norm_opt[@]+"${norm_opt[@]}"} ${clobber_opt[@]+"${clobber_opt[@]}"} ${tmp_opt[@]+"${tmp_opt[@]}"} ${nproc_opt[@]+"${nproc_opt[@]}"}
            ;;
        *)
            echo "INPUT_FORMAT must be 'scrattch' or 'abc'" >&2; exit 1
            ;;
    esac
fi

# 2. reference markers
if run_step "${reference_markers}"; then
    python -m cell_type_mapper.cli.reference_markers \
        --precomputed_path_list="['${precomputed_stats}']" \
        --output_dir="${OUT_DIR}" \
        ${clobber_opt[@]+"${clobber_opt[@]}"} ${tmp_opt[@]+"${tmp_opt[@]}"} ${nproc_opt[@]+"${nproc_opt[@]}"}
fi

# 3. query markers
if run_step "${query_markers}"; then
    python -m cell_type_mapper.cli.query_markers \
        --reference_marker_path_list="['${reference_markers}']" \
        --output_path="${query_markers}" \
        ${tmp_opt[@]+"${tmp_opt[@]}"} ${nproc_opt[@]+"${nproc_opt[@]}"}
fi

echo "Reference files written to ${OUT_DIR}"
