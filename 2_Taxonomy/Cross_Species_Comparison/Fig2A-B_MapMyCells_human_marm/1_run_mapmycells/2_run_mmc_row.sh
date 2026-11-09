#!/usr/bin/env bash
#SBATCH -J mmc
#SBATCH -o logs/mmc_%A_%a.log
#SBATCH -e logs/mmc_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=4

module purge
module load anacondapy/2024.02
module load cudatoolkit/12.4.1
conda activate /scratch/sdan/tools/cell_type_mapper

set -euo pipefail

proj_dir="${1:?Usage: $0 PROJ_DIR [CSV_ROW_NUMBER]}"
reference_data_table="${proj_dir}/pipeline_rna/mmc/mapmycells_reference.csv"

# CSV row to process: the Slurm array task id, or the second argument for manual runs
row_number="${2:-${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}}"

data="${proj_dir}/data/rna/rna_raw_combined.h5ad"
ortho_output_dir="${proj_dir}/data/rna/ortho"
out_dir="${proj_dir}/data/rna/mapmycells"

mkdir -p "${out_dir}" logs

# comma-separated row (fields must not contain commas); strip CR and BOM
line="$(sed -n "${row_number}p" "${reference_data_table}" | tr -d '\r' | sed 's/^\xef\xbb\xbf//')"

if [[ -z "${line}" ]]; then
    echo "Error: CSV row ${row_number} is empty or does not exist" >&2
    exit 1
fi

IFS=',' read -r reference_name species precomputed_stats_path query_markers_path gene_type <<< "${line}"

# trim surrounding whitespace
reference_name="$(xargs <<< "${reference_name:-}")"
species="$(xargs <<< "${species:-}")"
precomputed_stats_path="$(xargs <<< "${precomputed_stats_path:-}")"
query_markers_path="$(xargs <<< "${query_markers_path:-}")"
gene_type="$(xargs <<< "${gene_type:-}")"

if [[ -z "${reference_name}" || -z "${species}" || \
      -z "${precomputed_stats_path}" ]]; then
    echo "Error: Missing required field(s) in CSV row ${row_number}: ${line}" >&2
    exit 1
fi

result_csv="${out_dir}/${reference_name}_result_table_hierachical.csv"
extended_json="${out_dir}/${reference_name}_extended_result_hierachical.json"

if [[ -s "${result_csv}" ]]; then
    echo "Output already exists; skipping: ${result_csv}"
    exit 0
fi

if [[ "${species}" == "marmoset" ]]; then
    query_path="${data}"
else
    if [[ -z "${gene_type}" ]]; then
        echo "Error: gene_type is required for non-marmoset row ${row_number}" >&2
        exit 1
    fi
    query_path="${ortho_output_dir}/rna_${species}_ortho_${gene_type}.h5ad"
fi

for required_file in "${query_path}" "${precomputed_stats_path}"; do
    if [[ ! -s "${required_file}" ]]; then
        echo "Error: Required file is missing or empty: ${required_file}" >&2
        exit 1
    fi
done

# drop the supertype level during mapping if the reference taxonomy has one: a level whose label or name
# is 'supertype', or whose label is 'SUPT' / ends in '_SUPT' (case-insensitive); the leaf level is never dropped
drop_level="$(python - "${precomputed_stats_path}" <<'PYLEVEL'
import json
import sys

import h5py

with h5py.File(sys.argv[1], "r") as src:
    tree = json.loads(src["taxonomy_tree"][()].decode("utf-8"))
level_names = tree.get("hierarchy_mapper", {})
for level in tree["hierarchy"][:-1]:
    label = level.lower()
    name = str(level_names.get(level, "")).lower()
    if "supertype" in (label, name) or label == "supt" or label.endswith("_supt"):
        print(level)
        break
PYLEVEL
)"
drop_args=()
[[ -n "${drop_level}" ]] && drop_args=(--drop_level="${drop_level}")

# use the query marker JSON if given; otherwise find markers on the fly
if [[ -n "${query_markers_path}" ]]; then
    if [[ ! -s "${query_markers_path}" ]]; then
        echo "Error: Marker JSON was specified but is missing or empty: ${query_markers_path}" >&2
        exit 1
    fi
    mapping_module="cell_type_mapper.cli.from_specified_markers"
    mapping_mode="specified markers"
    marker_args=("--query_markers.serialized_lookup=${query_markers_path}")
else
    mapping_module="cell_type_mapper.cli.map_to_on_the_fly_markers"
    mapping_mode="on-the-fly markers"
    marker_args=()
fi

echo "[$(date)] CSV row ${row_number}: ${reference_name}"
echo "Species: ${species}"
echo "Query: ${query_path}"
echo "Reference stats: ${precomputed_stats_path}"
echo "Mapping mode: ${mapping_mode}"
echo "Dropped level: ${drop_level:-none}"

if [[ "${mapping_mode}" == "on-the-fly markers" ]]; then
    # the nested marker-finder settings can only be set in a JSON config (empty dicts use defaults)
    config_json="$(mktemp "${out_dir}/.${reference_name}.otf.XXXXXX.json")"
    trap 'rm -f "${config_json}"' EXIT
    python - "${config_json}" "${query_path}" "${extended_json}" \
        "${result_csv}" "${precomputed_stats_path}" "${drop_level}" <<'PYCONFIG'
import json
import sys
config_path, query_path, extended_json, result_csv, stats_path, drop_level = sys.argv[1:]
config = {
    "query_path": query_path,
    "extended_result_path": extended_json,
    "csv_result_path": result_csv,
    "precomputed_stats": {"path": stats_path},
    "reference_markers": {},
    "query_markers": {},
    "type_assignment": {"normalization": "raw"},
    "n_processors": int(__import__("os").environ.get("SLURM_CPUS_PER_TASK", "4")),
}
if drop_level:
    config["drop_level"] = drop_level
with open(config_path, "w") as dst:
    json.dump(config, dst, indent=2)
PYCONFIG
    python -m "${mapping_module}" --input_json "${config_json}"
else
    python -m "${mapping_module}" \
        --query_path="${query_path}" \
        --extended_result_path="${extended_json}" \
        --csv_result_path="${result_csv}" \
        --precomputed_stats.path="${precomputed_stats_path}" \
        "${marker_args[@]}" \
        ${drop_args[@]+"${drop_args[@]}"} \
        --type_assignment.normalization="raw"
fi

echo "[$(date)] Finished ${reference_name} (${mapping_mode})"
