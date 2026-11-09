#!/usr/bin/env bash

###############################################################################
# Configuration
###############################################################################

proj_dir="${1:-/jukebox/krienen/marm_hmba_integration/261007_reprocess_and_recluster}"
script_dir="${proj_dir}/pipeline_rna/mmc"
reference_data_table="${script_dir}/mapmycells_reference.csv"
ortho_output_dir="${proj_dir}/data/rna/ortho"
out_dir="${proj_dir}/data/rna/mapmycells"

# Set to 1 if mapmycells_reference.csv has a header row.
csv_has_header=0

# Maximum number of MapMyCells array tasks that can run simultaneously.
max_concurrent_tasks=10

# Ortholog preparation resources (can override via environment variables).
ortho_mem="${ORTHO_MEM:-256G}"
ortho_cpus="${ORTHO_CPUS:-4}"

mmc_job_name="mmc_$(basename "${proj_dir}")"

###############################################################################
# Validate inputs
###############################################################################

cd "${script_dir}" || exit 1
mkdir -p logs "${out_dir}"

if [[ ! -s "${reference_data_table}" ]]; then
    echo "Error: Reference table does not exist or is empty:"
    echo "  ${reference_data_table}"
    exit 1
fi

if [[ ! -f "${script_dir}/2_run_mmc_row.sh" ]]; then
    echo "Error: MapMyCells row worker script was not found:"
    echo "  ${script_dir}/2_run_mmc_row.sh"
    exit 1
fi

###############################################################################
# Determine CSV row range
###############################################################################

# Physical CSV line numbers are used directly as SLURM_ARRAY_TASK_ID values.
total_lines="$(awk 'END {print NR}' "${reference_data_table}")"
first_row=$((csv_has_header + 1))

if (( total_lines < first_row )); then
    echo "Error: No data rows were found in:"
    echo "  ${reference_data_table}"
    exit 1
fi

echo "Project directory: ${proj_dir}"
echo "Reference table:   ${reference_data_table}"
echo "CSV data rows:     ${first_row}-${total_lines}"
echo "Slurm job name:    ${mmc_job_name}"
echo

###############################################################################
# Find rows already represented by active Slurm array tasks
###############################################################################

declare -A active_rows=()

echo "Checking for existing active jobs named ${mmc_job_name}..."

while IFS= read -r task_id; do
    [[ "${task_id}" =~ ^[0-9]+$ ]] && active_rows["${task_id}"]=1
done < <(
    squeue \
        --noheader \
        --array \
        --user="${USER}" \
        --name="${mmc_job_name}" \
        --states="PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED" \
        --format="%a"
)

if (( ${#active_rows[@]} > 0 )); then
    echo "Found ${#active_rows[@]} active array task(s)."
else
    echo "No matching active array tasks were found."
fi

echo

###############################################################################
# Identify rows to submit and the ortholog files they require
###############################################################################

rows_to_submit=()
declare -A required_ortho_files=()

for ((row = first_row; row <= total_lines; row++)); do
    line="$(sed -n "${row}p" "${reference_data_table}" | tr -d '\r')"

    if [[ -z "${line//[[:space:]]/}" ]]; then
        echo "CSV row ${row} is blank; skipping."
        continue
    fi

    if [[ -n "${active_rows[$row]+x}" ]]; then
        echo "CSV row ${row} already has an active job; skipping."
        continue
    fi

    IFS=',' read -r reference_name species _ _ gene_type <<< "${line}"

    reference_name="$(xargs <<< "${reference_name:-}")"
    species="$(xargs <<< "${species:-}")"
    gene_type="$(xargs <<< "${gene_type:-}")"

    if [[ -z "${reference_name}" ]]; then
        echo "CSV row ${row} has no reference name; skipping."
        continue
    fi

    result_csv="${out_dir}/${reference_name}_result_table_hierachical.csv"
    if [[ -s "${result_csv}" ]]; then
        echo "CSV row ${row} is already complete; skipping:"
        echo "  ${result_csv}"
        continue
    fi

    rows_to_submit+=("${row}")

    # Marmoset rows use the original h5ad directly; all others need an
    # ortholog query file determined by species and gene_type.
    if [[ "${species}" != "marmoset" ]]; then
        if [[ -z "${species}" || -z "${gene_type}" ]]; then
            echo "Error: CSV row ${row} requires ortholog mapping but species/gene_type is missing:"
            echo "  ${line}"
            exit 1
        fi

        required_ortho_files["${ortho_output_dir}/rna_${species}_ortho_${gene_type}.h5ad"]=1
    fi
done

echo

if (( ${#rows_to_submit[@]} == 0 )); then
    echo "Nothing to submit."
    echo "All CSV rows are blank, already active, or already complete."
    exit 0
fi

row_list="$(IFS=,; echo "${rows_to_submit[*]}")"
array_spec="${row_list}%${max_concurrent_tasks}"

echo "Rows that will be submitted:"
echo "  ${row_list}"
echo
echo "Slurm array specification:"
echo "  ${array_spec}"
echo

###############################################################################
# Submit ortholog preparation only if required outputs are missing
###############################################################################

orthologs_ready=1

if (( ${#required_ortho_files[@]} == 0 )); then
    echo "No submitted rows require ortholog query files."
else
    echo "Checking ortholog files required by submitted rows..."
    for ortho_file in "${!required_ortho_files[@]}"; do
        if [[ -s "${ortho_file}" ]]; then
            echo "Found:   ${ortho_file}"
        else
            echo "Missing: ${ortho_file}"
            orthologs_ready=0
        fi
    done
fi

echo

prep_job_id=""
dependency_args=()

if (( orthologs_ready )); then
    echo "All required ortholog files already exist and are nonempty."
    echo "Skipping ortholog preparation and submitting mapping directly."
else
    prep_script="${script_dir}/1_prepare_orthologs.sh"
    if [[ ! -f "${prep_script}" ]]; then
        echo "Error: Required ortholog outputs are incomplete and preparation script was not found:"
        echo "  ${prep_script}"
        exit 1
    fi

    echo "Ortholog resources: ${ortho_mem} RAM, ${ortho_cpus} CPUs"
    prep_job_id="$(sbatch --parsable \
        --mem="${ortho_mem}" \
        --cpus-per-task="${ortho_cpus}" \
        "${prep_script}" "${proj_dir}")"
    prep_job_id="${prep_job_id%%;*}"

    echo "Submitted ortholog preparation job: ${prep_job_id}"
    dependency_args=(--dependency="afterok:${prep_job_id}")
fi

###############################################################################
# Submit MapMyCells array
###############################################################################

mmc_job_id="$(
    sbatch \
        --parsable \
        --job-name="${mmc_job_name}" \
        "${dependency_args[@]}" \
        --array="${array_spec}" \
        "${script_dir}/2_run_mmc_row.sh" \
        "${proj_dir}"
)"
mmc_job_id="${mmc_job_id%%;*}"

echo "Submitted MapMyCells array job: ${mmc_job_id}"
echo
echo "Job name:                ${mmc_job_name}"
echo "Submitted CSV rows:       ${row_list}"
echo "Maximum concurrent tasks: ${max_concurrent_tasks}"
echo
echo "Monitor these jobs with:"
if [[ -n "${prep_job_id}" ]]; then
    echo "  squeue -j ${prep_job_id},${mmc_job_id}"
else
    echo "  squeue -j ${mmc_job_id}"
fi

echo
echo "Monitor all jobs for this project with:"
echo "  squeue --array --user=\"${USER}\" --name=\"${mmc_job_name}\""
