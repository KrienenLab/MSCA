#!/usr/bin/env bash
#SBATCH -J mmc_ortho
#SBATCH -o logs/mmc_ortho_%j.log
#SBATCH -e logs/mmc_ortho_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4

module purge
module load anacondapy/2024.02
module load cudatoolkit/12.4.1
conda activate /scratch/sdan/tools/cell_type_mapper

proj_dir="${1:?Usage: $0 PROJ_DIR}"
script_dir="${proj_dir}/pipeline_rna/mmc"
ortho_table_dir="/jukebox/krienen/genomic_annotations/ortholog_tables"

data="${proj_dir}/data/rna/rna_raw_combined.h5ad"
ortho_output_dir="${proj_dir}/data/rna/ortho"
out_dir="${proj_dir}/data/rna/mapmycells"

mkdir -p "${ortho_output_dir}" "${out_dir}" logs

echo "[$(date)] Preparing ortholog query files"
echo "Data path: ${data}"
echo "Ortholog output: ${ortho_output_dir}"

for species in mouse human; do
    for gene_type in gene_name gene_id; do
        output_path="${ortho_output_dir}/rna_${species}_ortho_${gene_type}.h5ad"

        if [[ -s "${output_path}" ]]; then
            echo "Already exists; skipping: ${output_path}"
            continue
        fi

        echo "Converting ${species}, ${gene_type}"
        python "${script_dir}/1_marm_h5ad_to_orthologs.py" \
            "${ortho_table_dir}" \
            "${species}" \
            "${data}" \
            "${ortho_output_dir}" \
            "${gene_type}"
    done
done

echo "[$(date)] Ortholog preparation complete"
