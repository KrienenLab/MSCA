#!/usr/bin/env bash
PDIR=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster
cd $PDIR
JOB_CSV=${PDIR}/qc/fragment_jobs_to_run.csv
[ -e slurm_logs_atac ] || mkdir slurm_logs_atac
cd slurm_logs_atac
tail -n +2 "$JOB_CSV" | while IFS=',' read -r sample_name fragment_path barcoded_name; do
    sbatch ${PDIR}/pipeline_atac/make_snap_fragments_batch_script.sh "$sample_name"
done
