#!/usr/bin/env bash
PDIR=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster
# [ -e ${PDIR}/data/rna/rna_merged.rds ] && rm ${PDIR}/data/rna/rna_merged.rds
[ -e ${PDIR}/data/rna/rna_merged_nojoin.rds ] && rm ${PDIR}/data/rna/rna_merged_nojoin.rds
[ -e ${PDIR}/data/rna/rna_merged_chunkjoin.rds ] && rm ${PDIR}/data/rna/rna_merged_chunkjoin.rds
rm ${PDIR}/data/rna/rna_counts_*.mtx
