# #!/usr/bin/env bash
# 250625
PDIR=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/
OUTDIR=/jukebox/krienen/abc_atlas/241004_download_hmba_wg/hmba-marmoset-wg-802451596237-us-west-2/Subcortex/RNA/current
cp ${PDIR}/data/rna/rna_raw_v4_freeze.h5ad $OUTDIR
cp ${PDIR}/data/rna/metacell/subcortex/marm_Cluster_v4_metacell.h5ad $OUTDIR
