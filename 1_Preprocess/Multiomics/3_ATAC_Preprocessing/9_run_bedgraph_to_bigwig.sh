#!/usr/bin/env bash

INDIR=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/atac/bw/Subcortex_Group_v4/cpm_bedgraph
OUTDIR=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/atac/bw/Subcortex_Group_v4/cpm

[ -e $OUTDIR ] || mkdir -p $OUTDIR

BGTOBW=/jukebox/krienen/software/ucsc_source/bedGraphToBigWig
CHROMSIZES=/jukebox/krienen/genomic_annotations/marmoset_genome/mCalJa1.2.pat.X/custom_mito/chrom.sizes

for FILE in "$INDIR"/*.bedgraph; do
    base="${FILE##*/}"
    outname="${base%.bedgraph}.bw"
    "$BGTOBW" "$FILE" "$CHROMSIZES" "$OUTDIR/$outname"
    [[ -e "$OUTDIR/$outname" ]] || { echo "Failed to create $OUTDIR/$outname" >&2; exit 1; }
done
