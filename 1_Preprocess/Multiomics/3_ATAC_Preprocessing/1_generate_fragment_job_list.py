#!/usr/bin/env python3
#
# This script generates a list of samples/fragment files to use for
# downstream ATAC analysis. This list is used as input for creating
# the snapatac2 h5ad files with the fragments.
import pandas as pd
from pathlib import Path

# Setup
pdir = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_library_qc_table = Path(pdir, 'qc/library_qc_pass_fail_table.csv')
path_added_fails = Path(pdir, 'qc/additional_qc_fails.csv')

library_qc = pd.read_csv(path_library_qc_table).set_index('sample_name')
added_fails = pd.read_csv(path_added_fails).set_index('sample_name')
library_qc = library_qc.loc[library_qc['all_qc_pass'] == 'PASS', :]
library_qc = library_qc.loc[~library_qc.index.isin(added_fails.index), :]

# Add fragment path
library_qc['fragment_path'] = [
    str(Path(adir, sname, sname, 'outs/atac_fragments.tsv.gz')) 
    for sname, adir in zip(library_qc.index, library_qc['align_dir'])
]

# Save job list
out_csv = Path(pdir, 'qc/fragment_jobs_to_run.csv')
library_qc[['fragment_path', 'barcoded_cell_sample_name']].to_csv(out_csv, index_label='sample_name')
print(f"Wrote job list to {out_csv}")
