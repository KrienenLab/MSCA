#!/usr/bin/env python
#
# [!] Note the 'NT' columns is surprisingly empty for the neurons, so will skip for now
import pandas as pd
path_in = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/metadata/subcortex_anno_table.csv'
path_out = '/jukebox/krienen/abc_atlas/241004_download_hmba_wg/hmba-marmoset-wg-802451596237-us-west-2/Subcortex/250625_subcortex_anno_table.csv'

anno_table = pd.read_csv(path_in)
# anno_table = anno_table.loc[:, ['marm_cluster_id', 'Subcortex_Class_v4', 'Subcortex_Group_v4', 'NT']]
anno_table = anno_table.loc[:, ['marm_cluster_id', 'Subcortex_Class_v4', 'Subcortex_Group_v4']]
anno_table = anno_table.rename(columns={'marm_cluster_id': 'Cluster_v4'})
# anno_table['NT'] = anno_table['NT'].fillna('NN')
anno_table.to_csv(path_out, index=False)
