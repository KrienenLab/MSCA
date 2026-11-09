#!/usr/bin/env python
"""SAMap mapping between marmoset and mouse THM/SC GABA neurons (Fig. S16D).

Inputs (in ./data, next to this script):
    marm_THM_SC_GABA.h5ad               marmoset cells, groups in obs["Subcortex_Group_v4"]
    mouse_abc_select_TH_SC_GABA.h5ad    mouse cells, subclasses in obs["subclass"]
Gene maps: ./maps_by_gene (see https://github.com/mdeber/Corrigan2024/tree/main/SAMap/pipeline
for how to generate this directory).
Outputs (in ./output): mapping score table (csv) and hierarchically clustered heatmap (pdf).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import seaborn as sns
from samalg import SAM
from samap.analysis import get_mapping_scores
from samap.mapping import SAMAP

# --------------------------------------------------
# Paths
# --------------------------------------------------
script_dir = Path(__file__).parent.resolve()
data_dir = script_dir / "data"
maps_dir = script_dir / "maps_by_gene"
outdir = script_dir / "output"
outdir.mkdir(parents=True, exist_ok=True)

# species ids: 'cj' = marmoset, 'mm' = mouse
h5ads = {
    "mm": data_dir / "mouse_abc_select_TH_SC_GABA.h5ad",
    "cj": data_dir / "marm_THM_SC_GABA.h5ad",
}
anno_keys = {"cj": "Subcortex_Group_v4", "mm": "subclass"}

# --------------------------------------------------
# Run SAMap
# --------------------------------------------------
# use pre-processed SAM objects (<name>_pr.h5ad) if they exist for both species
pr_h5ads = {k: p.with_name(p.stem + "_pr.h5ad") for k, p in h5ads.items()}
if all(p.is_file() for p in pr_h5ads.values()):
    sams = {}
    for k, p in pr_h5ads.items():
        sams[k] = SAM()
        sams[k].load_data(str(p))
else:
    sams = {k: str(p) for k, p in h5ads.items()}

sm = SAMAP(sams, f"{maps_dir}/", keys=anno_keys)
sm.run(pairwise=True)

# --------------------------------------------------
# Cell type mapping scores
# --------------------------------------------------
# n_top=0 uses all cells from each cluster (n_top > 0 uses only the top n aligning cells)
_, mt = get_mapping_scores(sm, anno_keys, n_top=0)

# >>> mt
#                                          mm_203 LGv-SPFp-SPFm Nkx2-2 Tcf7l2 Gaba  mm_205 SC-PAG Lef1 Emx2 Gaba  mm_206 SCm-PAG Cdh23 Gaba  mm_207 SCs Dmbx1 Gaba  mm_209 SCs Pax7 Nfia Gaba  mm_211 SC Tnnt1 Gli3 Gaba  mm_212 SCs Lef1 Gli3 Gaba  mm_213 SCsg Gabrr2 Gaba  cj_MB-THM GABA          cj_PG NKX2-2-TCF7L2 GABA  cj_SC DMBX1 GABA        cj_SC GLI3-EBF1 GABA   cj_SC-PAG LEF1-EMX2 GABA  cj_SC-PAG PAX7-CDH23 GABA  cj_SCde GLI3-PAX7-SGCG GABA  cj_SCsu PAX7-COL5A2 GABA  cj_SCsu PAX7-CXCL14 GABA  cj_SCzo GABRR2 GABA    cj_THM GABA             cj_THM-MB GABA
# mm_203 LGv-SPFp-SPFm Nkx2-2 Tcf7l2 Gaba  0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.002725966393995082    0.8506818853100618        0.021607245009999457    0.04195972400345528    0.011927990489917066      0.2135921631322366         0.002216237290524274         0.0023180142745266434     0.0019090285515795274     0.0011978807734879392  0.004537045623301149    0.004979425726642258
# mm_205 SC-PAG Lef1 Emx2 Gaba             0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.0011244515578070703   0.00012756284293988373    0.04315802969761069     0.023753139265898      0.7432931869073033        0.12134415455209373        0.08826622856742257          0.4367527260186533        0.005619062092614575      0.0008902471002719744  0.00012063054747729737  0.0010942649033654417
# mm_206 SCm-PAG Cdh23 Gaba                0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.00017280483722529818  0.0                       0.014240093263246942    0.01473064153558828    0.15976847811505074       0.5239618041420585         0.061035090964584135         0.003105335158645392      0.02450458421986783       0.0004066636692015476  2.203899539495708e-05   0.0005258824789694203
# mm_207 SCs Dmbx1 Gaba                    0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.004094003237592668    0.0                       0.4447042474676863      0.0017999288226701042  0.003849370007159234      0.2212600109287969         0.0012063477022540151        0.0033850631502983673     0.0008502820768658689     0.004501371266278387   0.011280094316549676    0.003594497707451739
# mm_209 SCs Pax7 Nfia Gaba                0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.0004968514455740635   0.0                       0.003939928322768824    0.018530310497335695   0.00048179477929219644    0.3339330164566463         0.0017934239920129411        0.0002758721635149083     0.728442269589384         0.00783136213514125    0.0005267574576511813   0.001250468195768476
# mm_211 SC Tnnt1 Gli3 Gaba                0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.00014847176334130389  0.0                       0.0005225718251032195   0.7811761661761085     0.0                       0.016202269506269282       0.0                          0.0                       0.0032001397148410174     0.0035656568664648495  0.0003769253871433848   0.00032138868906000196
# mm_212 SCs Lef1 Gli3 Gaba                0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.0                     0.0                       0.0001966578932311929   0.0005193247860622856  0.0011052932934035764     0.0038908156292132         0.6013010990177688           0.0005087176075721396     0.0009785169638833115     0.0                    0.0012327594664476387   8.279145023056019e-05
# mm_213 SCsg Gabrr2 Gaba                  0.0                                      0.0                           0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.6751409397116938      0.0                       0.00028308418292773076  0.0007253814545967964  1.2870038835186667e-05    6.812882594822118e-05      0.002671530978635626         0.0                       0.0001873842393395267     0.572735087674715      0.6888868351151498      0.6804496226592397
# cj_MB-THM GABA                           0.002725966393995082                     0.0011244515578070703         0.00017280483722529818     0.004094003237592668   0.0004968514455740635      0.00014847176334130389     0.0                        0.6751409397116938       0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_PG NKX2-2-TCF7L2 GABA                 0.8506818853100618                       0.00012756284293988373        0.0                        0.0                    0.0                        0.0                        0.0                        0.0                      0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SC DMBX1 GABA                         0.021607245009999457                     0.04315802969761069           0.014240093263246942       0.4447042474676863     0.003939928322768824       0.0005225718251032195      0.0001966578932311929      0.00028308418292773076   0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SC GLI3-EBF1 GABA                     0.04195972400345528                      0.023753139265898             0.01473064153558828        0.0017999288226701042  0.018530310497335695       0.7811761661761085         0.0005193247860622856      0.0007253814545967964    0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SC-PAG LEF1-EMX2 GABA                 0.011927990489917066                     0.7432931869073033            0.15976847811505074        0.003849370007159234   0.00048179477929219644     0.0                        0.0011052932934035764      1.2870038835186667e-05   0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SC-PAG PAX7-CDH23 GABA                0.2135921631322366                       0.12134415455209373           0.5239618041420585         0.2212600109287969     0.3339330164566463         0.016202269506269282       0.0038908156292132         6.812882594822118e-05    0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SCde GLI3-PAX7-SGCG GABA              0.002216237290524274                     0.08826622856742257           0.061035090964584135       0.0012063477022540151  0.0017934239920129411      0.0                        0.6013010990177688         0.002671530978635626     0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SCsu PAX7-COL5A2 GABA                 0.0023180142745266434                    0.4367527260186533            0.003105335158645392       0.0033850631502983673  0.0002758721635149083      0.0                        0.0005087176075721396      0.0                      0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SCsu PAX7-CXCL14 GABA                 0.0019090285515795274                    0.005619062092614575          0.02450458421986783        0.0008502820768658689  0.728442269589384          0.0032001397148410174      0.0009785169638833115      0.0001873842393395267    0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_SCzo GABRR2 GABA                      0.0011978807734879392                    0.0008902471002719744         0.0004066636692015476      0.004501371266278387   0.00783136213514125        0.0035656568664648495      0.0                        0.572735087674715        0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_THM GABA                              0.004537045623301149                     0.00012063054747729737        2.203899539495708e-05      0.011280094316549676   0.0005267574576511813      0.0003769253871433848      0.0012327594664476387      0.6888868351151498       0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0
# cj_THM-MB GABA                           0.004979425726642258                     0.0010942649033654417         0.0005258824789694203      0.003594497707451739   0.001250468195768476       0.00032138868906000196     8.279145023056019e-05      0.6804496226592397       0.0                     0.0                       0.0                     0.0                    0.0                       0.0                        0.0                          0.0                       0.0                       0.0                    0.0                     0.0

mt.to_csv(outdir / "mapping_table_THM_GABA.csv")

# --------------------------------------------------
# Hierarchical clustering of mapping scores
# --------------------------------------------------
# similarity -> distance (diagonal of the score matrix is 0, so force the distance diagonal to 0)
distance_matrix = 1 - mt.values
np.fill_diagonal(distance_matrix, 0)
linkage = sch.linkage(ssd.squareform(distance_matrix), method="average")
ordered_labels = mt.index[sch.dendrogram(linkage, no_plot=True)["leaves"]]

# plot with the diagonal set to 1 (self-similarity)
scores = mt.to_numpy(copy=True)
np.fill_diagonal(scores, 1)
mt_plot = pd.DataFrame(scores, index=mt.index, columns=mt.columns).loc[ordered_labels, ordered_labels]

fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(mt_plot, cmap="viridis", vmin=0, vmax=1, ax=ax)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.8)
    spine.set_color("black")
ax.collections[0].colorbar.outline.set_linewidth(0.8)
ax.set_title("SAMap")
fig.tight_layout()
fig.savefig(outdir / "mapping_table_THM_heatmap_hclust.pdf")
plt.close(fig)
