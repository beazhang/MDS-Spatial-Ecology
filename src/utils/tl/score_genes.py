"""
Gene signature scoring.

Scales a copy of the AnnData, scores each signature in `markers_dict` with
`sc.tl.score_genes`, then writes z-scored `zscore_<name>` columns back onto the
original object's `.obs`.

Used by: Figure 5h, Extended Data Figure 4b, Extended Data Figure 12i.
"""
import pandas as pd
import scanpy as sc
from copy import deepcopy
from scipy.stats import zscore


def score_genes(markers_dict, adata):
    adata_tmp = deepcopy(adata)
    sc.pp.scale(adata_tmp)
    
    for m in markers_dict.keys():
        sc.tl.score_genes(adata_tmp, gene_list = set(markers_dict[m]),score_name = 'zscore_' + m, use_raw = False)

    zscore_names = adata_tmp.obs.columns[adata_tmp.obs.columns.str.contains('^zscore_')]
    adata.obs = adata.obs.drop(adata.obs.columns[adata.obs.columns.str.contains('^zscore_')], axis = 1)
    adata.obs = pd.concat([adata.obs, adata_tmp.obs.loc[:,zscore_names].apply(zscore)], axis=1)
