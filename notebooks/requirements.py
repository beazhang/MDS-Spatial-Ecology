"""
Standardised imports and paths used throughout the figure notebooks for the MDS bone
marrow Xenium study.

Mirrors the `notebooks/requirements.py` pattern of Moorman et al., Nature 2024
(dpeerlab/progressive-plasticity-crc-metastasis).

Usage, from a notebook in this directory:

    from requirements import *

Paths
-----
`DATA_ROOT` is the one place to repoint if you are not on the originating cluster.
The three primary dataset constants below derive from it. Note that individual cells
also read intermediate .csv/.txt tables by absolute path; those were left inline and
verbatim, and are not covered by DATA_ROOT.

The three primary objects are available from Zenodo, doi:10.5281/zenodo.24020398.
Download them and point DATA_ROOT (or the MDS_DATA_ROOT environment variable) at the
directory that holds them. Raw data: GEO GSE328696 (CITE-seq), GSE326970 (Xenium).
"""

# ---------------------------------------------------------------------------
# Consolidated imports: the union of what the figure notebooks use.
# ---------------------------------------------------------------------------
from anndata import AnnData
import anndata as ad
import ast
import cmasher as cmr
from collections import Counter
import colorsys
from copy import deepcopy
import decoupler as dc
import gzip
from itertools import chain
import itertools
import json
import matplotlib
import squidpy as sq
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import random
import scanpy as sc
from scipy import stats
import scipy
from scipy.stats import zscore
import seaborn as sns; sns.set_theme(color_codes=True)
from sklearn.neighbors import NearestNeighbors
import statistics
import sys
#import ternary
from typing import List
import warnings
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Local modules (src/utils)
# ---------------------------------------------------------------------------
from pathlib import Path

module_path = os.path.abspath('../src')
if module_path not in sys.path:
    sys.path.append(module_path)

# ---------------------------------------------------------------------------
# Standard directories
# ---------------------------------------------------------------------------
# Repoint DATA_ROOT (or set the MDS_DATA_ROOT environment variable) to run elsewhere.
DATA_ROOT = Path(os.environ.get(
    'MDS_DATA_ROOT',
    '/data1/greenbab/users/zhangb2',
))

data_dir = DATA_ROOT
media_dir = Path(f"{os.getcwd()}/../media")
media_dir.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Primary datasets. Each notebook loads only the one(s) its own panels need.
# ---------------------------------------------------------------------------
XEN_V1     = f'{data_dir}/xenium_files/MDS_final_h5ad/MDS_xenium_v1.h5ad'
CITE       = f'{data_dir}/xenium_files/MDS_final_h5ad/MDS_CITE.h5ad'
XEN_V2     = f'{data_dir}/xenium_files/MDS_final_h5ad/MDS_xenium_v2.h5ad'

