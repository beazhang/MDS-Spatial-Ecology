# MDS bone marrow Xenium spatial transcriptomics and scRNA-seq -- figure source code

Code behind the figures of "Ecologic determinants of disease and immunity in myelodysplastic syndromes":
Figures 1-5 and Extended Data Figures 1-12.

The repository is organised **by published figure**. Each figure has at most one Python
notebook and at most one R script, both named after the figure in the manuscript.

## Repository structure

```
notebooks/                   Python: one notebook per figure
  requirements.py            Standardised imports, dataset paths, data_dir
  Figure_01 .. 05.ipynb      Main figures
  ED_Figure_01 .. 12.ipynb   Extended Data figures
  DUNE_script.ipynb          DUNE -- niche enrichment method
  SAND.ipynb                 SAND -- spatial cluster detection method

R_scripts/                   R: one script per figure
  _setup.R                   Shared library() calls
  Figure_01 .. 05.R          Main figures
  ED_Figure_01 .. 12.R       Extended Data figures
  data/                      Input data table for each R panel, named per panel
    Fig1a.csv ... EDFig12i.csv

src/utils/                   helpers shared across notebooks
  pl/plot_ternary.py         ternary density plots
  tl/deg.py                  PyDESeq2 pseudobulk DE, paired or unpaired
  tl/score_genes.py          gene signature scoring, z-scored onto .obs
  tl/runSAND_transcript.py   SAND adapted to transcript coordinates

pyproject.toml               Python dependencies
data/                        not distributed -- see Data below
media/                       figure output written by the notebooks (gitignored)
```

## Contents

### Main figures

| Figure | Python notebook | Panels | R script | Panels |
|---|---|---|---|---|
| Figure 1 | [Figure_01.ipynb](notebooks/Figure_01.ipynb) | 1b (Xenium UMAP), 1b (CITE-seq UMAP), 1c | [Figure_01.R](R_scripts/Figure_01.R) | 1a, 1c |
| Figure 2 | [Figure_02.ipynb](notebooks/Figure_02.ipynb) | 2g; 2b-e via [DUNE](notebooks/DUNE_script.ipynb) | [Figure_02.R](R_scripts/Figure_02.R) | 2e, 2f, 2g |
| Figure 3 | [Figure_03.ipynb](notebooks/Figure_03.ipynb) | 3c, 3d, 3e, 3f (DE), 3g (DE); aggregates via [SAND](notebooks/SAND.ipynb) | [Figure_03.R](R_scripts/Figure_03.R) | 3b, 3f, 3g |
| Figure 4 | [Figure_04.ipynb](notebooks/Figure_04.ipynb) | 4a, 4b, 4c, 4g | [Figure_04.R](R_scripts/Figure_04.R) | 4h, 4i |
| Figure 5 | [Figure_05.ipynb](notebooks/Figure_05.ipynb) | 5b, 5h | [Figure_05.R](R_scripts/Figure_05.R) | 5a, 5b, 5c, 5e, 5h |

### Extended Data figures

| Figure | Python notebook | Panels | R script | Panels |
|---|---|---|---|---|
| ED Figure 1 | [ED_Figure_01.ipynb](notebooks/ED_Figure_01.ipynb) | 1c, 1d, 1e | [ED_Figure_01.R](R_scripts/ED_Figure_01.R) | 1a, 1b |
| ED Figure 2 | [ED_Figure_02.ipynb](notebooks/ED_Figure_02.ipynb) | 2a, 2d, 2e | [ED_Figure_02.R](R_scripts/ED_Figure_02.R) | 2e |
| ED Figure 3 | [ED_Figure_03.ipynb](notebooks/ED_Figure_03.ipynb) | 3a, 3e; 3f via [DUNE](notebooks/DUNE_script.ipynb) | [ED_Figure_03.R](R_scripts/ED_Figure_03.R) | 3b, 3c |
| ED Figure 4 | [ED_Figure_04.ipynb](notebooks/ED_Figure_04.ipynb) | 4b, 4c, 4e, 4j; 4a via [DUNE](notebooks/DUNE_script.ipynb) | [ED_Figure_04.R](R_scripts/ED_Figure_04.R) | 4b, 4d, 4g, 4h, 4i, 4k |
| ED Figure 5 | [ED_Figure_05.ipynb](notebooks/ED_Figure_05.ipynb) | 5a, 5d, 5e, 5i (DE) | [ED_Figure_05.R](R_scripts/ED_Figure_05.R) | 5a, 5g, 5h, 5i |
| ED Figure 6 | -- | -- | -- | -- |
| ED Figure 7 | [ED_Figure_07.ipynb](notebooks/ED_Figure_07.ipynb) | 7a, 7b, 7c | [ED_Figure_07.R](R_scripts/ED_Figure_07.R) | 7c |
| ED Figure 8 | [ED_Figure_08.ipynb](notebooks/ED_Figure_08.ipynb) | 8a, 8b | [ED_Figure_08.R](R_scripts/ED_Figure_08.R) | 8c, 8d |
| ED Figure 9 | [ED_Figure_09.ipynb](notebooks/ED_Figure_09.ipynb) | 9c, 9d, 9f | [ED_Figure_09.R](R_scripts/ED_Figure_09.R) | 9a, 9b, 9e |
| ED Figure 10 | [ED_Figure_10.ipynb](notebooks/ED_Figure_10.ipynb) | 10a, 10e | [ED_Figure_10.R](R_scripts/ED_Figure_10.R) | 10c, 10g, 10i, 10j, 10l, 10m |
| ED Figure 11 | [ED_Figure_11.ipynb](notebooks/ED_Figure_11.ipynb) | 11a, 11l-n (DE) | [ED_Figure_11.R](R_scripts/ED_Figure_11.R) | 11b, 11c, 11d, 11e, 11g, 11h, 11j, 11k, 11l, 11m, 11n |
| ED Figure 12 | [ED_Figure_12.ipynb](notebooks/ED_Figure_12.ipynb) | 12i (score) | [ED_Figure_12.R](R_scripts/ED_Figure_12.R) | 12i |

"(DE)" and "(score)" mark notebook sections that compute the values for a panel plotted in
R rather than drawing the panel themselves.


### Method notebooks

| Notebook | Method |
|---|---|
| [DUNE_script.ipynb](notebooks/DUNE_script.ipynb) | **DUNE** (Discovery by Unsupervised Niche Enrichment) -- Produces Figures 2b-e, ED 3f and ED 4a. |
| [SAND.ipynb](notebooks/SAND.ipynb) | **SAND** (Spatial Adjacency Network Detection) -- demarcates discrete local structures ("clusters") in spatial data. Used for immune aggregate identification in Figure 3. The transcript-level adaptation used for cytokine transcript clusters is `src/utils/tl/runSAND_transcript.py`. |


## Installation

```bash
pip install .
```

## Quickstart

Python:

1. Point `MDS_DATA_ROOT` at your copy of the data, or edit `DATA_ROOT` in
   `notebooks/requirements.py`.
2. Open the notebook for the figure you want, e.g. `notebooks/Figure_03.ipynb`.
3. Run the import cell and the dataset-load cell at the top, then the panel you want.
   Within a notebook the panels are independent of each other.

`notebooks/requirements.py` holds the consolidated imports, defines `DATA_ROOT` and the
three dataset constants, and puts `../src` on
`sys.path`, so notebooks can import helpers as `from utils.tl.deg import runDEseq`.

R -- no setup needed, the input tables ship with the repository:

```r
setwd("R_scripts")
source("_setup.R")
source("Figure_03.R")
```

## Data

**The R panels are self-contained.** Every data table required to generate R figures is located in
`R_scripts/data/`, which includes one file per panel, named for the panel it feeds -- `Fig3b.csv`,
`EDFig11l.csv`, etc. 

**The Python notebooks are not.** They read Xenium and CITE-seq AnnData objects, plus
intermediate `.csv`/`.txt` tables, that are too large to distribute here. Three primary
datasets are named in `notebooks/requirements.py` and resolved under `DATA_ROOT`:

| Constant | Dataset |
|---|---|
| `XEN_V1` | Xenium v1, resegmented (Core MDS Spatial Cohort) |
| `CITE` | CITE-seq: scRNA + scProtein + scTCR |
| `XEN_V2` | Xenium v2, resegmented (MDS Spatial Cohort with TCR and mutation probe panels) |

### Data availability

| | Accession |
|---|---|
| Processed data (the AnnData objects above) | Zenodo [10.5281/zenodo.24020398](https://doi.org/10.5281/zenodo.24020398) |
| Raw CITE-seq data | GEO [GSE328696](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328696) |
| Raw Xenium data | GEO [GSE326970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326970) |

To run the notebooks, download the processed objects from Zenodo and point `DATA_ROOT` (or
the `MDS_DATA_ROOT` environment variable) at the directory holding them.

## Things to know before running anything

- **Only the Python side still reads absolute paths.** Anything other than the three
  `DATA_ROOT` datasets is read by an absolute `/data1/...` path from inside a notebook cell
  and will need rewriting to run elsewhere. The R scripts no longer have this problem --
  all 61 of their read calls point into `R_scripts/data/`.
- **Python feeds R, but you do not have to rerun it.** For a panel listed under both
  languages in table, the R block reads the output of the python notebook which is included
  in `R_scripts/data/`, so the R script reproduces the panel on its own. Rerun
  the notebook section to regenerate the table from the raw data.
- **Several R blocks keep superseded commented-out variants** of the same plot, and
  commented-out `read.delim` lines for the other inputs the same block was run on (for
  example mutant vs control megakaryocytes in `Figure_05.R`, and the myeloid/erythroid
  variants in `ED_Figure_11.R`). Uncomment the input you want; the plotting code below it
  is shared, and all four alternative inputs ship in `R_scripts/data/` too.

## License

MIT -- see [`LICENSE`](LICENSE).
