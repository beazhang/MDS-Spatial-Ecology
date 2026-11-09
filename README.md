# MDS bone marrow Xenium spatial transcriptomics -- figure source code

Code behind the figures of the MDS bone marrow Xenium spatial transcriptomics study:
Figures 1-5 and Extended Data Figures 1-12.

The repository is organised **by published figure**. Each figure has at most one Python
notebook and at most one R script, both named after the figure, and every panel inside
them carries its own `Figure <n><panel>` heading.

## Repository structure

```
notebooks/                   Python: one notebook per figure
  requirements.py            standardised imports, dataset paths, data_dir / media_dir
  Figure_01 .. 05.ipynb      main figures
  ED_Figure_01 .. 12.ipynb   Extended Data figures
  DUNE_script.ipynb          DUNE -- niche enrichment method
  SAND.ipynb                 SAND -- spatial cluster detection method

R_scripts/                   R: one script per figure
  _setup.R                   shared library() calls
  Figure_01 .. 05.R          main figures
  ED_Figure_01 .. 12.R       Extended Data figures
  data/                      the input table for each R panel, named per panel
    Fig1a.csv ... EDFig12i.csv

src/utils/                   helpers shared across notebooks
  pl/plot_ternary.py         ternary density plots
  pp/segmentation.py         segmentation-mask -> polygon helpers
  tl/deg.py                  PyDESeq2 pseudobulk DE, paired or unpaired
  tl/score_genes.py          gene signature scoring, z-scored onto .obs
  tl/runSAND_transcript.py   SAND adapted to transcript coordinates

pyproject.toml               Python dependencies
data/                        not distributed -- see Data below
media/                       figure output written by the notebooks (gitignored)
```

## How the Python and R code divide

Several panels appear in **both** a notebook and an R script. The division is consistent:
Python computes the underlying values -- differential expression, gene signature scores,
neighbourhood enrichment scores, per-sample mean expression -- and R renders the published
panel from the table Python wrote. Volcano plots, alluvials, and most of the per-sample
boxplots are therefore R; UMAPs, dotplots, heatmaps, and the ternary plots are Python.

Where a figure's heading says *"see DUNE_script.ipynb"* or *"see SAND.ipynb"*, the panel
comes from one of the two method notebooks rather than from that figure's own notebook.

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

Panels not listed above are not produced by code in this repository. They are spatial
image renders exported from Xenium Explorer, hand-drawn schematics and summary models,
wet-lab experiments (flow cytometry, transwell assays, gels, westerns, Sanger traces),
immunofluorescence quantification, probe-threshold FDR simulations, and splicing and
protein-structure analyses.

### Method notebooks

| Notebook | Method |
|---|---|
| [DUNE_script.ipynb](notebooks/DUNE_script.ipynb) | **DUNE** (Discovery by Unsupervised Niche Enrichment) -- niche enrichment scores and niche entropy, Leiden clustering of niches, cell-type composition per niche, and differentially enriched niches. Produces Figures 2b-e, ED 3f and ED 4a. |
| [SAND.ipynb](notebooks/SAND.ipynb) | **SAND** (Spatial Adjacency Network Detection) -- identifies local structures ("clusters") in spatial data from a radius graph plus a disjoint-set join. Used for immune aggregate identification in Figure 3. The transcript-level adaptation used for mutant/isoform transcript clusters is `src/utils/tl/runSAND_transcript.py`. |

`SAND.ipynb` demonstrates the method on a small toy AnnData rather than on a primary
dataset. That toy object is not distributed here; point its `read_h5ad` call at any
AnnData with spatial coordinates in `.obs`.

## Installation

```bash
pip install .
```

Versions are deliberately unpinned in `pyproject.toml`: the analysis ran over an extended
period and no single resolved environment reproduces every notebook. Pin to your own
working set.

The R scripts need `RColorBrewer`, `ggplot2`, `reshape2`, `ggpubr` and `ggalluvial`
(loaded by `_setup.R`), plus `ggbreak`, `ggrepel`, `ggfortify`, `tidyr`, `dplyr` and
`tibble`, which individual scripts load where they are used. With those installed, all
16 figure scripts run to completion on R 4.3 and draw every panel from the tables in
`R_scripts/data/`.

## Quickstart

Python:

1. Point `MDS_DATA_ROOT` at your copy of the data, or edit `DATA_ROOT` in
   `notebooks/requirements.py`.
2. Open the notebook for the figure you want, e.g. `notebooks/Figure_03.ipynb`.
3. Run the import cell and the dataset-load cell at the top, then the panel you want.
   Within a notebook the panels are independent of each other.

Every figure notebook has the same shape: cell 0 is `from requirements import *` plus any
`src/utils` helpers, cell 1 the figure title, cell 2 the dataset load, and from cell 3 on
one `## Figure 3c. <title>` heading per panel followed by that panel's code. A heading
reading *"DEG analysis for volcano plot: ..."* marks code that computes values for a panel
plotted in R; one reading *"see DUNE_script.ipynb"* is a pointer to a method notebook.

`notebooks/requirements.py` holds the consolidated imports, defines `DATA_ROOT` and the
three dataset constants, creates `media_dir` (`../media`), and puts `../src` on
`sys.path`, so notebooks import helpers as `from utils.tl.deg import runDEseq`.

R -- no setup needed, the input tables ship with the repository:

```r
setwd("R_scripts")
source("_setup.R")
source("Figure_03.R")
```

## Data

**The R panels are self-contained.** Every table the R scripts read lives in
`R_scripts/data/`, one file per panel, named for the panel it feeds -- `Fig3b.csv`,
`EDFig11l.csv`, and so on. The scripts reference them by relative path, so
`setwd("R_scripts")` is the only requirement. Four of the larger tables are gzipped, which
`read.delim()` handles transparently.

Three of those tables were reduced to the columns their figures actually use; the other 53
are verbatim copies of the pipeline output. The dropped columns are never read by the
plotting code, so the rendered panels are unchanged -- verified by rendering each affected
panel before and after and comparing the images pixel for pixel.

| File | Kept | Dropped |
|---|---|---|
| `EDFig3c.txt.gz` | `diagnosis`, `study_id`, `ct_general`, `ct_refined`, `tls`, `in_tls` | cell barcode, `sample_name` -- 462 MiB down to 18 MiB |
| `Fig1a.csv` | `manuscript_id`, `Diagnosis`, `IPSS-M`, `Age`, `IMPACT`, `SF3B1` | `Study_ID`, `ISSUES`, `slide_id`, `region`, `Sex` |
| `EDFig1b.csv` | `Type`, `IPSS-M`, `Timepoint`, `Patient`, `TP53 Mutant`, `SF Mutant`, `Age`, `Treatment` | `Study_ID`, `YrsSinceTp`, `Blastpct`, `Diagnosis WHO`, the three `IPSS-R`/`IPSS-M` risk-score columns, `Disease Status at BM`, `Status` |

`Age` in `Fig1a.csv` is read only to derive the `Under 40` / `40-60` / `Over 60` category
that the panel plots; the exact value never appears in a figure.

**The Python notebooks are not.** They read Xenium and CITE-seq AnnData objects, plus
intermediate `.csv`/`.txt` tables, that are too large to distribute here. Three primary
datasets are named in `notebooks/requirements.py` and resolved under `DATA_ROOT`:

| Constant | Dataset |
|---|---|
| `XEN_V1` | Xenium v1, resegmented, all cell types |
| `CITE` | CITE-seq: scRNA + scProtein + scTCR |
| `XEN_V2` | Xenium v2, resegmented (TCR and mutation probe panels) |

### Data availability

| | Accession |
|---|---|
| Processed data (the AnnData objects above) | Zenodo [10.5281/zenodo.24020398](https://doi.org/10.5281/zenodo.24020398) |
| Raw CITE-seq data | GEO [GSE328696](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328696) |
| Raw Xenium data | GEO [GSE326970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326970) |

To run the notebooks, download the processed objects from Zenodo and point `DATA_ROOT` (or
the `MDS_DATA_ROOT` environment variable) at the directory holding them.

## Things to know before running anything

- **The notebooks are not runnable end to end as a batch.** Each notebook loads only the
  dataset(s) its own panels need, but individual cells still depend on intermediate
  `.csv`/`.txt` tables written elsewhere in the pipeline.
- **Only the Python side still reads absolute paths.** Anything other than the three
  `DATA_ROOT` datasets is read by an absolute `/data1/...` path from inside a notebook cell
  and will need rewriting to run elsewhere. The R scripts no longer have this problem --
  all 61 of their read calls point into `R_scripts/data/`.
- **Python feeds R, but you do not have to rerun it.** For a panel listed under both
  languages, the table the R block reads is the output of the notebook section -- and that
  table ships in `R_scripts/data/`, so the R script reproduces the panel on its own. Rerun
  the notebook section only if you want to regenerate the table from the raw data.
- **Some helpers exist in several diverging versions and are deliberately inline.**
  `makeneighborgraph` and `joinsets` differ between `SAND.ipynb` (cell-level) and
  `src/utils/tl/runSAND_transcript.py` (transcript-level), and `p_to_star` is defined in
  the cells that use it. Only helpers with a single unambiguous definition were hoisted
  into `src/utils/`. Do not deduplicate the inline ones without checking which version
  each caller needs.
- **Several R blocks keep superseded commented-out variants** of the same plot, and
  commented-out `read.delim` lines for the other inputs the same block was run on (for
  example mutant vs control megakaryocytes in `Figure_05.R`, and the myeloid/erythroid
  variants in `ED_Figure_11.R`). Uncomment the input you want; the plotting code below it
  is shared, and all four alternative inputs ship in `R_scripts/data/` too.

## License

MIT -- see [`LICENSE`](LICENSE).
