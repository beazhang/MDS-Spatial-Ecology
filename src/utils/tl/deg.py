"""
PyDESeq2-based pseudobulk DE, paired (patient-blocked) or unpaired.
"""

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def get_pseudobulk_by_patient_condition(
    adata,
    patient_col="patient_id",
    condition_col="condition",
    min_cells_per_sample=10,
):
    """
    Aggregates counts per (patient x condition) pair, filtering out pairs
    with fewer than `min_cells_per_sample` cells. This is the PAIRED path:
    a patient contributing cells to both conditions gets two separate
    pseudobulk rows (one per condition), correctly preserving the paired
    within-patient structure for a patient-blocked design.
    """
    sub_adata = adata.copy()

    # Unique sample grouping key combining patient and condition
    sub_adata.obs["sample_group"] = (
        sub_adata.obs[patient_col].astype(str)
        + "_"
        + sub_adata.obs[condition_col].astype(str)
    )

    # Filter out patient-condition pairs with too few cells
    group_cell_counts = sub_adata.obs["sample_group"].value_counts()
    valid_groups = group_cell_counts[group_cell_counts >= min_cells_per_sample].index

    print(group_cell_counts)
    print(
        f"\nRetaining {len(valid_groups)} / {len(group_cell_counts)} "
        f"patient-condition pairs with >= {min_cells_per_sample} cells."
    )

    sub_adata = sub_adata[sub_adata.obs["sample_group"].isin(valid_groups)].copy()

    # Extract raw counts and sum per patient-condition pair
    counts = sub_adata.layers["raw"]
    if hasattr(counts, "toarray"):
        counts = counts.toarray()

    df_counts = pd.DataFrame(counts, index=sub_adata.obs_names, columns=sub_adata.var_names)
    df_counts["sample_group"] = sub_adata.obs["sample_group"].values

    pb_counts = df_counts.groupby("sample_group").sum()

    pb_metadata = (
        sub_adata.obs[["sample_group", patient_col, condition_col]]
        .drop_duplicates()
        .set_index("sample_group")
    )
    pb_metadata = pb_metadata.loc[pb_counts.index]

    # Sanity check for downstream saturation warning in runDEseq()
    n_unique_patients = pb_metadata[patient_col].nunique()
    n_rows = len(pb_metadata)
    print(
        f"{n_unique_patients} unique patients across {n_rows} pseudobulk rows "
        f"(patient-blocked design will spend ~{n_unique_patients} degrees of "
        f"freedom on patient identity alone)."
    )

    return pb_counts, pb_metadata


def get_unpaired_pseudobulk(
    adata,
    patient_col="study_id",
    condition_col="condition",
    min_cells_per_sample=10,
):
    """
    Aggregates ALL of a patient's cells into a single pseudobulk row,
    regardless of condition. This is only valid when each patient has
    exactly one condition — a patient contributing cells to BOTH
    conditions has no well-defined single condition label, so such
    patients are detected and dropped (with a warning) rather than
    silently mislabeled by whichever condition happened to appear first
    in the data (the bug in the original version).

    If your data has patients spanning both conditions, use
    get_pseudobulk_by_patient_condition (the paired path) instead — this
    function's one-row-per-patient premise doesn't apply to them.
    """
    sub_adata = adata.copy()

    # --- Fix: detect and drop patients with >1 distinct condition ---
    conditions_per_patient = sub_adata.obs.groupby(patient_col)[condition_col].nunique()
    mixed_patients = conditions_per_patient[conditions_per_patient > 1].index.tolist()
    if mixed_patients:
        print(
            f"WARNING: {len(mixed_patients)} patient(s) have cells in more than "
            f"one '{condition_col}' — dropping them, since get_unpaired_pseudobulk "
            f"aggregates one row per patient and their condition label would be "
            f"undefined: {mixed_patients}"
        )
        sub_adata = sub_adata[~sub_adata.obs[patient_col].isin(mixed_patients)].copy()

    # Filter patients with too few cells
    patient_counts = sub_adata.obs[patient_col].value_counts()
    valid_patients = patient_counts[patient_counts >= min_cells_per_sample].index
    n_dropped_low_cell = sub_adata.obs[patient_col].nunique() - len(valid_patients)
    if n_dropped_low_cell > 0:
        print(
            f"Dropping {n_dropped_low_cell} patient(s) with < {min_cells_per_sample} cells."
        )

    sub_adata = sub_adata[sub_adata.obs[patient_col].isin(valid_patients)].copy()

    # Aggregate all counts directly by patient
    counts = sub_adata.layers["raw"]
    if hasattr(counts, "toarray"):
        counts = counts.toarray()

    df_counts = pd.DataFrame(counts, index=sub_adata.obs_names, columns=sub_adata.var_names)
    df_counts[patient_col] = sub_adata.obs[patient_col].values

    pb_counts = df_counts.groupby(patient_col).sum()

    # Now safe: every remaining patient has exactly one condition
    pb_metadata = (
        sub_adata.obs[[patient_col, condition_col]].drop_duplicates().set_index(patient_col)
    )
    pb_metadata = pb_metadata.loc[pb_counts.index]

    return pb_counts, pb_metadata


def runDEseq(
    adata,
    patient_col,
    condition_col,
    cond_1,
    cond_2,
    paired=True,
    min_cells_per_sample=10,
    saturation_warn_ratio=0.7,
    gene_filter = 0.5,
):
    """
    Run PyDESeq2 DE, either paired (patient-blocked, design_factors =
    [patient_col, condition_col]) or unpaired (design_factors =
    [condition_col] only, one row per patient).

    saturation_warn_ratio : float
        If paired=True and (n_unique_patients / n_pseudobulk_rows) exceeds
        this, print a warning — the patient-blocking design is close to
        saturated (nearly one parameter per row), so residual degrees of
        freedom for the condition Wald test may be too thin to trust,
        independent of how well-specified the NB GLM otherwise is.
    """

    ### first filter out samples with no pairs
    if paired:
        df_tmp = pd.crosstab(index=adata.obs[patient_col], columns=adata.obs[condition_col])
        removeidx = [x for x in df_tmp[(df_tmp[cond_1] < min_cells_per_sample) | (df_tmp[cond_2] < min_cells_per_sample)].index]
        print('Removing ' + str(len(removeidx)) + ' unpaired samples')
        adata = adata[~adata.obs[patient_col].isin(removeidx)]
    # ==========================================
    # 1. Run Aggregation
    # ==========================================
    if paired:
        pb_counts, pb_metadata = get_pseudobulk_by_patient_condition(
            adata,
            patient_col=patient_col,
            condition_col=condition_col,
            min_cells_per_sample=min_cells_per_sample,
        )

    else:
        pb_counts, pb_metadata = get_unpaired_pseudobulk(
            adata,
            patient_col=patient_col,
            condition_col=condition_col,
            min_cells_per_sample=min_cells_per_sample,
        )

    # Convert to integer counts and drop low-expression genes
    pb_counts = pb_counts.astype(int)
    min_nonzero = int(gene_filter * pb_counts.shape[0])
    pb_counts = pb_counts.loc[:, (pb_counts > 0).sum(axis=0) >= min_nonzero]
    
    print((pb_counts > 0).sum(axis=0) )
    print(pb_counts.shape)
    print(pb_metadata[condition_col].value_counts())

    # ==========================================
    # 2. PyDESeq2 Model Setup
    # ==========================================
    # Include patient_id FIRST to block baseline patient variance,
    # and condition LAST to make it the target test variable.
    if paired:
        dds = DeseqDataSet(
            counts=pb_counts,
            metadata=pb_metadata,
            design_factors=[patient_col, condition_col],
            refit_cooks=True,
           # n_cpus=4,
        )
    else:
        dds = DeseqDataSet(
            counts=pb_counts,
            metadata=pb_metadata,
            design_factors=[condition_col],
            refit_cooks=True,
           # n_cpus=4,
        )

    dds.deseq2()

    # ==========================================
    # 3. Extract Differential Expression
    # ==========================================
    stat_res = DeseqStats(dds, contrast=[dds.design_factors[-1], cond_2, cond_1])
    stat_res.summary()
    results_df = stat_res.results_df

    return results_df
