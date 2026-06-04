# validation.py

import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

from shared import test_cols


def test_per_base_file(
    processed_data: pd.DataFrame,
    selected_means_df: pd.DataFrame,
    full_mean_df: pd.DataFrame,
) -> pd.DataFrame:
    """Function to test the metrics and return a dataframe with summary pass/fail results"""

    # If the dataframes are empty, return an empty dataframe
    if selected_means_df.empty or full_mean_df.empty:
        return pd.DataFrame()

    return test_mean_values(processed_data, selected_means_df, full_mean_df)


def test_mean_values(
    processed_data: pd.DataFrame,
    selected_means_df: pd.DataFrame,
    full_mean_df: pd.DataFrame,
) -> pd.DataFrame:
    """Function to test the mean values and return a dataframe with summary pass/fail results"""

    # Run t-tests for all columns and collect results
    p_values = []
    t_stats = []
    for column in test_cols:
        selected_values = processed_data[processed_data["is_selected"]][column]
        unselected_values = processed_data[~processed_data["is_selected"]][column]

        if len(selected_values) < 2 or len(unselected_values) < 2:
            t_stats.append(np.nan)
            p_values.append(np.nan)
            continue

        t_test_result = stats.ttest_ind(
            selected_values,
            unselected_values,
            equal_var=False,
        )
        t_stats.append(t_test_result.statistic)
        p_values.append(t_test_result.pvalue)

    # Apply Benjamini-Hochberg FDR correction (filter NaN p-values first)
    p_arr = np.array(p_values)
    valid_mask = ~np.isnan(p_arr)
    corrected_p = np.full_like(p_arr, np.nan)
    rejected = np.full(len(p_arr), False)

    if valid_mask.any():
        rej, cor, _, _ = multipletests(p_arr[valid_mask], method="fdr_bh")
        corrected_p[valid_mask] = cor
        rejected[valid_mask] = rej

    rows = []
    for i, column in enumerate(test_cols):
        rows.append({
            "Metric": column,
            "Test": "Compare means",
            "Result": "Skip" if np.isnan(p_values[i]) else ("Pass" if rejected[i] else "Fail"),
            "p_value": p_values[i],
            "p_adjusted": corrected_p[i],
            "t_stat": t_stats[i],
        })

    return pd.DataFrame(rows).set_index("Metric")
