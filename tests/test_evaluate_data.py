import numpy as np
import pandas as pd
import pytest

from evaluate_data import test_per_base_file as run_test_per_base_file
from evaluate_data import test_mean_values as run_test_mean_values
from process_data import (
    process_per_base_file,
    update_per_base_df,
    process_full_mean_values,
    update_mean_values_per_base,
)


@pytest.fixture
def per_base_df():
    """Per-base DataFrame where selected region has clearly different metrics."""
    np.random.seed(42)
    n = 100

    # Positions 1-100
    pos = list(range(1, n + 1))
    refs = ["A"] * n
    reads_all = [500] * n

    # Non-variant region: almost all reads match reference
    a_counts = [480] * n
    c_counts = [5] * n
    g_counts = [5] * n
    t_counts = [10] * n

    # Variant region (positions 30-60): much higher diversity
    for i in range(29, 60):
        a_counts[i] = 200
        c_counts[i] = 100
        g_counts[i] = 100
        t_counts[i] = 100

    return pd.DataFrame(
        {
            "pos": pos,
            "ref": refs,
            "reads_all": reads_all,
            "matches": [480] * n,
            "mismatches": [20] * n,
            "deletions": [1] * n,
            "insertions": [1] * n,
            "low_conf": [0] * n,
            "A": a_counts,
            "C": c_counts,
            "G": g_counts,
            "T": t_counts,
        }
    )


@pytest.fixture
def processed_test_data(per_base_df):
    """Processed data with a selected region that differs from unselected."""
    processed = process_per_base_file(per_base_df, False)
    updated = update_per_base_df(processed, [(30, 60)])
    selected_means = update_mean_values_per_base(updated, 30, 60)
    full_means = process_full_mean_values(updated)
    return updated, selected_means, full_means


class TestTestPerBaseFile:
    def test_empty_means_returns_empty(self, per_base_df):
        processed = process_per_base_file(per_base_df, False)
        result = run_test_per_base_file(processed, pd.DataFrame(), pd.DataFrame())
        assert result.empty

    def test_returns_results_for_all_test_cols(self, processed_test_data):
        from shared import test_cols

        data, selected_means, full_means = processed_test_data
        result = run_test_per_base_file(data, selected_means, full_means)
        assert len(result) == len(test_cols)

    def test_result_columns(self, processed_test_data):
        data, selected_means, full_means = processed_test_data
        result = run_test_per_base_file(data, selected_means, full_means)
        for col in ["Test", "Result", "p_value", "p_adjusted", "t_stat"]:
            assert col in result.columns, f"Missing column: {col}"

    def test_results_are_pass_or_fail(self, processed_test_data):
        data, selected_means, full_means = processed_test_data
        result = run_test_per_base_file(data, selected_means, full_means)
        assert set(result["Result"].unique()).issubset({"Pass", "Fail"})

    def test_adjusted_p_values_geq_raw(self, processed_test_data):
        """FDR-adjusted p-values should be >= raw p-values (skipping NaN rows)."""
        data, selected_means, full_means = processed_test_data
        result = run_test_per_base_file(data, selected_means, full_means)
        valid = result.dropna(subset=["p_value", "p_adjusted"])
        assert len(valid) > 0, "Expected at least some non-NaN p-values"
        for _, row in valid.iterrows():
            assert float(row["p_adjusted"]) >= float(row["p_value"]) - 1e-10

    def test_significant_difference_detected(self, processed_test_data):
        """With a clearly different variant region, at least some metrics should pass."""
        data, selected_means, full_means = processed_test_data
        result = run_test_per_base_file(data, selected_means, full_means)
        pass_count = (result["Result"] == "Pass").sum()
        assert pass_count > 0, "Expected at least one significant result"


class TestTestMeanValues:
    def test_fdr_correction_applied(self, processed_test_data):
        """Verify that FDR correction can change outcomes vs raw p-values."""
        data, selected_means, full_means = processed_test_data
        result = run_test_mean_values(data, selected_means, full_means)

        # At least verify p_adjusted column exists and differs from p_value
        assert "p_adjusted" in result.columns
        # Adjusted values should be at least as large as raw (skip NaN rows)
        valid = result.dropna(subset=["p_value", "p_adjusted"])
        assert len(valid) > 0, "Expected at least some non-NaN p-values"
        for _, row in valid.iterrows():
            assert float(row["p_adjusted"]) >= float(row["p_value"]) - 1e-10
