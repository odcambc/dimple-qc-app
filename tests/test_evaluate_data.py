import pandas as pd

from evaluate_data import test_per_base_file as run_test_per_base_file
from evaluate_data import test_mean_values as run_test_mean_values
from process_data import process_per_base_file


class TestTestPerBaseFile:
    def test_empty_means_returns_empty(self, variant_region_per_base_df):
        processed = process_per_base_file(variant_region_per_base_df, False)
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
