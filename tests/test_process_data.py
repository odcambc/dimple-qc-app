import numpy as np
import pandas as pd

from process_data import (
    compute_entropy,
    compute_effective_entropy,
    compute_max_non_ref_base,
    compute_n_variants,
    process_per_base_file,
    update_per_base_df,
    process_full_mean_values,
    update_mean_values_per_base,
)


class TestComputeNVariants:
    def test_basic(self, base_counts_df):
        result = compute_n_variants(base_counts_df)
        # Row 0: max=A(100), variants = 5+3+2 = 10
        # Row 1: max=C(200), variants = 5+4+6 = 15
        # Row 2: max=G(150), variants = 10+7+3 = 20
        # Row 3: max=T(180), variants = 8+12+6 = 26
        expected = np.array([10, 15, 20, 26])
        np.testing.assert_array_equal(result, expected)

    def test_all_same_base(self):
        df = pd.DataFrame({"A": [100], "C": [0], "G": [0], "T": [0]})
        result = compute_n_variants(df)
        assert result[0] == 0

    def test_even_distribution(self):
        df = pd.DataFrame({"A": [50], "C": [50], "G": [50], "T": [50]})
        result = compute_n_variants(df)
        # All bases are tied for max, so mask zeros all of them
        assert result[0] == 0


class TestComputeEntropy:
    def test_zero_entropy_single_base(self):
        df = pd.DataFrame({"A": [100], "C": [0], "G": [0], "T": [0]})
        result = compute_entropy(df)
        assert result[0] == 0.0

    def test_max_entropy_even_distribution(self):
        df = pd.DataFrame({"A": [100], "C": [100], "G": [100], "T": [100]})
        result = compute_entropy(df)
        np.testing.assert_almost_equal(result[0], np.log(4))

    def test_multiple_rows(self, base_counts_df):
        result = compute_entropy(base_counts_df)
        assert len(result) == 4
        # All rows have a dominant base, so entropy should be low but > 0
        assert all(0 < e < np.log(4) for e in result)


class TestComputeEffectiveEntropy:
    def test_single_base_returns_nan(self):
        # When all non-max bases are zero, entropy is undefined (nan)
        df = pd.DataFrame({"A": [100], "C": [0], "G": [0], "T": [0]})
        result = compute_effective_entropy(df)
        assert np.isnan(result[0])

    def test_counts_below_threshold_returns_nan(self):
        # Non-max bases all below 3 are zeroed → all-zero input → nan
        df = pd.DataFrame({"A": [100], "C": [2], "G": [1], "T": [2]})
        result = compute_effective_entropy(df)
        assert np.isnan(result[0])

    def test_counts_above_threshold(self):
        # Non-max bases above threshold should produce nonzero entropy
        df = pd.DataFrame({"A": [100], "C": [10], "G": [10], "T": [10]})
        result = compute_effective_entropy(df)
        # Three equal non-ref bases → entropy = ln(3)
        np.testing.assert_almost_equal(result[0], np.log(3))


class TestComputeMaxNonRefBase:
    def test_basic(self, base_counts_df):
        result = compute_max_non_ref_base(base_counts_df)
        # Row 0: ref=A, non-ref max = max(C=5, G=3, T=2) = 5
        # Row 1: ref=C, non-ref max = max(A=5, G=4, T=6) = 6
        # Row 2: ref=G, non-ref max = max(A=10, C=7, T=3) = 10
        # Row 3: ref=T, non-ref max = max(A=8, C=12, G=6) = 12
        np.testing.assert_array_equal(result, [5, 6, 10, 12])

    def test_all_zeros_except_ref(self):
        df = pd.DataFrame({"ref": ["A"], "A": [100], "C": [0], "G": [0], "T": [0]})
        result = compute_max_non_ref_base(df)
        assert result[0] == 0


class TestProcessPerBaseFile:
    def test_empty_input(self):
        result = process_per_base_file(pd.DataFrame(), False)
        assert result.empty

    def test_does_not_mutate_input(self, minimal_per_base_df):
        original = minimal_per_base_df.copy()
        process_per_base_file(minimal_per_base_df, False)
        pd.testing.assert_frame_equal(minimal_per_base_df, original)

    def test_output_columns(self, minimal_per_base_df):
        result = process_per_base_file(minimal_per_base_df, False)
        expected_cols = [
            "n_variants", "n_indels", "n_total", "variant_fraction",
            "indel_fraction", "entropy", "effective_entropy",
            "percent_of_max_entropy", "max_variant_base",
            "is_selected", "aligned_ref", "alignment_mismatch",
        ]
        for col in expected_cols:
            assert col in result.columns, f"Missing column: {col}"

    def test_all_selected_by_default(self, minimal_per_base_df):
        result = process_per_base_file(minimal_per_base_df, False)
        assert result["is_selected"].all()

    def test_reverse_complement_positions(self, minimal_per_base_df):
        result = process_per_base_file(minimal_per_base_df, True)
        # Positions should be reversed: max_pos - pos + 1
        max_pos = minimal_per_base_df["pos"].max()
        expected_positions = sorted(max_pos - minimal_per_base_df["pos"] + 1)
        assert sorted(result["pos"].tolist()) == expected_positions

    def test_reverse_complement_bases(self, minimal_per_base_df):
        result_fwd = process_per_base_file(minimal_per_base_df, False)
        result_rev = process_per_base_file(minimal_per_base_df, True)

        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

        # First position forward should match last position reverse (complemented)
        fwd_first_ref = result_fwd.iloc[0]["ref"]
        rev_last_ref = result_rev.iloc[-1]["ref"]
        assert rev_last_ref == complement[fwd_first_ref]


class TestUpdatePerBaseDf:
    def test_does_not_mutate_input(self, minimal_per_base_df):
        processed = process_per_base_file(minimal_per_base_df, False)
        original = processed.copy()
        update_per_base_df(processed, [(2, 4)])
        pd.testing.assert_frame_equal(processed, original)

    def test_returns_dataframe(self, minimal_per_base_df):
        processed = process_per_base_file(minimal_per_base_df, False)
        result = update_per_base_df(processed, [(2, 4)])
        assert isinstance(result, pd.DataFrame)

    def test_selection_range(self, minimal_per_base_df):
        processed = process_per_base_file(minimal_per_base_df, False)
        result = update_per_base_df(processed, [(2, 4)])
        # Positions 2, 3 should be selected; 1, 4, 5, 6 should not
        selected_pos = result[result["is_selected"]]["pos"].tolist()
        assert 2 in selected_pos
        assert 3 in selected_pos
        assert 1 not in selected_pos
        assert 4 not in selected_pos

    def test_multiple_ranges(self, minimal_per_base_df):
        processed = process_per_base_file(minimal_per_base_df, False)
        result = update_per_base_df(processed, [(1, 2), (5, 6)])
        selected_pos = result[result["is_selected"]]["pos"].tolist()
        assert 1 in selected_pos
        assert 5 in selected_pos
        assert 3 not in selected_pos

    def test_empty_input(self):
        result = update_per_base_df(pd.DataFrame(), [(1, 5)])
        assert result.empty


class TestProcessFullMeanValues:
    def test_empty_input(self):
        result = process_full_mean_values(pd.DataFrame())
        assert not result.empty  # Should return structured NaN DataFrame
        assert all(result.isna().all())

    def test_with_data(self, minimal_per_base_df):
        processed = process_per_base_file(minimal_per_base_df, False)
        result = process_full_mean_values(processed)
        assert not result.empty


class TestUpdateMeanValuesPerBase:
    def test_empty_input(self):
        result = update_mean_values_per_base(pd.DataFrame(), 0, 100)
        assert not result.empty  # Structured NaN DataFrame

    def test_with_data(self, minimal_per_base_df):
        processed = process_per_base_file(minimal_per_base_df, False)
        # Mark some as selected, some not
        updated = update_per_base_df(processed, [(2, 4)])
        result = update_mean_values_per_base(updated, 2, 4)
        assert "selected" in result.index
        assert "unselected" in result.index
