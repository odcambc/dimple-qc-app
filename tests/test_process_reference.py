import pandas as pd
import pytest

from process_reference import align_ref_to_variants


class TestAlignRefToVariants:
    def test_empty_dataframe(self):
        result = align_ref_to_variants(pd.DataFrame(), "ACGT")
        assert result.empty

    def test_none_reference_returns_unchanged(self, minimal_per_base_df):
        result = align_ref_to_variants(minimal_per_base_df, None)
        assert "aligned_ref" not in result.columns

    def test_matching_reference_produces_columns(self, minimal_per_base_df):
        ref_seq = "".join(minimal_per_base_df["ref"])
        result = align_ref_to_variants(minimal_per_base_df.copy(), ref_seq)
        assert "aligned_ref" in result.columns
        assert "alignment_mismatch" in result.columns
        assert len(result) == len(minimal_per_base_df)

    def test_matching_reference_no_mismatches(self, minimal_per_base_df):
        ref_seq = "".join(minimal_per_base_df["ref"])
        result = align_ref_to_variants(minimal_per_base_df.copy(), ref_seq)
        assert result["alignment_mismatch"].sum() == 0

    def test_reverse_complement_reference_does_not_crash(self, minimal_per_base_df):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        ref_seq = "".join(minimal_per_base_df["ref"])
        rc_seq = "".join(complement.get(b, b) for b in reversed(ref_seq))
        result = align_ref_to_variants(minimal_per_base_df.copy(), rc_seq)
        assert "aligned_ref" in result.columns
        assert "alignment_mismatch" in result.columns
        assert len(result) == len(minimal_per_base_df)

    def test_short_region_aligns_to_reference_prefix(self):
        # A short sequenced window aligns to the matching prefix of a longer
        # reference; the trailing reference bases are insertions-in-reference
        # and are dropped, leaving a per-row mapping of the right length. The
        # matched bases are filled in rather than left as "-" placeholders.
        df = pd.DataFrame({
            "pos": [1, 2, 3],
            "ref": ["A", "C", "G"],
            "A": [10, 10, 10],
            "C": [5, 5, 5],
            "G": [5, 5, 5],
            "T": [1, 1, 1],
        })
        df["aligned_ref"] = ["-"] * len(df)
        df["alignment_mismatch"] = [0] * len(df)
        ref_seq = "ACGTACGTACGT"
        result = align_ref_to_variants(df.copy(), ref_seq)
        assert result["aligned_ref"].tolist() == ["A", "C", "G"]
        assert result["alignment_mismatch"].tolist() == [0, 0, 0]
