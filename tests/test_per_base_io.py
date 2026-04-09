"""Tests for per-base file reading."""

from pathlib import Path

import pandas as pd
import pytest

from per_base_io import read_per_base_table
from validation import expected_columns


@pytest.fixture
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


class TestReadPerBaseTable:
    def test_reads_comma_separated_sample(self, repo_root: Path) -> None:
        path = repo_root / "per_base_df.csv"
        if not path.is_file():
            pytest.skip("per_base_df.csv not in repo")

        df = read_per_base_table(path)
        assert not df.empty
        assert set(expected_columns).issubset(set(df.columns))

    def test_reads_tsv_roundtrip(self, tmp_path: Path) -> None:
        df_in = pd.DataFrame(
            {
                "pos": [1, 2],
                "ref": ["A", "C"],
                "reads_all": [100, 100],
                "matches": [90, 95],
                "mismatches": [10, 5],
                "deletions": [0, 0],
                "insertions": [0, 0],
                "low_conf": [0, 0],
                "A": [90, 5],
                "C": [5, 90],
                "G": [3, 3],
                "T": [2, 2],
            }
        )
        tsv_path = tmp_path / "sample.tsv"
        df_in.to_csv(tsv_path, sep="\t", index=False)
        df_out = read_per_base_table(tsv_path)
        assert not df_out.empty
        assert set(expected_columns).issubset(set(df_out.columns))
        assert len(df_out) == 2

    def test_reads_csv_roundtrip(self, tmp_path: Path) -> None:
        df_in = pd.DataFrame(
            {
                "pos": [1, 2],
                "ref": ["A", "C"],
                "reads_all": [100, 100],
                "matches": [90, 95],
                "mismatches": [10, 5],
                "deletions": [0, 0],
                "insertions": [0, 0],
                "low_conf": [0, 0],
                "A": [90, 5],
                "C": [5, 90],
                "G": [3, 3],
                "T": [2, 2],
            }
        )
        csv_path = tmp_path / "sample.csv"
        df_in.to_csv(csv_path, index=False)
        df_out = read_per_base_table(csv_path)
        assert not df_out.empty
        assert set(expected_columns).issubset(set(df_out.columns))
