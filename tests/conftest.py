"""Shared pytest fixtures for the test suite."""

from collections.abc import Generator
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

from process_data import (
    process_full_mean_values,
    process_per_base_file,
    update_mean_values_per_base,
    update_per_base_df,
)
from validation import expected_columns


@pytest.fixture
def patch_validation_ui() -> Generator[MagicMock, None, None]:
    """Patch Shiny ``ui.notification_show`` when validating without a session."""
    with patch("validation.ui.notification_show") as mock:
        yield mock


@pytest.fixture
def valid_df() -> pd.DataFrame:
    """DataFrame with all required per-base columns."""
    return pd.DataFrame({col: [1] for col in expected_columns})


@pytest.fixture
def base_counts_df() -> pd.DataFrame:
    """DataFrame with known base counts for testing helper functions."""
    return pd.DataFrame(
        {
            "ref": ["A", "C", "G", "T"],
            "A": [100, 5, 10, 8],
            "C": [5, 200, 7, 12],
            "G": [3, 4, 150, 6],
            "T": [2, 6, 3, 180],
        }
    )


@pytest.fixture
def minimal_per_base_df() -> pd.DataFrame:
    """Minimal valid per-base DataFrame matching the expected input format."""
    return pd.DataFrame(
        {
            "pos": [1, 2, 3, 4, 5, 6],
            "ref": ["A", "C", "G", "T", "A", "C"],
            "reads_all": [200, 200, 200, 200, 200, 200],
            "matches": [190, 195, 185, 192, 188, 196],
            "mismatches": [10, 5, 15, 8, 12, 4],
            "deletions": [1, 0, 2, 1, 0, 1],
            "insertions": [0, 1, 0, 0, 2, 0],
            "low_conf": [0, 0, 0, 0, 0, 0],
            "A": [180, 5, 10, 8, 170, 4],
            "C": [5, 185, 7, 12, 10, 186],
            "G": [3, 4, 175, 6, 8, 4],
            "T": [2, 6, 3, 170, 7, 6],
        }
    )


@pytest.fixture
def variant_region_per_base_df() -> pd.DataFrame:
    """Per-base DataFrame where selected region has clearly different metrics."""
    np.random.seed(42)
    n = 100

    pos = list(range(1, n + 1))
    refs = ["A"] * n
    reads_all = [500] * n

    a_counts = [480] * n
    c_counts = [5] * n
    g_counts = [5] * n
    t_counts = [10] * n

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
def processed_test_data(
    variant_region_per_base_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Processed data with a selected region that differs from unselected."""
    processed = process_per_base_file(variant_region_per_base_df, False)
    updated = update_per_base_df(processed, [(30, 60)])
    selected_means = update_mean_values_per_base(updated, 30, 60)
    full_means = process_full_mean_values(updated)
    return updated, selected_means, full_means
