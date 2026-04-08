"""Helpers for reading Plasmidsaurus-style per-base CSV/TSV uploads."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from validation import expected_columns


def _has_required_columns(df: pd.DataFrame) -> bool:
    """Return True if ``df`` has all columns required for downstream processing."""
    if df.empty:
        return False
    return not (set(expected_columns) - set(df.columns))


def _strip_leading_junk_column(df: pd.DataFrame) -> pd.DataFrame:
    """Drop a leading empty/unused column (e.g. from CSV with leading comma)."""
    if "pos" in df.columns:
        return df
    if len(df.columns) > 1 and df.columns[1] == "pos":
        return df.iloc[:, 1:].copy()
    first = str(df.columns[0])
    if first in ("", "Unnamed: 0") or first.startswith("Unnamed"):
        out = df.iloc[:, 1:].copy()
        if "pos" in out.columns:
            return out
    return df


def read_per_base_table(path: str | Path) -> pd.DataFrame:
    """
    Read a per-base table from CSV or TSV.

    Tries auto-detected delimiter first, then tab, then comma.

    Args:
        path: Path to the uploaded file on disk.

    Returns:
        Parsed DataFrame with required columns, or empty if parsing fails.
    """
    path = Path(path)
    seps: list[str | None] = [None, "\t", ","]
    last_error: Exception | None = None

    for sep in seps:
        try:
            if sep is None:
                df = pd.read_csv(path, sep=None, engine="python", header=0)
            else:
                df = pd.read_csv(path, sep=sep, header=0)
        except Exception as exc:
            last_error = exc
            continue

        df = _strip_leading_junk_column(df)
        if _has_required_columns(df):
            return df

    if last_error is not None:
        raise last_error
    return pd.DataFrame()
