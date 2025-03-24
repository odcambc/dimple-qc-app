# validation.py

import pandas as pd
from shiny.express import ui

# Expected columns in the per-base file
expected_columns = [
    "pos",
    "ref",
    "reads_all",
    "matches",
    "mismatches",
    "deletions",
    "insertions",
    "low_conf",
    "A",
    "C",
    "G",
    "T",
]


def validate_per_base_file(per_base_file: pd.DataFrame) -> bool:
    """
    Validate that the uploaded per-base file has the required columns.

    Args:
        per_base_file: DataFrame containing sequencing data

    Returns:
        bool: True if file is valid, False otherwise
    """
    if per_base_file.empty:
        return False

    # Check that the file has the required columns defined in shared.py
    missing_cols = set(expected_columns) - set(per_base_file.columns)
    if missing_cols:
        ui.notification_show(
            f"Input per-base file is missing required data. Check format. Missing columns: {', '.join(missing_cols)}",
        )
        return False

    return True
