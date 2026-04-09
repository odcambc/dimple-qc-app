import numpy as np
import pandas as pd
import scipy.stats as stats

# Aggregation functions
aggregation_functions = {
    "n_total": ["mean", "std"],
    "reads_all": ["mean", "std"],
    "n_variants": ["mean", "std"],
    "variant_fraction": ["mean", "std"],
    "variant_fraction_percent": ["mean", "std"],
    "indel_fraction": ["mean", "std"],
    "indel_substitution_ratio": ["mean", "std"],
    "max_variant_base": ["mean", "std"],
    "entropy": ["mean", "std"],
    "effective_entropy": ["mean", "std"],
    "percent_of_max_entropy": ["mean", "std"],
    "expected_variant_codons": ["mean", "std"],
    "expected_ref_n": ["mean", "std"],
    "insertions": ["mean", "std"],
    "deletions": ["mean", "std"],
}

aggregation_columns = [
    "n_total",
    "reads_all",
    "n_variants",
    "variant_fraction",
    "variant_fraction_percent",
    "indel_fraction",
    "indel_substitution_ratio",
    "max_variant_base",
    "entropy",
    "effective_entropy",
    "percent_of_max_entropy",
    "expected_variant_codons",
    "expected_ref_n",
    "insertions",
    "deletions",
]


# Helper functions for data processing
def compute_n_variants(df: pd.DataFrame, bases: np.ndarray | None = None) -> pd.Series:
    if bases is None:
        bases = df[["A", "C", "G", "T"]].to_numpy()
    max_vals = bases.max(axis=1, keepdims=True)
    mask = bases == max_vals
    return (bases * (~mask)).sum(axis=1)


def compute_entropy(df: pd.DataFrame, bases: np.ndarray | None = None) -> np.ndarray:
    """Vectorized Shannon entropy across all rows."""
    if bases is None:
        bases = df[["A", "C", "G", "T"]].to_numpy()
    return stats.entropy(bases.astype(float), axis=1)


def compute_effective_entropy(df: pd.DataFrame, bases: np.ndarray | None = None) -> np.ndarray:
    """Vectorized effective entropy: entropy of non-reference bases only,
    with counts below 3 zeroed out."""
    if bases is None:
        bases = df[["A", "C", "G", "T"]].to_numpy()
    bases_f = bases.astype(float)
    max_vals = bases_f.max(axis=1, keepdims=True)
    # Zero out the reference (max) base
    filtered = np.where(bases_f == max_vals, 0, bases_f)
    # Zero out counts below 3
    filtered = np.where(filtered < 3, 0, filtered)
    return stats.entropy(filtered, axis=1)


def compute_max_non_ref_base(df: pd.DataFrame, bases: np.ndarray | None = None) -> np.ndarray:
    """Vectorized max non-reference base count."""
    if bases is None:
        bases = df[["A", "C", "G", "T"]].to_numpy()
    base_cols = np.array(["A", "C", "G", "T"])
    ref_vals = df["ref"].to_numpy()

    # Build a mask where each row's reference base column is True
    ref_mask = ref_vals[:, None] == base_cols[None, :]

    # Zero out the reference base, then take the max
    non_ref = np.where(ref_mask, 0, bases)
    return non_ref.max(axis=1)


def process_per_base_file(
    per_base_df: pd.DataFrame,
    reverse_complement: bool,
    origin_shift: int = 0,
) -> pd.DataFrame:
    if per_base_df.empty:
        return pd.DataFrame()

    per_base_df = per_base_df.copy()

    # Set the sequence length
    sequence_length = per_base_df["pos"].max()

    selected_codon_range = range(0, sequence_length // 3)

    subpool_codon_fraction = 1 / (len(selected_codon_range) + 1)

    per_base_df["codon_number"] = per_base_df["pos"] // 3

    per_base_df["is_selected"] = [True] * len(per_base_df)

    bases = per_base_df[["A", "C", "G", "T"]].to_numpy()

    per_base_df["n_variants"] = compute_n_variants(per_base_df, bases)

    per_base_df["n_indels"] = per_base_df[["insertions", "deletions"]].sum(axis=1)

    per_base_df["n_total"] = bases.sum(axis=1)

    # Calculate ratios while avoiding division by zero

    per_base_df["variant_fraction"] = (
        per_base_df["n_variants"] / per_base_df["reads_all"]
    ).replace([np.inf, -np.inf], np.nan)

    per_base_df["variant_fraction_percent"] = (
        (4 / 3) * per_base_df["variant_fraction"] / subpool_codon_fraction
    ).replace([np.inf, -np.inf], np.nan)

    per_base_df["indel_fraction"] = (
        per_base_df["n_indels"] / per_base_df["reads_all"]
    ).replace([np.inf, -np.inf], np.nan)

    per_base_df["indel_substitution_ratio"] = (
        per_base_df["n_indels"] / per_base_df["n_variants"]
    ).replace([np.inf, -np.inf], np.nan)

    # Calculate the counts of the most common non-reference base
    per_base_df["max_variant_base"] = compute_max_non_ref_base(per_base_df, bases)

    per_base_df["entropy"] = compute_entropy(per_base_df, bases)
    per_base_df["effective_entropy"] = compute_effective_entropy(per_base_df, bases)

    per_base_df["expected_variant_codons"] = (
        subpool_codon_fraction * per_base_df["n_total"]
    )

    per_base_df["expected_ref_n"] = per_base_df["n_total"] * (
        1 - subpool_codon_fraction
    )

    per_base_df["percent_of_max_entropy"] = per_base_df["effective_entropy"] / np.log(3)

    # Create empty columns for alignment to start
    per_base_df["aligned_ref"] = ["-"] * len(per_base_df)
    per_base_df["alignment_mismatch"] = [0] * len(per_base_df)

    # Apply origin shift (before reverse complement so they compose correctly)
    if origin_shift > 0:
        per_base_df["pos"] = ((per_base_df["pos"] - origin_shift - 1) % sequence_length) + 1
        per_base_df = per_base_df.sort_values("pos").reset_index(drop=True)

    if reverse_complement:
        # Reverse the index
        per_base_df = per_base_df.iloc[::-1]

        # Reverse the position values
        per_base_df["pos"] = sequence_length - per_base_df["pos"] + 1

        # Reverse complement the bases
        per_base_df[["A", "C", "G", "T"]] = per_base_df[["T", "G", "C", "A"]]

        # Reverse complement reference base
        per_base_df["ref"] = per_base_df["ref"].apply(
            lambda x: {"A": "T", "C": "G", "G": "C", "T": "A"}[x]
        )

    return per_base_df


def update_per_base_df(
    per_base_df: pd.DataFrame,
<<<<<<< HEAD
    selected_range_low: int,
    selected_range_high: int,
    parsed_reference: dict[str, str | list | None] | None,
    origin_shift: int = 0,
) -> None:
=======
    selected_range_list: list[tuple[int, int]],
) -> pd.DataFrame:
>>>>>>> main

    if per_base_df.empty:
        return per_base_df

<<<<<<< HEAD
    sequence_length = per_base_df["pos"].max()

    selected_positions = list(range(selected_range_low, selected_range_high))

    if parsed_reference:
        for feature in parsed_reference["features"]:
            if feature:
                if feature.type != "source":
                    start = int(feature.location.start)
                    end = int(feature.location.end)

                    # Apply origin shift to feature coordinates
                    if origin_shift > 0:
                        shifted_start = (start - origin_shift) % sequence_length
                        shifted_end = (end - origin_shift) % sequence_length

                        # Handle features that straddle the new origin
                        if shifted_start < shifted_end:
                            selected_positions += list(range(shifted_start, shifted_end))
                        else:
                            # Feature wraps around the origin
                            selected_positions += list(range(shifted_start, sequence_length))
                            selected_positions += list(range(0, shifted_end))
                    else:
                        selected_positions += list(range(start, end))
=======
    df = per_base_df.copy()

    # Build a vectorized boolean mask from the list of (start, end) ranges
    mask = pd.Series(False, index=df.index)
    total_selected = 0
    for start, end in selected_range_list:
        mask |= (df["pos"] >= start) & (df["pos"] < end)
        total_selected += end - start
>>>>>>> main

    subpool_codon_fraction = 3 / (total_selected + 1)

    df["is_selected"] = mask

    df["expected_variant_codons"] = (
        subpool_codon_fraction * df["n_total"]
    )

    df["variant_fraction_percent"] = (
        (4 / 3) * df["variant_fraction"] / subpool_codon_fraction
    ).replace([np.inf, -np.inf], np.nan)

    return df


def process_full_mean_values(processed_per_base_df: pd.DataFrame) -> pd.DataFrame:
    # Define the expected columns and multi-index structure
    multi_cols = pd.MultiIndex.from_product([aggregation_columns, ["mean", "std"]])

    # If input is empty, return a DataFrame with the expected structure filled with NaNs
    if processed_per_base_df.empty:
        return pd.DataFrame(
            np.nan,
            index=pd.Index(["selected", "unselected", "full"]),
            columns=multi_cols,
        )

    # Perform the groupby aggregation
    full_mean = processed_per_base_df.agg(aggregation_functions)

    # Now make sure all series are present
    full_mean = full_mean.reindex(aggregation_columns, fill_value=np.nan)

    return pd.DataFrame(full_mean)


def update_mean_values_per_base(
    processed_per_base_df: pd.DataFrame,
) -> pd.DataFrame:
    # Define the expected columns and multi-index structure
    multi_cols = pd.MultiIndex.from_product([aggregation_columns, ["mean", "std"]])

    # If input is empty, return a DataFrame with the expected structure filled with NaNs
    if processed_per_base_df.empty:
        return pd.DataFrame(
            np.nan,
            index=pd.Index(["selected", "unselected", "full"]),
            columns=multi_cols,
        )

    # Perform the groupby aggregation
    means = processed_per_base_df.groupby("is_selected").agg(aggregation_functions)
    means.columns = means.columns.map("_".join)

    # Flatten
    means.index = means.index.map({True: "selected", False: "unselected"})

    # Now make sure all series are present
    means = means.reindex(["selected", "unselected"], fill_value=np.nan)

    return pd.DataFrame(means)
