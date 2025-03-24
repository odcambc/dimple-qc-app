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
def effective_entropy(x: pd.Series) -> np.float64:
    max_val = x.max()

    filtered = x[x != max_val].copy()

    filtered[filtered < 3] = 0

    return np.float64(stats.entropy(filtered))


def compute_n_variants(df: pd.DataFrame) -> pd.Series:
    bases = df[["A", "C", "G", "T"]].to_numpy()

    max_vals = bases.max(axis=1, keepdims=True)

    # Try using a mask instead of a loop
    mask = bases == max_vals

    return (bases * (~mask)).sum(axis=1)


def n_variants(x: pd.Series) -> int:
    return sum([i for i in x if i != max(x)])


def n_observations(x: pd.Series) -> int:
    return sum(x)


def max_non_ref_base(x):
    ref = x["ref"]
    non_ref_bases = [i for i in x[["A", "C", "G", "T"]] if i != x[ref]]
    if non_ref_bases:
        return max(non_ref_bases)
    return 0


def process_per_base_file(
    per_base_df: pd.DataFrame,
    reverse_complement: bool,
) -> pd.DataFrame:
    if per_base_df.empty:
        return pd.DataFrame()

    # Set the sequence length
    sequence_length = per_base_df["pos"].max()

    selected_codon_range = range(0, sequence_length // 3)

    subpool_codon_fraction = 1 / (len(selected_codon_range) + 1)

    per_base_df["codon_number"] = per_base_df["pos"] // 3

    per_base_df["is_selected"] = [True] * len(per_base_df)

    per_base_df["n_variants"] = compute_n_variants(per_base_df)

    per_base_df["n_indels"] = per_base_df[["insertions", "deletions"]].sum(axis=1)

    per_base_df["n_total"] = per_base_df[["A", "C", "G", "T"]].sum(axis=1)

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

    per_base_df["max_variant_base"] = per_base_df[["ref", "A", "C", "G", "T"]].apply(
        max_non_ref_base, axis=1
    )

    per_base_df["entropy"] = per_base_df[["A", "C", "G", "T"]].apply(
        lambda x: stats.entropy(x), axis=1
    )
    per_base_df["effective_entropy"] = per_base_df[["A", "C", "G", "T"]].apply(
        effective_entropy, axis=1
    )

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
    selected_range_low: int,
    selected_range_high: int,
    parsed_reference: dict[str, str | list | None] | None,
) -> None:

    if per_base_df.empty:
        return

    selected_positions = list(range(selected_range_low, selected_range_high))

    if parsed_reference:
        for feature in parsed_reference["features"]:
            if feature:
                if feature.type != "source":
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    selected_positions += list(range(start, end))

    selected_codon_range = range(selected_range_low // 3, selected_range_high // 3)
    subpool_codon_fraction = 3 / (len(selected_positions) + 1)

    per_base_df["is_selected"] = per_base_df["pos"].isin(selected_positions)

    per_base_df["expected_variant_codons"] = (
        subpool_codon_fraction * per_base_df["n_total"]
    )

    per_base_df["variant_fraction_percent"] = (
        (4 / 3) * per_base_df["variant_fraction"] / subpool_codon_fraction
    ).replace([np.inf, -np.inf], np.nan)

    return


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
    processed_per_base_df: pd.DataFrame, min_pos: int, max_pos: int
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
