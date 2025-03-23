from Bio import pairwise2
import numpy as np
import pandas as pd
import scipy.stats as stats


# Helper functions for data processing
def effective_entropy(x: pd.Series) -> np.float64:
    row_list = [i for i in x if i != max(x)]

    # Set values below 2 to 0
    row_list = [0 if i < 3 else i for i in row_list]

    return np.float64(stats.entropy(row_list))


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

    per_base_df["n_variants"] = per_base_df[["A", "C", "G", "T"]].apply(
        n_variants, axis=1
    )

    per_base_df["n_indels"] = per_base_df[["insertions", "deletions"]].apply(
        n_observations, axis=1
    )
    per_base_df["n_total"] = per_base_df[["A", "C", "G", "T"]].apply(
        n_observations, axis=1
    )

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
    parsed_reference: dict[str, str | list | None],
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
                    print(f"Feature: {feature.type} {start} {end}")

    selected_codon_range = range(selected_range_low // 3, selected_range_high // 3)
    subpool_codon_fraction = 3 / (len(selected_positions) + 1)

    per_base_df["is_selected"] = per_base_df["pos"].isin(selected_positions)
    print(len(selected_positions))

    per_base_df["expected_variant_codons"] = (
        subpool_codon_fraction * per_base_df["n_total"]
    )

    per_base_df["variant_fraction_percent"] = (
        (4 / 3) * per_base_df["variant_fraction"] / subpool_codon_fraction
    ).replace([np.inf, -np.inf], np.nan)

    return


def update_mean_values_per_base(
    processed_per_base_df: pd.DataFrame, min_pos: int, max_pos: int
) -> pd.DataFrame:
    if processed_per_base_df.empty:
        return pd.DataFrame()

    means = processed_per_base_df.groupby("is_selected").agg(
        {
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
    )

    return means
