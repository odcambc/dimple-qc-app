from pathlib import Path

app_dir = Path(__file__).parent
# example_df = pd.read_csv(app_dir / "FKYSRV_1_PXR2.tsv")
example_per_base_tsv = app_dir / "FKYSRV_1_PXR2.tsv"
example_reference_fasta = app_dir / "FKYSRV_1_PXR2.fasta"

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

# Define the columns to show in the plot
column_names_dict = {
    "entropy": "Entropy",
    "effective_entropy": "Effective entropy",
    "percent_of_max_entropy": "Entropy % max",
    "n_total": "Total reads",
    "n_variants": "Variant reads",
    "variant_fraction": "Variant fraction",
    "variant_fraction_percent": "Variant fraction % expected",
    "insertions": "Insertion count",
    "deletions": "Deletion count",
    "indel_fraction": "Indel fraction",
    "indel_substitution_ratio": "Indel to substitution ratio",
    "alignment_mismatch": "red",
    "max_variant_base": "Max counts non-ref base",
}

column_colors_dict = {
    "entropy": "#009E73",
    "effective_entropy": "#56B4E9",
    "percent_of_max_entropy": "#D55E00",
    "n_total": "#000000",
    "n_variants": "#F0E442",
    "variant_fraction": "#0072B2",
    "variant_fraction_percent": "green",
    "insertions": "#E69F00",
    "deletions": "black",
    "indel_fraction": "blue",
    "indel_substitution_ratio": "purple",
    "max_variant_base": "green",
}

# Define the tooltips for each option
column_tooltips = {
    "entropy": "Shannon entropy, measuring sequence diversity. Higher values indicate a more even distribution of bases.",
    "effective_entropy": "Entropy of non-reference (i.e., variant) reads. This value more accurately reflects the diversity of the variant library: it excludes the reference base counts, which are always expected to be the majority at any given position.",
    "percent_of_max_entropy": "Fraction of maximum possible entropy at position. Values significantly below 1 indicate that the sequence diversity is lower than expected.",
    "n_total": "Total number of reads covering each position.",
    "n_variants": "Number of variant reads at each position.",
    "variant_fraction": "Fraction of variant reads at each position.",
    "variant_fraction_percent": "Fraction of expected variant reads at each position. Values below than 1 may suggest that the library contains a large amount of non-mutated sequences.",
    "insertions": "Count of insertions at each position.",
    "deletions": "Count of deletions at each position.",
    "indel_fraction": "Fraction of reads with indels at each position.",
    "indel_substitution_ratio": "Ratio of indels to substitutions at each position.",
    "max_variant_base": "Counts of the most common non-reference base.",
}

# Columns to show in tabular form
tabular_cols = [
    "pos",
    "is_selected",
    "ref",
    "aligned_ref",
    "A",
    "C",
    "G",
    "T",
    "reads_all",
    "n_variants",
    "variant_fraction",
    "entropy",
    "effective_entropy",
    "percent_of_max_entropy",
    "insertions",
    "deletions",
    "indel_fraction",
    "indel_substitution_ratio",
    "alignment_mismatch",
    "max_variant_base",
]
