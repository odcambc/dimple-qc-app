import seaborn as sns

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

from shiny import reactive
from shiny.express import input, render, ui
from shiny.types import FileInfo

from Bio import SeqIO
from Bio import pairwise2

ui.page_opts(title="DIMPLE quick QC", fillable=True)
ui.include_css("./styles.css")

column_names_dict = {
    "entropy": "Entropy",
    "entropy_without_max_value": "Effective entropy",
    "n_total": "Total reads",
    "n_variants": "Variant reads",
    "variant_fraction": "Variant fraction",
    "percent_of_max_entropy": "Entropy % max",
    "percent_of_max_entropy_estimated": "Est. entropy % max",
    "insertions": "Insertion count",
    "deletions": "Deletion count",
    "indel_fraction": "Indel fraction",
    "indel_missense_ratio": "Indel to missense ratio",
    "alignment_mismatch": "Alignment mismatches",
    "max_variant_base": "Max counts non-ref base",
}

column_colors_dict = {
    "entropy": "#009E73",
    "entropy_without_max_value": "#56B4E9",
    "n_total": "#000000",
    "n_variants": "#F0E442",
    "variant_fraction": "#0072B2",
    "percent_of_max_entropy": "#D55E00",
    "percent_of_max_entropy_estimated": "#CC79A7",
    "insertions": "#E69F00",
    "deletions": "black",
    "indel_fraction": "blue",
    "indel_missense_ratio": "purple",
    "alignment_mismatch": "red",
    "max_variant_base": "green",
}

# Store the selected upper and lower range for plotting
plot_range_low = reactive.value(0)
plot_range_high = reactive.value(100)

# Store the max position (i.e., the length of the sequencing data)
sequence_length = reactive.value(100)

# Store the selected upper and lower selection range
selected_range_low = reactive.value(0)
selected_range_high = reactive.value(100)

# Store the last selected series for violin plots
last_selected_series = reactive.value("entropy")

# Sidebar layout
with ui.sidebar(title="Settings"):

    # Range selection
    @render.ui
    def range_slider():
        return ui.input_slider(
            "selected_range",
            "Select range",
            min=0,
            max=sequence_length.get(),
            value=[selected_range_low.get(), selected_range_high.get()],
        )

    with ui.layout_columns():

        @render.ui
        def min_pos_input():
            return ui.input_numeric(
                "min_pos", "Lower", min=0, value=selected_range_low.get()
            )

        @render.ui
        def max_pos_input():
            return ui.input_numeric(
                "max_pos",
                "Upper",
                max=sequence_length.get(),
                value=selected_range_high.get(),
            )

    # Data series selection
    ui.input_checkbox_group(
        "data_series",
        "Display",
        column_names_dict,
    )

    # Input files
    ui.input_file(
        "per_base_file",
        "Upload sequencing data",
        accept=[".csv", ".tsv"],
        multiple=False,
    )
    ui.input_file(
        "reference_fasta",
        "Upload reference sequence",
        accept=[".fa:", ".fasta"],
        multiple=False,
    )


# Main content area
with ui.layout_columns(columns=2, col_widths=[9, 3]):
    with ui.card():
        with ui.navset_card_pill():
            with ui.nav_panel("Plots"):

                @render.plot
                def render_scatterplot():
                    if processed_per_base_file().empty:
                        return plt.figure()
                    pos_plot = sns.scatterplot()
                    pos_plot.axvline(
                        selected_range_low.get(), color="black", linestyle="--"
                    )
                    pos_plot.axvline(
                        selected_range_high.get(), color="black", linestyle="--"
                    )
                    for series in input.data_series():
                        sns.scatterplot(
                            data=processed_per_base_file(),
                            x="pos",
                            y=series,
                            color=column_colors_dict[series],
                            edgecolors="black",
                        )
                    pos_plot.set_xlim(plot_range_low.get(), plot_range_high.get())
                    pos_plot.set_xlabel("Position")
                    pos_plot.set_ylabel(column_names_dict[last_selected_series.get()])

                    # If show_means is selected, draw horizontal lines at means of last
                    # selected series for the full range and selected range

                    if input.show_means():
                        if input.data_series():
                            full_mean = processed_per_base_file()[
                                last_selected_series.get()
                            ].mean(skipna=True)
                            selected_mean = (
                                processed_per_base_file()
                                .loc[
                                    processed_per_base_file()["pos"].between(
                                        selected_range_low.get(),
                                        selected_range_high.get(),
                                    ),
                                    last_selected_series.get(),
                                ]
                                .mean(skipna=True)
                            )

                            pos_plot.axhline(
                                float(full_mean),
                                color=column_colors_dict[last_selected_series.get()],
                                linestyle="solid",
                            )
                            pos_plot.axhline(
                                float(selected_mean),
                                color=column_colors_dict[last_selected_series.get()],
                                linestyle="dashed",
                            )
                            # Add legend with correct line styles displayed
                            from matplotlib.lines import Line2D

                            legend_handles = [
                                Line2D(
                                    [0],
                                    [0],
                                    color=column_colors_dict[
                                        last_selected_series.get()
                                    ],
                                    linestyle="solid",
                                    label="Full range",
                                ),
                                Line2D(
                                    [0],
                                    [0],
                                    color=column_colors_dict[
                                        last_selected_series.get()
                                    ],
                                    linestyle="dashed",
                                    label="Selected range",
                                ),
                            ]
                            pos_plot.legend(
                                handles=legend_handles,
                                loc="upper right",
                            )

                        return pos_plot

                with ui.layout_columns():
                    ui.input_action_button("zoom", "Zoom")
                    ui.input_action_button("reset_zoom", "Reset")

            with ui.nav_panel("Tabular data"):

                @render.data_frame
                def sequencing_data():
                    if processed_per_base_file().empty:
                        return pd.DataFrame()
                    cols = [
                        "pos",
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
                        "entropy_without_max_value",
                        "estimated_entropy",
                        "percent_of_max_entropy",
                        "percent_of_max_entropy_estimated",
                    ]

                    return render.DataGrid(
                        pd.DataFrame(processed_per_base_file()[cols]), filters=False
                    )

        # Bottom 2D plots
        with ui.card(height=400):

            ui.card_header("Distributions")

            @render.plot
            @reactive.event(
                selected_range_low, selected_range_high, last_selected_series
            )
            def render_value_violins():
                if processed_per_base_file().empty:
                    return plt.figure()
                if last_selected_series.get() not in processed_per_base_file():
                    return plt.figure()

                # Now refactor the df to allow for easy plotting.
                # We want to plot the distributions of the selected series
                # for the full range of the data vs the selected range.
                # We can use the "is_selected" column to split the data.

                plotting_df = processed_per_base_file().melt(
                    id_vars=["is_selected"],
                    value_vars=last_selected_series.get(),
                    var_name="series",
                    value_name="value",
                )

                violin_plot = sns.violinplot(
                    data=plotting_df,
                    x="value",
                    hue="is_selected",
                    split=True,
                    linewidth=1,
                    linecolor="black",
                    fill=True,
                )
                violin_plot.set_xlabel(column_names_dict[last_selected_series.get()])
                violin_plot.set_ylabel("Density")
                violin_plot.set_title(
                    f"Distribution of {column_names_dict[last_selected_series.get()]}"
                )
                violin_plot.legend(
                    title="Selected range",
                    loc="upper right",
                )

                return violin_plot

    # Top summary fields
    with ui.card():
        ui.input_switch("show_means", "Show means")

        with ui.value_box():
            "Average reads per base"

            @render.text
            def count():
                df = processed_per_base_file()
                if df.empty:
                    return "0"

                selected_range = range(
                    selected_range_low.get(), selected_range_high.get()
                )

                full_avg = df["reads_all"].mean(skipna=True)
                range_avg = df.loc[df["pos"].isin(selected_range), "reads_all"].mean(
                    skipna=True
                )

                return f"{int(full_avg)} ({int(range_avg)} in selected)"

        with ui.value_box():

            @render.text
            def last_selected_series_text():
                return f"Average {column_names_dict[last_selected_series.get()]}"

            @render.text
            def avg_last_selected():
                df = processed_per_base_file()
                if df.empty:
                    return "0"

                selected_range = range(
                    selected_range_low.get(), selected_range_high.get()
                )

                selected_series = last_selected_series.get()

                full_avg = df[selected_series].mean(skipna=True)
                range_avg = df.loc[df["pos"].isin(selected_range)][
                    selected_series
                ].mean(skipna=True)

                return f"{full_avg:.2f} ({range_avg:.2f} in selected)"

        with ui.value_box():
            "Average effective entropy"

            @render.text
            def avg_entropy_without_max_value():
                df = processed_per_base_file()
                if df.empty:
                    return "0"

                selected_range = range(
                    selected_range_low.get(), selected_range_high.get()
                )

                full_avg = df["entropy_without_max_value"].mean(skipna=True)
                range_avg = df.loc[df["pos"].isin(selected_range)][
                    "entropy_without_max_value"
                ].mean(skipna=True)

                return f"{full_avg:.2f} ({range_avg:.2f} in selected)"


# Reactive calcs and effects
@reactive.calc
def parsed_reference_fasta():
    file = input.reference_fasta()
    if file is None:
        return None
    try:
        fasta_record = list(SeqIO.parse(file[0]["datapath"], "fasta"))
        if len(fasta_record) != 1:
            return None  # Only handle single-sequence FASTA files
        return str(fasta_record[0].seq)
    except Exception:
        return None


@reactive.calc
def parsed_per_base_file():
    file: list[FileInfo] | None = input.per_base_file()
    if file is None:
        return pd.DataFrame()

    # Plot entropy by default
    ui.update_checkbox_group("data_series", selected=["entropy"])

    return pd.read_csv(
        file[0]["datapath"], sep="\t"
    )  # pyright: ignore[reportUnknownMemberType]


@reactive.calc
def processed_per_base_file():

    if parsed_per_base_file().empty:
        return pd.DataFrame()

    data = parsed_per_base_file()

    # Set the sequence length
    sequence_length.set(max(data["pos"]))

    # Update the selected range and plotting range
    selected_range_high.set(sequence_length.get())
    plot_range_high.set(sequence_length.get())

    selected_codon_range = range(0, sequence_length.get() // 3)

    subpool_codon_fraction = 3 / (len(selected_codon_range) + 1)

    data["codon_number"] = data["pos"] // 3

    data["n_variants"] = data[["A", "C", "G", "T"]].apply(n_variants, axis=1)
    data["n_indels"] = data[["insertions", "deletions"]].apply(n_observations, axis=1)
    data["n_total"] = data[["A", "C", "G", "T"]].apply(n_observations, axis=1)

    # Calculate ratios while avoiding division by zero

    data["variant_fraction"] = (data["n_variants"] / data["reads_all"]).replace(
        [np.inf, -np.inf], np.nan
    )

    data["indel_fraction"] = (data["n_indels"] / data["reads_all"]).replace(
        [np.inf, -np.inf], np.nan
    )

    data["indel_missense_ratio"] = (data["n_indels"] / data["n_variants"]).replace(
        [np.inf, -np.inf], np.nan
    )

    # Calculate the counts of the most common non-reference base

    data["max_variant_base"] = data[["ref", "A", "C", "G", "T"]].apply(
        max_non_ref_base, axis=1
    )

    data["expected_variant_codons"] = subpool_codon_fraction * data["n_total"]
    data["estimated_wt_n"] = (
        data["expected_variant_codons"] - data["n_variants"]
    ).clip(lower=0)

    data["entropy"] = data[["A", "C", "G", "T"]].apply(
        lambda x: stats.entropy(x), axis=1
    )
    data["entropy_without_max_value"] = data[["A", "C", "G", "T"]].apply(
        entropy_without_max_value, axis=1
    )

    data["estimated_entropy"] = data[["A", "C", "G", "T", "estimated_wt_n"]].apply(
        entropy_without_max_value, axis=1
    )
    data["entropy_diff"] = data["entropy"] - data["estimated_entropy"]
    data["percent_of_max_entropy_estimated"] = data["estimated_entropy"] / np.log(4)
    data["percent_of_max_entropy"] = data["entropy_without_max_value"] / np.log(3)

    # Create empty columns for alignment to start
    data["aligned_ref"] = ["-"] * len(data)
    data["alignment_mismatch"] = [0] * len(data)

    return data


# Helper functions for data processing
def entropy_without_max_value(x):
    row_list = [i for i in x if i != max(x)]

    # Set values below 2 to 0
    row_list = [0 if i < 3 else i for i in row_list]

    return stats.entropy(row_list)


def n_variants(x):
    return sum([i for i in x if i != max(x)])


def n_observations(x):
    return sum(x)


def max_non_ref_base(x):
    ref = x["ref"]
    non_ref_bases = [i for i in x[["A", "C", "G", "T"]] if i != x[ref]]
    if non_ref_bases:
        return max(non_ref_bases)
    return 0


# Reactive effects
@reactive.effect
def update_data_selected_range():
    # Calculate data that is dependent on the range selection

    if processed_per_base_file().empty:
        return

    selected_codon_range = range(
        selected_range_low.get() // 3, selected_range_high.get() // 3
    )
    subpool_codon_fraction = 3 / (len(selected_codon_range) + 1)

    parsed_per_base_file()["is_selected"] = parsed_per_base_file()["codon_number"].isin(
        selected_codon_range
    )

    parsed_per_base_file()["codon_fraction"] = parsed_per_base_file()[
        "is_selected"
    ].apply(lambda x: subpool_codon_fraction if x else 1 - subpool_codon_fraction)


# Update alignment when a reference sequence is uploaded
@reactive.effect
@reactive.event(input.reference_fasta)
def update_alignment():
    if not parsed_reference_fasta():
        return
    # Align reference sequence optimally
    ref_seq = parsed_reference_fasta()
    if ref_seq:
        df_ref_sequence = "".join(
            processed_per_base_file()["ref"]
        )  # Convert ref column to a string

        # Perform global alignment (Needleman-Wunsch)
        alignments = pairwise2.align.globalxx(
            ref_seq, df_ref_sequence, one_alignment_only=True
        )

        if alignments:
            best_alignment = alignments[0]  # Best alignment result
            aligned_seq = list(best_alignment.seqA)  # Extract aligned reference

            # Handle insertions/deletions (gap handling)
            adjusted_aligned_seq = []
            mismatch_count = 0

            for i, (ref_base, df_base) in enumerate(zip(aligned_seq, df_ref_sequence)):
                if ref_base == "-":  # Reference has an insertion
                    adjusted_aligned_seq.append("-")  # Keep the gap
                elif ref_base != df_base:
                    adjusted_aligned_seq.append(f"[{ref_base}]")  # Mark mismatch
                    mismatch_count += 1
                else:
                    adjusted_aligned_seq.append(ref_base)

            processed_per_base_file()["aligned_ref"] = adjusted_aligned_seq

            processed_per_base_file()["alignment_mismatch"] = [
                1 if ref_base != df_base or ref_base == "-" else 0
                for ref_base, df_base in zip(aligned_seq, df_ref_sequence)
            ]


# Dynamically update the max_pos input field when sequence_length updates
@reactive.effect
def update_max_pos():
    if not processed_per_base_file().empty:
        ui.update_numeric("max_pos", value=sequence_length.get())
        ui.update_slider("pos_range", max=sequence_length.get())


@reactive.effect
@reactive.event(input.selected_range)
def updated_selected_range():
    selected_range_low.set(input.selected_range()[0])
    selected_range_high.set(input.selected_range()[1])
    ui.update_numeric("min_pos", value=selected_range_low.get())
    ui.update_numeric("max_pos", value=selected_range_high.get())


@reactive.effect
@reactive.event(input.min_pos)
def updated_min_pos():
    selected_range_low.set(input.min_pos())
    ui.update_slider(
        "selected_range", value=(selected_range_low.get(), selected_range_high.get())
    )
    ui.update_numeric("min_pos", value=selected_range_low.get())


@reactive.effect
@reactive.event(input.max_pos)
def updated_max_pos():
    selected_range_high.set(input.max_pos())
    ui.update_slider(
        "selected_range", value=(selected_range_low.get(), selected_range_high.get())
    )
    ui.update_numeric("max_pos", value=selected_range_high.get())


# Update zoom when button is clicked
@reactive.effect
@reactive.event(input.zoom)
def zoom():
    # ui.update_numeric("min_pos", value=input.min_pos())
    # ui.update_numeric("max_pos", value=input.max_pos())
    plot_range_low.set(input.min_pos())
    plot_range_high.set(input.max_pos())


# Reset zoom when button is clicked
@reactive.effect
@reactive.event(input.reset_zoom)
def reset_zoom():
    # ui.update_slider("pos_range", value=[0, sequence_length])
    # ui.update_numeric("min_pos", value=0)
    plot_range_low.set(0)
    plot_range_high.set(sequence_length.get())


# Update the last selected series from the checkbox group
@reactive.effect
def update_last_selected_series():
    if input.data_series():
        last_selected_series.set(input.data_series()[-1])
