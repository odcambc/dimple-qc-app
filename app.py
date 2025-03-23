import pandas as pd

from shiny import reactive
from shiny.express import input, render, ui
from shiny.types import FileInfo

from shinywidgets import render_plotly

# from Bio import SeqIO

from faicons import icon_svg

from plotly_plots import base_position_vs_value_plot_plotly

from plots import violin_plot
from process_data import (
    process_per_base_file,
    update_per_base_df,
    update_mean_values_per_base,
)
import process_reference
from shared import (
    expected_columns,
    column_colors_dict,
    column_names_dict,
    column_tooltips,
    tabular_cols,
)

from process_reference import (
    process_reference_fasta,
    process_reference_genbank,
)

# from input_checkbox_group_tooltips import input_checkbox_group_tooltips

ui.page_opts(title="DIMPLE quick QC", fillable=True)
ui.include_css("styles.css")


# Store the max position (i.e., the length of the sequencing data)
sequence_length = reactive.value(100)

# Store the last selected series for violin plots
last_selected_series = reactive.value("entropy")


# Create a checkbox with a tooltip at the end
def checkbox_with_tooltip(key, names, tooltips):
    return ui.tags.label(
        names[key],
        ui.tags.span(" ", style="margin-right: 5px;"),
        ui.tags.span(
            icon_svg("circle-info", title=tooltips[key]),
            ui.tags.span(tooltips[key], class_="tooltip-text"),
            class_="tooltip-container",
        ),
        style="cursor: pointer;",
    )


# Sidebar layout
with ui.sidebar(title="Settings"):
    "Selected range:"
    with ui.layout_columns():

        @render.ui
        def min_pos_input():
            return ui.input_numeric(
                "min_pos",
                "Lower",
                min=0,
                value=0,
                update_on="blur",
            )

        @render.ui
        def max_pos_input():
            return ui.input_numeric(
                "max_pos",
                "Upper",
                max=sequence_length.get(),
                value=100,
                update_on="blur",
            )

    # Checkbox for reverse complementing per-base file
    ui.input_switch("reverse_complement", "Reverse complement")

    # Render the checkbox group with tooltips
    ui.input_checkbox_group(
        "data_series",
        "Display",
        {
            key: checkbox_with_tooltip(key, column_names_dict, column_tooltips)
            for key in column_tooltips
        },
    )
    # Input files
    ui.input_file(
        "per_base_file",
        "Upload sequencing data",
        accept=[".csv", ".tsv"],
        multiple=False,
    )
    ui.input_file(
        "reference_file",
        "Upload reference sequence",
        accept=[".fa:", ".fasta", ".gb", ".genbank", ".gbk"],
        multiple=False,
    )

    # Add a selectize input for identified features in parsed genbank file.
    # Only show if a genbank file is uploaded.
    @render.ui
    def feature_select():
        feature_names = ["None"]
        features = None

        if parsed_reference():
            if parsed_reference()["features"]:
                features = parsed_reference()["features"]
                feature_names = [
                    feature.qualifiers["label"][0]
                    for feature in features
                    if "label" in feature.qualifiers
                ]

        return ui.input_selectize(
            "selected_features",
            "Select feature",
            choices=feature_names,
            multiple=True,
        )


# Main content area
with ui.layout_columns(columns=2, col_widths=[9, 3]):
    with ui.card():
        with ui.navset_card_pill():
            with ui.nav_panel("Plots"):

                @render_plotly
                def plotly_position_plot():
                    pos_plot = base_position_vs_value_plot_plotly(
                        processed_per_base_file(),
                        input.data_series(),
                        [0, sequence_length.get()],
                        input.min_pos(),
                        input.max_pos(),
                        last_selected_series.get(),
                        input.show_means(),
                    )

                    return pos_plot

            with ui.nav_panel("Tabular data"):

                @render.data_frame
                def sequencing_data():
                    if processed_per_base_file().empty:
                        return pd.DataFrame()

                    return render.DataGrid(
                        pd.DataFrame(processed_per_base_file()[tabular_cols]),
                        filters=False,
                    )

        # Bottom 2D plots
        with ui.card(height=400):

            @render.plot
            @reactive.event(
                input.min_pos,
                input.max_pos,
                last_selected_series,
                input.data_series,
            )
            def render_value_violins():
                distribution_plots = violin_plot(
                    processed_per_base_file(),
                    last_selected_series.get(),
                    column_names_dict,
                )

                return distribution_plots

    # Top summary fields
    with ui.card():
        ui.input_switch("show_means", "Show means")

        with ui.value_box():
            "Average reads per base"

            @render.text
            def count():
                if mean_values_per_base().empty:
                    return "0"

                if "reads_all" not in mean_values_per_base().columns:
                    return "1"

                if True not in mean_values_per_base().index:
                    return str(mean_values_per_base().index)
                    return "2"

                selected_mean_reads = mean_values_per_base().at[
                    True, ("reads_all", "mean")
                ]

                if False in mean_values_per_base().index:
                    unselected_mean_reads = mean_values_per_base().at[
                        False, ("reads_all", "mean")
                    ]
                else:
                    unselected_mean_reads = selected_mean_reads

                return f"{int(unselected_mean_reads)} ({int(selected_mean_reads)} in selected)"

        with ui.value_box():

            @render.text
            def last_selected_series_text():
                return f"Average {column_names_dict[last_selected_series.get()]}"

            @render.text
            def avg_last_selected():
                selected_series = last_selected_series.get()

                if mean_values_per_base().empty:
                    return "0"

                if True not in mean_values_per_base().index:
                    return "0"

                if selected_series not in mean_values_per_base().columns:
                    return "0"

                selected_mean_series = mean_values_per_base().at[
                    True, (selected_series, "mean")
                ]

                if False in mean_values_per_base().index:
                    unselected_mean_series = mean_values_per_base().at[
                        False, (selected_series, "mean")
                    ]
                else:
                    unselected_mean_series = selected_mean_series

                return f"{unselected_mean_series:.2f} ({selected_mean_series:.2f} in selected)"

        with ui.value_box():
            "Average effective entropy"

            @render.text
            def avg_effective_entropy():
                if mean_values_per_base().empty:
                    return "0"

                if "effective_entropy" not in mean_values_per_base().columns:
                    return "0"

                if True not in mean_values_per_base().index:
                    return "0"

                selected_mean_entropy = mean_values_per_base().at[
                    True, ("effective_entropy", "mean")
                ]

                if False in mean_values_per_base().index:
                    unselected_mean_entropy = mean_values_per_base().at[
                        False, ("effective_entropy", "mean")
                    ]
                else:
                    unselected_mean_entropy = selected_mean_entropy

                return f"{unselected_mean_entropy:.2f} ({selected_mean_entropy:.2f} in selected)"


# Reactive calcs and effects


# Parse the input reference FASTA. Only handle single-sequence FASTA files.
@reactive.calc
def parsed_reference() -> dict[str, str | list | None] | None:
    file = input.reference_file()
    if file is None:
        return None

    # If it's fasta, parse as such
    if file[0]["name"].endswith((".fa", ".fasta")):
        return process_reference_fasta(file)
    # If it's genbank, parse as such
    if file[0]["name"].endswith((".gb", ".genbank", ".gbk")):
        return process_reference_genbank(file)

    # If it's neither, return None
    return None


@reactive.calc
def parsed_per_base_file():
    file: list[FileInfo] | None = input.per_base_file()
    if file is None:
        return pd.DataFrame()

    df = pd.read_csv(file[0]["datapath"], sep="\t")

    if not validate_per_base_file(df):
        return pd.DataFrame()

    # Plot entropy by default
    ui.update_checkbox_group("data_series", selected=["entropy"])

    return df


@reactive.calc
def processed_per_base_file():

    if parsed_per_base_file().empty:
        return pd.DataFrame()

    data = process_per_base_file(parsed_per_base_file(), input.reverse_complement())

    # Set the sequence length
    sequence_length.set(max(data["pos"]))

    return data


# TODO: this df should also include the means for the full series, not just selected vs. non-selected
@reactive.calc
def mean_values_per_base():
    if processed_per_base_file().empty:
        return pd.DataFrame()

    return update_mean_values_per_base(
        processed_per_base_file(), input.min_pos(), input.max_pos()
    )


# Reactive effects


# TODO: this is being performed _after_ mean_values_per_base is calculated, which is not ideal
@reactive.effect
def update_data_selected_range():
    # Calculate data that is dependent on the range selection

    if processed_per_base_file().empty:
        return

    update_per_base_df(
        processed_per_base_file(), input.min_pos(), input.max_pos(), parsed_reference()
    )


# Update alignment when a reference sequence is uploaded
@reactive.effect
@reactive.event(input.reference_file)
def update_alignment():
    return


# Dynamically update the max_pos input field when sequence_length updates
@reactive.effect
def update_max_pos():
    if not processed_per_base_file().empty:
        ui.update_numeric("max_pos", value=sequence_length.get())


# Update the last selected series from the checkbox group
@reactive.effect
def update_last_selected_series():
    if input.data_series():
        last_selected_series.set(input.data_series()[-1])


# Validate range inputs
@reactive.effect
def validate_range():
    # Minimum should be strictly less than maximum
    if input.min_pos() >= input.max_pos():
        ui.update_numeric("min_pos", value=input.max_pos() - 1)
    # Minimum should be greater than or equal to 0
    if input.min_pos() < 0:
        ui.update_numeric("min_pos", value=0)
    # Maximum should be less than the sequence length
    if input.max_pos() > sequence_length.get():
        ui.update_numeric("max_pos", value=sequence_length.get())


def validate_per_base_file(per_base_file):
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
