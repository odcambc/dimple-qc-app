from numpy import mean
import pandas as pd

from shiny import reactive
from shiny.express import input, render, ui
from shiny.types import FileInfo

from shinywidgets import render_plotly

from faicons import icon_svg

from plotly_plots import (
    base_position_vs_value_plot_plotly,
    distribution_violin_plot_plotly,
)

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

from validation import validate_per_base_file


class AppState:
    def __init__(self):
        self.sequence_length = reactive.value(100)
        self.last_selected_series = reactive.value("entropy")
        self.processed_data = reactive.value(pd.DataFrame())
        self.reference_data = reactive.value(None)


app_state = AppState()

# from input_checkbox_group_tooltips import input_checkbox_group_tooltips

ui.page_opts(title="DIMPLE quick QC", fillable=True)
ui.include_css("styles.css")


# Reactive calcs and effects
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
                max=app_state.sequence_length.get(),
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
                        mean_values_per_base(),
                        input.data_series(),
                        [0, app_state.sequence_length.get()],
                        input.min_pos(),
                        input.max_pos(),
                        app_state.last_selected_series.get(),
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

            @render_plotly
            @reactive.event(
                input.min_pos,
                input.max_pos,
                input.data_series,
                app_state.last_selected_series,
            )
            def render_value_violins():
                distribution_plots = distribution_violin_plot_plotly(
                    processed_per_base_file(),
                    app_state.last_selected_series.get(),
                    column_names_dict,
                    column_colors_dict,
                )

                return distribution_plots

    # Top summary fields
    with ui.card():
        ui.input_switch("show_means", "Show means")

        with ui.value_box():
            "Average reads per base"

            @render.text
            def count():
                means = mean_values_per_base()
                # If the DataFrame is empty or contains NaNs, return a default value
                if means.empty or pd.isna(means.at[True, ("reads_all", "mean")]):
                    return "0"

                selected_mean_reads = means.at[True, ("reads_all", "mean")]
                if False in means.index and not pd.isna(
                    means.at[False, ("reads_all", "mean")]
                ):
                    unselected_mean_reads = means.at[False, ("reads_all", "mean")]
                else:
                    unselected_mean_reads = selected_mean_reads

                return f"{int(unselected_mean_reads)} ({int(selected_mean_reads)} in selected)"

        with ui.value_box():

            @render.text
            def last_selected_series_text():
                return (
                    f"Average {column_names_dict[app_state.last_selected_series.get()]}"
                )

            @render.text
            def avg_last_selected():
                selected_series = app_state.last_selected_series.get()

                means = mean_values_per_base()
                # If the DataFrame is empty or contains NaNs, return a default value
                if means.empty or pd.isna(means.at[True, (selected_series, "mean")]):
                    return "0"

                selected_mean_reads = means.at[True, (selected_series, "mean")]
                if False in means.index and not pd.isna(
                    means.at[False, (selected_series, "mean")]
                ):
                    unselected_mean_reads = means.at[False, (selected_series, "mean")]
                else:
                    unselected_mean_reads = selected_mean_reads

                return f"{unselected_mean_reads:.2f} ({selected_mean_reads:.2f} in selected)"

        with ui.value_box():
            "Average effective entropy"

            @render.text
            def avg_effective_entropy():
                means = mean_values_per_base()
                # If the DataFrame is empty or contains NaNs, return a default value
                if means.empty or pd.isna(
                    means.at[True, ("effective_entropy", "mean")]
                ):
                    return "0"

                selected_mean_reads = means.at[True, ("effective_entropy", "mean")]
                if False in means.index and not pd.isna(
                    means.at[False, ("effective_entropy", "mean")]
                ):
                    unselected_mean_reads = means.at[
                        False, ("effective_entropy", "mean")
                    ]
                else:
                    unselected_mean_reads = selected_mean_reads

                return f"{unselected_mean_reads:.2f} ({selected_mean_reads:.2f} in selected)"


# Reactive calcs
# Parse the input reference FASTA. Only handle single-sequence FASTA files.
@reactive.calc
def parsed_reference() -> dict[str, str | list | None] | None:
    """Parse reference file (FASTA or GenBank) and return dict with sequence string and
    feature objects (if present)."""
    file = input.reference_file()
    if file is None:
        return None

    # If it's fasta, parse as such
    if file[0]["name"].endswith((".fa", ".fasta")):
        result = process_reference_fasta(file)
    # If it's genbank, parse as such
    if file[0]["name"].endswith((".gb", ".genbank", ".gbk")):
        result = process_reference_genbank(file)
    else:
        return None

    app_state.reference_data.set(result)
    return result


@reactive.calc
def parsed_per_base_file():
    """Parse input per-base sequencing file."""
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

    # Update the shared state
    app_state.sequence_length.set(max(data["pos"]))
    app_state.processed_data.set(data)

    return data


# TODO: this df should also include the means for the full series, not just selected vs. non-selected
@reactive.calc
def mean_values_per_base():
    if processed_per_base_file().empty:
        return pd.DataFrame()

    return update_mean_values_per_base(
        app_state.processed_data.get(), input.min_pos(), input.max_pos()
    )


# Reactive effects


# TODO: this is being performed _after_ mean_values_per_base is calculated, which is not ideal
@reactive.effect
@reactive.event(input.min_pos, input.max_pos)
def update_data_selected_range():
    """Update data selection when range changes."""

    if app_state.processed_data.get().empty:
        return

    update_per_base_df(
        app_state.processed_data.get(),
        input.min_pos(),
        input.max_pos(),
        app_state.reference_data.get(),
    )


# Update alignment when a reference sequence is uploaded
@reactive.effect
@reactive.event(input.reference_file)
def update_alignment():
    return


@reactive.effect
@reactive.event(app_state.sequence_length)
def update_max_pos():
    """Update max position UI when sequence length changes."""
    if not app_state.processed_data.get().empty:
        ui.update_numeric("max_pos", value=app_state.sequence_length.get())


@reactive.effect
@reactive.event(input.data_series)
def update_last_selected_series():
    """Update last selected series when checkbox selection changes."""
    if input.data_series():
        app_state.last_selected_series.set(input.data_series()[-1])


# Validate range inputs
@reactive.effect
@reactive.event(input.min_pos, input.max_pos)
def validate_range():
    # Minimum should be strictly less than maximum
    if input.min_pos() >= input.max_pos():
        ui.update_numeric("min_pos", value=input.max_pos() - 1)
    # Minimum should be greater than or equal to 0
    if input.min_pos() < 0:
        ui.update_numeric("min_pos", value=0)
    # Maximum should be less than the sequence length
    if input.max_pos() > app_state.sequence_length.get():
        ui.update_numeric("max_pos", value=app_state.sequence_length.get())
