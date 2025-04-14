import pandas as pd

from shiny import reactive
from shiny.express import input, render, ui
from shiny.types import FileInfo

from shinywidgets import render_plotly

from faicons import icon_svg

from evaluate_data import test_per_base_file

from plotly_plots import (
    base_position_vs_value_plot_plotly,
    distribution_violin_plot_plotly,
)

from process_data import (
    process_full_mean_values,
    process_per_base_file,
    update_per_base_df,
    update_mean_values_per_base,
)
from shared import (
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


# Class to store the app state and keep track of reactive values
class AppState:
    def __init__(self):
        self.sequence_length = reactive.value(100)
        self.last_selected_series = reactive.value("entropy")
        self.processed_data = reactive.value(pd.DataFrame())
        self.selected_mean_values = reactive.value(pd.DataFrame())
        self.full_mean_values = reactive.value(pd.DataFrame())
        self.selected_features = reactive.value([])
        self.reference_data = reactive.value(None)
        self.summary_metrics = reactive.value(pd.DataFrame())


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

        if parsed_reference() is not None:
            if parsed_reference()["features"]:
                features = parsed_reference()["features"]
                print(f"Type of features: {type(features)}")
                feature_names = (
                    list(features.keys()) if isinstance(features, dict) else []
                )
            else:
                return None
        else:
            return None

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

            with ui.nav_panel("Statistical tests"):

                @render.data_frame
                def test_results_table():
                    if test_results().empty:
                        return pd.DataFrame()

                    return render.DataGrid(test_results(), filters=False)

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
                if means.empty or pd.isna(means.at["unselected", "reads_all_mean"]):
                    return "0"

                selected_mean_reads = means.at["selected", "reads_all_mean"]
                if "unselected" in means.index and not pd.isna(
                    means.at["unselected", "reads_all_mean"]
                ):
                    unselected_mean_reads = means.at["unselected", "reads_all_mean"]
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
                if means.empty or pd.isna(
                    means.at["selected", f"{selected_series}_mean"]
                ):
                    return "0"

                selected_mean_reads = means.at["selected", f"{selected_series}_mean"]
                if "unselected" in means.index and not pd.isna(
                    means.at["unselected", f"{selected_series}_mean"]
                ):
                    unselected_mean_reads = means.at[
                        "unselected", f"{selected_series}_mean"
                    ]

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
                    means.at["selected", "effective_entropy_mean"]
                ):
                    return "0"

                selected_mean_reads = means.at["selected", "effective_entropy_mean"]
                if "unselected" in means.index and not pd.isna(
                    means.at["unselected", "effective_entropy_mean"]
                ):
                    unselected_mean_reads = means.at[
                        "unselected", "effective_entropy_mean"
                    ]
                else:
                    unselected_mean_reads = selected_mean_reads

                return f"{unselected_mean_reads:.2f} ({selected_mean_reads:.2f} in selected)"


# Reactive calcs
# Parse the input reference FASTA. Only handle single-sequence FASTA files.
@reactive.calc
def parsed_reference() -> dict[str, dict | None] | None:
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
    app_state.selected_features.set(input.selected_features())

    app_state.full_mean_values.set(process_full_mean_values(data))

    return data


# TODO: this df should also include the means for the full series, not just selected vs. non-selected
@reactive.calc
def mean_values_per_base():
    if processed_per_base_file().empty:
        return pd.DataFrame()
    app_state.selected_mean_values.set(
        update_mean_values_per_base(
            app_state.processed_data.get(), input.min_pos(), input.max_pos()
        )
    )

    app_state.summary_metrics.set(
        test_per_base_file(
            app_state.processed_data.get(),
            app_state.selected_mean_values.get(),
            app_state.full_mean_values.get(),
        )
    )

    return app_state.selected_mean_values()


@reactive.calc
def test_results():
    if app_state.processed_data.get().empty:
        return pd.DataFrame()

    return app_state.summary_metrics.get()


# Reactive effects


@reactive.effect
@reactive.event(input.min_pos, input.max_pos)
def update_data_selected_range():
    """Update data selection when range changes."""

    print("updating data selected range")
    print(f"Selected features: {input.selected_features()}")
    print(f"Selected features state: {app_state.selected_features.get()}")

    if app_state.processed_data.get().empty:
        return

    selected_range = [(input.min_pos(), input.max_pos())]

    update_per_base_df(
        app_state.processed_data.get(),
        selected_range,
    )


@reactive.effect
@reactive.event(input.selected_features, app_state.selected_features)
def update_data_selected_range_features():
    """Update selected range when features are selected."""
    print("updating data selected range features")

    if app_state.processed_data.get().empty:
        return
    app_state.selected_features.set(input.selected_features())
    print(f"Type of input.selected_features: {type(input.selected_features())}")
    print(f"Selected features input: {input.selected_features()}")

    print(f"Type of app_state.processed_data: {type(app_state.processed_data.get())}")
    print(
        f"Type of app_state.selected_features: {type(app_state.selected_features.get())}"
    )
    print(f"Type of app_state.reference_data: {type(app_state.reference_data.get())}")
    print(f"Keys of app_state.reference_data: {app_state.reference_data.get().keys()}")
    print(
        f'Keys of features in app_state.reference_data: {app_state.reference_data.get()["features"].keys()}'
    )
    try:
        print(
            f'Type of app_state.reference_data: {type(app_state.reference_data.get()["features"])}'
        )
    except TypeError:
        print("app_state.reference_data is None")
    print(f"Selected features: {app_state.selected_features.get()}")

    selected_range = []
    if not app_state.reference_data.get():
        print("No reference data")
        return
    if not app_state.selected_features.get():
        print("No selected features")
        return
    if not app_state.reference_data.get():
        print("No features in reference data")
        return

    # Add the feature ranges to selected range
    for feature in app_state.selected_features.get():
        print(f"Selected feature: {feature}")
        if feature in app_state.reference_data.get()["features"]:
            print(
                f'Feature start: {app_state.reference_data.get()["features"][feature].location.start}, feature end: {app_state.reference_data.get()["features"][feature].location.end}'
            )
            selected_range.append(
                (
                    int(
                        app_state.reference_data.get()["features"][
                            feature
                        ].location.start
                    ),
                    int(
                        app_state.reference_data.get()["features"][feature].location.end
                    ),
                )
            )

    print(selected_range)

    update_per_base_df(
        app_state.processed_data.get(),
        selected_range,
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
    print("updating range")
    print(f"Min pos:", input.min_pos())
    # Minimum should be strictly less than maximum
    if input.min_pos() >= input.max_pos():
        ui.update_numeric("min_pos", value=input.max_pos() - 1)
    # Minimum should be greater than or equal to 0
    if input.min_pos() < 0:
        ui.update_numeric("min_pos", value=0)
    # Maximum should be less than the sequence length
    if input.max_pos() > app_state.sequence_length.get():
        ui.update_numeric("max_pos", value=app_state.sequence_length.get())
