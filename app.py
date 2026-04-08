import pandas as pd

from shiny import reactive
from shiny.express import input, render, ui
from shiny.types import FileInfo

from shinywidgets import render_plotly

from faicons import icon_svg

from evaluate_data import test_per_base_file

from per_base_io import read_per_base_table

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




@reactive.calc
def last_selected_series() -> str:
    """Derive the last-checked series from the checkbox group.

    Using a calc (not a reactive.value set by an effect) means this resolves
    in the same flush as input.data_series changes, so downstream renders only
    fire once per user interaction.
    """
    series = input.data_series()
    return series[-1] if series else "entropy"


def format_mean_metric(means: pd.DataFrame, col: str, as_int: bool = False) -> str:
    """Format a metric from the means DataFrame as 'full (selected in selected)'."""
    if (
        means.empty
        or col not in means.columns
        or "selected" not in means.index
        or pd.isna(means.at["selected", col])
    ):
        return "0"

    selected = means.at["selected", col]
    if "unselected" in means.index and not pd.isna(means.at["unselected", col]):
        unselected = means.at["unselected", col]
    else:
        unselected = selected

    if as_int:
        return f"{int(unselected)} ({int(selected)} in selected)"
    return f"{unselected:.2f} ({selected:.2f} in selected)"

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
        # Static inputs: wrapping Upper in @render.ui with dynamic max=sequence_length
        # recreates the widget whenever length changes and resets value=100, which
        # fights ui.update_numeric and causes a rapid reset loop.
        ui.input_numeric(
            "min_pos",
            "Lower",
            min=0,
            value=0,
            update_on="blur",
        )
        ui.input_numeric(
            "max_pos",
            "Upper",
            min=0,
            max=10**9,
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
        accept=[".fa", ".fasta", ".gb", ".genbank", ".gbk"],
        multiple=False,
    )

    # Add a selectize input for identified features in parsed genbank file.
    # Only show if a genbank file is uploaded.
    @render.ui
    def feature_select():
        if parsed_reference() is None:
            return None
        if not parsed_reference()["features"]:
            return None

        features = parsed_reference()["features"]
        if not isinstance(features, dict):
            return None

        # Group features by type for the selectize optgroups
        grouped: dict[str, dict[str, str]] = {}
        for key, feat in features.items():
            group = feat.type if hasattr(feat, "type") else "Other"
            loc_start = int(feat.location.start)
            loc_end = int(feat.location.end)
            display = f"{key} ({loc_start}–{loc_end})"
            grouped.setdefault(group, {})[key] = display

        return ui.input_selectize(
            "selected_features",
            "Select feature",
            choices=grouped,
            multiple=True,
        )


# Main content area
with ui.layout_columns(columns=2, col_widths=[9, 3]):
    with ui.card():
        with ui.navset_card_pill():
            with ui.nav_panel("Plots"):

                @render_plotly
                def plotly_position_plot():
                    data = processed_per_base_file()
                    if data.empty:
                        return base_position_vs_value_plot_plotly(
                            data,
                            pd.DataFrame(),
                            [],
                            [0, 0],
                            0,
                            0,
                            last_selected_series(),
                            input.show_means(),
                        )
                    pos_plot = base_position_vs_value_plot_plotly(
                        data,
                        mean_values_per_base(),
                        input.data_series(),
                        [0, int(data["pos"].max())],
                        input.min_pos(),
                        input.max_pos(),
                        last_selected_series(),
                        input.show_means(),
                        feature_regions_for_plot(),
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
                last_selected_series,
            )
            def render_value_violins():
                distribution_plots = distribution_violin_plot_plotly(
                    processed_per_base_file(),
                    last_selected_series(),
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
                return format_mean_metric(
                    mean_values_per_base(), "reads_all_mean", as_int=True
                )

        with ui.value_box():

            @render.text
            def last_selected_series_text():
                return (
                    f"Average {column_names_dict[last_selected_series()]}"
                )

            @render.text
            def avg_last_selected():
                col = f"{last_selected_series()}_mean"
                if col not in mean_values_per_base().columns:
                    return "N/A"
                return format_mean_metric(mean_values_per_base(), col)

        with ui.value_box():
            "Average effective entropy"

            @render.text
            def avg_effective_entropy():
                return format_mean_metric(
                    mean_values_per_base(), "effective_entropy_mean"
                )


# Reactive calcs
# Parse the input reference file. Handles fasta or genbank by parsing extension.
@reactive.calc
def parsed_reference() -> dict[str, dict | None] | None:
    """Parse reference file (FASTA or GenBank) and return dict with sequence string and
    feature objects (if present). If no feature objects are present, the "feature" value
    is None. Returns None if no file is uploaded or if the file is not a valid fasta or
    genbank file."""

    file = input.reference_file()
    if file is None:
        return None

    # If it's fasta, parse as such
    if file[0]["name"].endswith((".fa", ".fasta")):
        result = process_reference_fasta(file)
    # If it's genbank, parse as such
    elif file[0]["name"].endswith((".gb", ".genbank", ".gbk")):
        result = process_reference_genbank(file)
    else:
        return None

    return result


@reactive.calc
def parsed_per_base_file():
    """Parse input per-base sequencing file."""
    file: list[FileInfo] | None = input.per_base_file()
    if file is None:
        return pd.DataFrame()

    try:
        df = read_per_base_table(file[0]["datapath"])
    except Exception:
        ui.notification_show("Could not parse the uploaded file. Check the format.")
        return pd.DataFrame()

    if df.empty:
        ui.notification_show("Could not parse the uploaded file. Check the format.")
        return pd.DataFrame()

    if not validate_per_base_file(df):
        return pd.DataFrame()

    return df


@reactive.calc
def sequence_length():
    """Derive sequence length from parsed file, independent of range inputs."""
    parsed = parsed_per_base_file()
    if parsed.empty:
        return 100
    return int(parsed["pos"].max())


@reactive.effect
@reactive.event(input.per_base_file)
def initialize_series_selection():
    """Initialize displayed metric selection once per upload."""
    if input.per_base_file() is None:
        return
    ui.update_checkbox_group("data_series", selected=["entropy"])


@reactive.calc
def base_processed_data():
    """Run expensive per-base processing. Independent of range inputs."""
    parsed = parsed_per_base_file()
    if parsed.empty:
        return pd.DataFrame()
    return process_per_base_file(parsed, input.reverse_complement())


@reactive.calc
def processed_per_base_file():
    """Apply range and feature selection to base processed data (cheap)."""
    data = base_processed_data()
    if data.empty:
        return pd.DataFrame()

    data = update_per_base_df(data, [(input.min_pos(), input.max_pos())])

    ref = parsed_reference()
    if ref and ref["features"] and input.selected_features():
        selected_range = []
        for feature in input.selected_features():
            if feature in ref["features"]:
                selected_range.append(
                    (
                        int(ref["features"][feature].location.start),
                        int(ref["features"][feature].location.end),
                    )
                )
        if selected_range:
            data = update_per_base_df(data, selected_range)

    return data


@reactive.calc
def feature_regions_for_plot() -> list[dict]:
    """Build a list of feature region dicts for plot highlighting.

    Selected features get a vivid color; unselected features get a muted gray.
    """
    ref = parsed_reference()
    if not ref or not ref["features"]:
        return []

    try:
        selected = set(input.selected_features() or [])
    except Exception:
        selected = set()

    selected_colors = {
        "CDS": "rgba(100, 149, 237, 0.35)",
        "gene": "rgba(100, 149, 237, 0.25)",
        "misc_feature": "rgba(255, 165, 0, 0.35)",
        "promoter": "rgba(0, 180, 0, 0.35)",
    }
    selected_default = "rgba(120, 120, 255, 0.30)"
    unselected_color = "rgba(180, 180, 180, 0.12)"

    regions = []
    for key, feat in ref["features"].items():
        feat_type = feat.type if hasattr(feat, "type") else "Other"
        if key in selected:
            color = selected_colors.get(feat_type, selected_default)
        else:
            color = unselected_color
        regions.append({
            "start": int(feat.location.start),
            "end": int(feat.location.end),
            "label": key,
            "color": color,
        })
    return regions


# TODO: this df should also include the means for the full series, not just selected vs. non-selected
@reactive.calc
def mean_values_per_base():
    data = processed_per_base_file()
    if data.empty:
        return pd.DataFrame()
    return update_mean_values_per_base(data)


@reactive.calc
def test_results():
    data = processed_per_base_file()
    if data.empty:
        return pd.DataFrame()
    selected_means = mean_values_per_base()
    full_means = process_full_mean_values(data)
    return test_per_base_file(data, selected_means, full_means)


# Reactive effects


@reactive.effect
@reactive.event(sequence_length)
def update_max_pos():
    """Update max position UI when sequence length changes."""
    if input.per_base_file() is not None:
        seq_len = sequence_length()
        ui.update_numeric("max_pos", value=seq_len, max=seq_len)



@reactive.effect
@reactive.event(input.min_pos, input.max_pos)
def validate_range():
    min_val = input.min_pos()
    max_val = input.max_pos()
    seq_len = sequence_length()

    new_min = min_val
    new_max = max_val

    if min_val < 0:
        new_min = 0
    if max_val > seq_len:
        new_max = seq_len
    if new_min >= new_max:
        new_min = max(new_max - 1, 0)

    if new_min != min_val:
        ui.update_numeric("min_pos", value=new_min)
    if new_max != max_val:
        ui.update_numeric("max_pos", value=new_max, max=seq_len)
