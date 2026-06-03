import pandas as pd
import plotly.graph_objects as go

from pathlib import Path

from shiny import reactive
from shiny.express import input, render, ui
from shiny.types import FileInfo
from shiny.ui import Theme

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
    align_ref_to_variants,
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

# Theme from the suite-wide brand (vendored _brand.yml; see dms-tools brand/).
ui.page_opts(
    title="DIMPLE quick QC",
    fillable=True,
    theme=Theme.from_brand(Path(__file__).parent / "_brand.yml"),
)
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


# These calcs must be defined before the UI block because they are referenced
# as arguments to @reactive.event decorators, which are evaluated at definition time.
@reactive.calc
def parsed_reference() -> dict[str, dict | None] | None:
    """Parse reference file (FASTA or GenBank) and return dict with sequence string and
    feature objects (if present). If no feature objects are present, the "feature" value
    is None. Returns None if no file is uploaded or if the file is not a valid fasta or
    genbank file."""
    file = input.reference_file()
    if file is None:
        return None
    if file[0]["name"].endswith((".fa", ".fasta")):
        result = process_reference_fasta(file)
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


@reactive.calc
def base_processed_data():
    """Run expensive per-base processing. Independent of range inputs."""
    parsed = parsed_per_base_file()
    if parsed.empty:
        return pd.DataFrame()
    data = process_per_base_file(parsed, input.reverse_complement(), input.origin_shift())
    ref = parsed_reference()
    if ref and ref.get("sequence"):
        data = align_ref_to_variants(data, ref["sequence"])
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
        color = selected_colors.get(feat_type, selected_default) if key in selected else unselected_color
        regions.append({
            "start": int(feat.location.start),
            "end": int(feat.location.end),
            "label": key,
            "color": color,
        })
    return regions


# Sidebar layout
with ui.sidebar(title="Settings"):
    ui.input_slider(
        "pos_range",
        "Selected range:",
        min=0,
        max=100,
        value=[0, 100],
        step=1,
    )

    # Checkbox for reverse complementing per-base file
    ui.input_switch("reverse_complement", "Reverse complement")

    # Numeric input for origin shift
    ui.input_numeric("origin_shift", "Origin shift (bp)", value=0, min=0, step=1)

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

    ui.hr()
    "Download:"

    @render.download(
        filename="per_base_data.csv", media_type="text/csv", label="Per-base data"
    )
    def download_per_base_csv():
        """Download per-base data as CSV."""
        df = processed_per_base_file()
        if df.empty:
            yield ""
            return
        yield df[tabular_cols].to_csv(index=False)

    @render.download(
        filename="test_results.csv", media_type="text/csv", label="Test results"
    )
    def download_test_csv():
        """Download statistical test results as CSV."""
        df = test_results()
        if df.empty:
            yield ""
            return
        yield df.to_csv()

    @render.download(
        filename="position_plot.png",
        media_type="image/png",
        label="Position plot (PNG)",
    )
    def download_plot_png():
        """Download position plot as PNG."""
        data = processed_per_base_file()
        fig = base_position_vs_value_plot_plotly(
            data,
            mean_values_per_base(),
            input.data_series(),
            [0, int(data["pos"].max()) if not data.empty else 0],
            input.pos_range()[0],
            input.pos_range()[1],
            last_selected_series(),
            input.show_means(),
            feature_regions_for_plot(),
        )
        yield fig.to_image(format="png")

    @render.download(
        filename="position_plot.html",
        media_type="text/html",
        label="Position plot (HTML)",
    )
    def download_plot_html():
        """Download position plot as interactive HTML."""
        data = processed_per_base_file()
        fig = base_position_vs_value_plot_plotly(
            data,
            mean_values_per_base(),
            input.data_series(),
            [0, int(data["pos"].max()) if not data.empty else 0],
            input.pos_range()[0],
            input.pos_range()[1],
            last_selected_series(),
            input.show_means(),
            feature_regions_for_plot(),
        )
        yield fig.to_html(include_plotlyjs="cdn").encode()


# Main content area
with ui.layout_columns(columns=2, col_widths=[9, 3]):
    with ui.card():
        with ui.navset_card_pill():
            with ui.nav_panel("Plots"):

                @render_plotly
                @reactive.event(
                    base_processed_data,
                    input.data_series,
                    last_selected_series,
                    feature_regions_for_plot,
                )
                def plotly_position_plot():
                    data = base_processed_data()
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
                        input.pos_range()[0],
                        input.pos_range()[1],
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
                input.pos_range,
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


@reactive.effect
@reactive.event(input.per_base_file)
def initialize_series_selection():
    """Initialize displayed metric selection once per upload."""
    if input.per_base_file() is None:
        return
    ui.update_checkbox_group("data_series", selected=["entropy"])


@reactive.calc
def processed_per_base_file():
    """Apply range and feature selection to base processed data (cheap)."""
    data = base_processed_data()
    if data.empty:
        return pd.DataFrame()

    data = update_per_base_df(data, [(input.pos_range()[0], input.pos_range()[1])])

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
def update_pos_range():
    """Update slider bounds and reset value when a new file is loaded."""
    if input.per_base_file() is not None:
        seq_len = sequence_length()
        ui.update_slider("pos_range", min=0, max=seq_len, value=[0, seq_len])


@reactive.effect
@reactive.event(input.pos_range, mean_values_per_base, input.show_means)
def update_position_plot_shapes():
    """Update vlines/hlines via widget delta — avoids re-rendering scatter traces."""
    w = plotly_position_plot.widget
    if w is None:
        return

    min_p, max_p = input.pos_range()

    # Build shapes on a temp figure using the same Plotly API as plotly_plots.py
    tmp = go.Figure()
    tmp.add_vline(x=min_p, line_width=1, line_dash="dash", line_color="black")
    tmp.add_vline(x=max_p, line_width=1, line_dash="dash", line_color="black")

    for region in feature_regions_for_plot():
        tmp.add_vrect(
            x0=region["start"],
            x1=region["end"],
            fillcolor=region.get("color", "rgba(100, 100, 255, 0.15)"),
            layer="below",
            line_width=0,
            annotation_text=region.get("label", ""),
            annotation_position="top left",
            annotation_font_size=10,
            annotation_font_color="gray",
        )

    if input.show_means():
        means = mean_values_per_base()
        if not means.empty:
            col = f"{last_selected_series()}_mean"
            if col in means.columns:
                sel_mean = means.at["selected", col]
                all_mean = (
                    means.at["unselected", col]
                    if "unselected" in means.index
                    and not pd.isna(means.at["unselected", col])
                    else sel_mean
                )
                if not pd.isna(sel_mean):
                    tmp.add_hline(
                        y=float(sel_mean),
                        line_width=1,
                        line_dash="dash",
                        line_color="black",
                        annotation_text=f"Selected mean: {float(sel_mean):.2f}",
                    )
                if not pd.isna(all_mean):
                    tmp.add_hline(
                        y=float(all_mean),
                        line_width=1,
                        line_dash="solid",
                        line_color="black",
                        annotation_text=f"Full mean: {float(all_mean):.2f}",
                    )

    w.update_layout(
        shapes=list(tmp.layout.shapes),
        annotations=list(tmp.layout.annotations),
    )


@reactive.effect
@reactive.event(input.origin_shift)
def validate_origin_shift():
    """Clamp origin shift to valid range [0, sequence_length - 1]."""
    seq_len = sequence_length()
    if seq_len <= 0:
        return
    if input.origin_shift() < 0:
        ui.update_numeric("origin_shift", value=0)
    elif input.origin_shift() >= seq_len:
        ui.update_numeric("origin_shift", value=seq_len - 1)
