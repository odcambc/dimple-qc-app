import time

import pandas as pd
import plotly.graph_objects as go

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
    example_per_base_tsv,
    example_reference_fasta,
    plottable_series,
    reverse_complement_sequence,
    tabular_cols,
)

from process_reference import (
    align_ref_to_variants,
    process_reference_fasta,
    process_reference_genbank,
)

from validation import validate_per_base_file




def reactive_debounce(value_fn, delay_secs):
    """Per-session debounce: a reactive.calc yielding value_fn()'s latest value,
    but only after value_fn has been stable for delay_secs.

    NOTE on the module-level-state convention: this creates reactive.value()s, which
    CLAUDE.md generally forbids at module scope. That rule guards against state shared
    across sessions, which happens for *imported* modules (import-cached once). app.py is
    different: Shiny Express re-executes this file's body per session
    (shiny/express/_run.py -> run_express), so these reactive.value()s are per-session and
    do NOT leak across users. Kept local to app.py for exactly that reason.
    """
    when = reactive.value(None)  # epoch time at which to commit, or None
    trigger = reactive.value(0)

    @reactive.effect
    def _on_change():
        value_fn()  # take a reactive dependency on the source
        when.set(time.time() + delay_secs)

    @reactive.effect
    def _on_timer():
        w = when()
        if w is None:
            return
        now = time.time()
        if now >= w:
            with reactive.isolate():
                when.set(None)
                trigger.set(trigger() + 1)
        else:
            reactive.invalidate_later(w - now)

    @reactive.calc
    def _debounced():
        trigger()  # re-evaluate only when a settled value is committed
        with reactive.isolate():
            return value_fn()

    return _debounced


# Debounced view of the range slider. The dashed selection lines still track the slider
# live (see update_position_plot_shapes); only the heavy recompute chain reads this.
pos_range_debounced = reactive_debounce(lambda: input.pos_range(), 0.3)


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
        return "—"

    selected = means.at["selected", col]
    if "full" in means.index and not pd.isna(means.at["full", col]):
        full = means.at["full", col]
    elif "unselected" in means.index and not pd.isna(means.at["unselected", col]):
        full = means.at["unselected", col]
    else:
        full = selected

    if as_int:
        return f"{int(full)} ({int(selected)} in selected)"
    return f"{full:.2f} ({selected:.2f} in selected)"

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


# These calcs must be defined before the UI block because they are referenced
# as arguments to @reactive.event decorators, which are evaluated at definition time.
@reactive.calc
def per_base_input() -> list[FileInfo] | None:
    """Resolve the active per-base file source.

    An explicit upload always wins; otherwise, once the user has clicked
    "Load example data", fall back to the bundled example TSV. Returning a
    FileInfo-shaped list keeps the rest of the pipeline agnostic about whether
    data came from an upload or the example. This is a per-session calc (no
    module-level state), preserving session isolation.
    """
    uploaded = input.per_base_file()
    if uploaded is not None:
        return uploaded
    if input.load_example() > 0:
        return [{"name": example_per_base_tsv.name, "datapath": str(example_per_base_tsv)}]
    return None


@reactive.calc
def reference_input() -> list[FileInfo] | None:
    """Resolve the active reference file source (upload wins, else the example)."""
    uploaded = input.reference_file()
    if uploaded is not None:
        return uploaded
    if input.load_example() > 0:
        return [
            {"name": example_reference_fasta.name, "datapath": str(example_reference_fasta)}
        ]
    return None


@reactive.calc
def parsed_reference() -> dict[str, dict | None] | None:
    """Parse reference file (FASTA or GenBank) and return dict with sequence string and
    feature objects (if present). If no feature objects are present, the "feature" value
    is None. Returns None if no file is uploaded or if the file is not a valid fasta or
    genbank file."""
    file = reference_input()
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
    file: list[FileInfo] | None = per_base_input()
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
        ref_seq = ref["sequence"]
        if input.reverse_complement():
            ref_seq = reverse_complement_sequence(ref_seq)
        data = align_ref_to_variants(data, ref_seq)
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
    # --- 1. Data input ---
    ui.input_file(
        "per_base_file",
        "Upload sequencing data",
        accept=[".csv", ".tsv"],
        multiple=False,
    )
    ui.input_action_button(
        "load_example",
        "Load example data",
        class_="btn-sm btn-outline-secondary",
    )
    ui.help_text("No file yet? Load a bundled example library to explore the app.")

    # --- 2. Reference input (optional) ---
    ui.input_file(
        "reference_file",
        "Upload reference sequence (optional)",
        accept=[".fa", ".fasta", ".gb", ".genbank", ".gbk"],
        multiple=False,
    )

    # Selectize for features parsed from a GenBank reference.
    # Only shown when a GenBank file provides features.
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

    # --- 3. Selection ---
    ui.input_slider(
        "pos_range",
        "Selected range:",
        min=0,
        max=100,
        value=[0, 100],
        step=1,
    )

    # --- 4. Metric display ---
    ui.input_checkbox_group(
        "data_series",
        "Display",
        {
            key: checkbox_with_tooltip(key, column_names_dict, column_tooltips)
            for key in plottable_series
        },
    )

    ui.hr()

    # --- 5. Advanced / reference alignment ---
    ui.input_switch("reverse_complement", "Reverse complement")
    ui.input_numeric("origin_shift", "Origin shift (bp)", value=0, min=0, step=1)

    ui.hr()
    "Export:"

    # Drives the conditional display of the download buttons (hidden output;
    # see #export_ready in styles.css). Downloads are gated until data exists.
    @render.text
    def export_ready():
        return "true" if not base_processed_data().empty else "false"

    with ui.panel_conditional("output.export_ready !== 'true'"):
        ui.help_text("Load data to enable downloads.")

    with ui.panel_conditional("output.export_ready === 'true'"):

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
                normalize=input.normalize_plot(),
            )
            try:
                yield fig.to_image(format="png")
            except Exception:
                ui.notification_show(
                    "PNG export requires the 'kaleido' package. Install with: pip install kaleido"
                )
                yield b""

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
                normalize=input.normalize_plot(),
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
                    input.normalize_plot,
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
                            normalize=input.normalize_plot(),
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
                        normalize=input.normalize_plot(),
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
                pos_range_debounced,
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
        ui.input_switch("normalize_plot", "Normalize series (0–1)")

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
@reactive.event(per_base_input)
def initialize_series_selection():
    """Initialize displayed metric selection once per upload (or example load)."""
    if per_base_input() is None:
        return
    ui.update_checkbox_group("data_series", selected=["entropy"])


@reactive.calc
def processed_per_base_file():
    """Apply range and feature selection to base processed data (cheap)."""
    data = base_processed_data()
    if data.empty:
        return pd.DataFrame()

    low, high = pos_range_debounced()
    data = update_per_base_df(data, [(low, high)])

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
    if per_base_input() is not None:
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
