import pandas as pd

from shiny import reactive
from shiny.express import input, render, ui
from shiny.types import FileInfo

from shinywidgets import render_plotly

from Bio import SeqIO

from faicons import icon_svg

from plotly_plots import base_position_vs_value_plot_plotly
from plots import position_vs_value_plot, violin_plot
from process_data import process_per_base_file, update_per_base_df
from shared import column_colors_dict, column_names_dict, column_tooltips, tabular_cols

# from input_checkbox_group_tooltips import input_checkbox_group_tooltips

ui.page_opts(title="DIMPLE quick QC", fillable=True)
ui.include_css("./styles.css")


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

                @render_plotly
                def plotly_position_plot():
                    pos_plot = base_position_vs_value_plot_plotly(
                        processed_per_base_file(),
                        input.data_series(),
                        [plot_range_low.get(), plot_range_high.get()],
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

                    return render.DataGrid(
                        pd.DataFrame(processed_per_base_file()[tabular_cols]),
                        filters=False,
                    )

        # Bottom 2D plots
        with ui.card(height=400):

            ui.card_header("Distributions")

            @render.plot
            @reactive.event(
                selected_range_low, selected_range_high, last_selected_series
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
            def avg_effective_entropy():
                df = processed_per_base_file()
                if df.empty:
                    return "0"

                selected_range = range(
                    selected_range_low.get(), selected_range_high.get()
                )

                full_avg = df["effective_entropy"].mean(skipna=True)
                range_avg = df.loc[df["pos"].isin(selected_range)][
                    "effective_entropy"
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

    data = process_per_base_file(parsed_per_base_file())

    # Set the sequence length
    sequence_length.set(max(data["pos"]))

    # Update the selected range and plotting range
    selected_range_high.set(sequence_length.get())
    plot_range_high.set(sequence_length.get())

    return data


# Reactive effects
@reactive.effect
@reactive.event(selected_range_low, selected_range_high)
def update_data_selected_range():
    # Calculate data that is dependent on the range selection

    if processed_per_base_file().empty:
        return

    update_per_base_df(
        processed_per_base_file(), selected_range_low.get(), selected_range_high.get()
    )


# Update alignment when a reference sequence is uploaded
@reactive.effect
@reactive.event(input.reference_fasta)
def update_alignment():
    return


# Dynamically update the max_pos input field when sequence_length updates
@reactive.effect
def update_max_pos():
    if not processed_per_base_file().empty:
        ui.update_numeric("max_pos", value=sequence_length.get())


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
    ui.update_numeric("min_pos", value=selected_range_low.get())


@reactive.effect
@reactive.event(input.max_pos)
def updated_max_pos():
    selected_range_high.set(input.max_pos())
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
