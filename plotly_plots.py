import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

# Default empty figure
empty_fig = go.Figure().update_layout(template="simple_white")


def base_position_vs_value_plot_plotly(
    per_base_df: pd.DataFrame,
    mean_values: pd.DataFrame,
    displayed_fields: list,
    range: list,
    selected_range_low: int,
    selected_range_high: int,
    last_selected_series: str,
    show_means: bool,
) -> go.Figure:
    if per_base_df.empty:
        return empty_fig
    if not displayed_fields:
        return empty_fig
    if not range:
        range = [0, per_base_df["pos"].max()]
    if "pos" not in per_base_df.columns:
        return empty_fig

    fig = px.scatter(template="simple_white")
    for field in displayed_fields:
        fig.add_trace(
            go.Scatter(
                x=per_base_df["pos"],
                y=per_base_df[field],
                mode="markers",
                name=field,
            )
        )
        fig.update_layout(
            title="Position vs Value",
            xaxis_title="Position",
            yaxis_title="Value",
            xaxis=dict(range=range),
        )

    fig.add_vline(
        x=selected_range_low, line_width=1, line_dash="dash", line_color="black"
    )
    fig.add_vline(
        x=selected_range_high, line_width=1, line_dash="dash", line_color="black"
    )

    if show_means:
        selected_mean = mean_values.at[True, (last_selected_series, "mean")]
        all_mean = mean_values.at[False, (last_selected_series, "mean")]

        fig.add_hline(
            y=selected_mean,
            line_width=1,
            line_dash="dash",
            line_color="black",
            annotation_text=f"Selected mean: {selected_mean:.2f}",
        )
        fig.add_hline(
            y=all_mean,
            line_width=1,
            line_dash="solid",
            line_color="black",
            annotation_text=f"Full mean: {all_mean:.2f}",
        )

    return fig


def distribution_violin_plot_plotly(
    per_base_df: pd.DataFrame,
    selected_series: str,
    column_names_dict: dict,
    column_colors_dict: dict = None,  # type: ignore
) -> go.Figure:
    """
    Create a split violin plot showing the distribution of values for selected vs unselected ranges.
    The violin is split horizontally with the top half showing one series and the bottom half showing the other.

    Args:
        per_base_df: DataFrame containing the per-base data
        selected_series: Column name to plot
        column_names_dict: Dictionary mapping column names to display names
        column_colors_dict: Optional dictionary mapping column names to colors

    Returns:
        Plotly figure object
    """
    # Handle empty or invalid data
    if per_base_df.empty:
        return empty_fig
    if selected_series not in per_base_df.columns:
        return empty_fig
    if selected_series not in column_names_dict:
        return empty_fig

    # Create a copy of the dataframe to avoid modifying the original

    fig = go.Figure()

    # Get data for selected and unselected ranges
    selected_data = per_base_df[per_base_df["is_selected"]][selected_series]

    unselected_data = per_base_df[~per_base_df["is_selected"]][selected_series]

    selected_color = "#1F77B4"  # Blue
    unselected_color = "#FF7F0E"  # Orange

    # Use custom colors if provided
    if column_colors_dict and selected_series in column_colors_dict:
        base_color = column_colors_dict[selected_series]
        # Create a lighter shade for the selected range
        selected_color = base_color
        # Keep the unselected color as is or create a complementary color

    # Add the top half (selected range)
    fig.add_trace(
        go.Violin(
            x=selected_data,
            y0="Distribution",
            name="Selected Range",
            side="positive",  # Only show the top half
            line_color="black",
            fillcolor=selected_color,
            opacity=0.6,
            meanline_visible=True,
            meanline_color="white",
            meanline_width=2,
            box_visible=True,
            box_width=0.1,
            box_fillcolor="rgba(255,255,255,0.5)",
            points="outliers",
            jitter=0,
            marker=dict(color=selected_color),
            alignmentgroup="selected",  # Ensure violins are aligned
            offsetgroup="selected",  # Ensure violins are aligned
            legendgroup="selected",
        )
    )

    # Add the bottom half (unselected/full sequence)
    fig.add_trace(
        go.Violin(
            x=unselected_data,
            y0="Distribution",
            name="Full Sequence",
            side="negative",  # Only show the bottom half
            line_color="black",
            fillcolor=unselected_color,
            opacity=0.6,
            meanline_visible=True,
            meanline_color="white",
            meanline_width=2,
            box_visible=True,
            box_width=0.1,
            box_fillcolor="rgba(255,255,255,0.5)",
            points="outliers",
            jitter=0,
            marker=dict(color=unselected_color),
            alignmentgroup="unselected",  # Ensure violins are aligned
            offsetgroup="unselected",  # Ensure violins are aligned
            legendgroup="unselected",
        )
    )

    # Need to set orientation to horizontal
    fig.update_traces(meanline_visible=False, orientation="h")

    fig.update_layout(
        title={
            "text": f"Distribution of {column_names_dict[selected_series]}",
            "y": 0.95,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
            "font": {"size": 18},
        },
        yaxis_title="",
        xaxis_title={"text": column_names_dict[selected_series], "font": {"size": 14}},
        legend_title_text="",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=50, r=50, t=80, b=50),
        template="simple_white",
        violingap=0.5,
        violingroupgap=0.5,
        violinmode="overlay",
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            showticklabels=False,
        ),
    )

    return fig


def distribution_histogram_plot_plotly(
    per_base_df: pd.DataFrame,
    selected_series: str,
    column_names_dict: dict,
) -> go.Figure:
    if per_base_df.empty:
        return go.Figure()
    if selected_series not in per_base_df.columns:
        return go.Figure()
    if selected_series not in column_names_dict:
        return go.Figure()

    fig = px.histogram(
        per_base_df,
        x=selected_series,
        color="is_selected",
        opacity=0.8,
        nbins=20,
        barmode="overlay",
    )

    fig.update_layout(
        title=f"Distribution of {column_names_dict[selected_series]}",
        xaxis_title=column_names_dict[selected_series],
        yaxis_title="Count",
    )
    return fig
