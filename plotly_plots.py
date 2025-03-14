import plotly.graph_objects as go
import plotly.express as px
import pandas as pd


def base_position_vs_value_plot_plotly(
    per_base_df: pd.DataFrame,
    displayed_fields: list,
    range: list,
    selected_range_low: int,
    selected_range_high: int,
    last_selected_series: str,
    show_means: bool,
) -> go.Figure:
    if per_base_df.empty:
        return go.Figure()
    if not displayed_fields:
        return go.Figure()
    if not range:
        range = [0, per_base_df["pos"].max()]
    if "pos" not in per_base_df.columns:
        return go.Figure()

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
        selected_mean = per_base_df.loc[
            per_base_df["pos"].between(selected_range_low, selected_range_high),
            last_selected_series,
        ].mean()
        all_mean = per_base_df[last_selected_series].mean()

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
) -> go.Figure:
    if per_base_df.empty:
        return go.Figure()
    if selected_series not in per_base_df.columns:
        return go.Figure()
    if selected_series not in column_names_dict:
        return go.Figure()

    fig = px.violin(
        per_base_df,
        x=selected_series,
        color="is_selected",
        box=False,
        template="simple_white",
        violinmode="overlay",
        points="all",
    )

    # Plot selected vs unselected ranges
    fig.update_layout(
        title=f"Distribution of {column_names_dict[selected_series]}",
        yaxis_title="",
        xaxis_title=column_names_dict[selected_series],
    )

    fig.update_traces(
        meanline_visible=True,
    )  # scale violin plot area with total count

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
