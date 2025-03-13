import plotly.graph_objects as go
import plotly.express as px
import pandas as pd


def base_position_vs_value_plot_plotly(
    per_base_df: pd.DataFrame, displayed_fields: list, range: list
) -> go.Figure:
    if per_base_df.empty:
        return go.Figure()
    if not displayed_fields:
        return go.Figure()
    if not range:
        range = [0, per_base_df["pos"].max()]
    if "pos" not in per_base_df.columns:
        return go.Figure()

    fig = px.scatter()
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
    return fig
