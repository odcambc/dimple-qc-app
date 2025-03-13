import plotly.express as px
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D


def position_vs_value_plot_plotly(per_base_df, displayed_fields, range):
    vals = per_base_df[stat]
    vals = vals[~vals.isnull()]
    fig = ff.create_distplot(
        [vals],
        ["Overall"],
        rug_text=[careers_df["player_name"]],
        colors=["black"],
        show_hist=False,
    )
    # Clean up some defaults (1st trace is the density plot, 2nd is the rug plot)
    fig.data[0].hoverinfo = "none"
    fig.data[0].showlegend = False
    fig.data[1].hoverinfo = "text+x"
    fig.data[1].customdata = careers_df["person_id"]
    # Use height of the density plot to inform the vertical lines
    ymax = fig.data[0].y.max()
    # Arrange rows from highest to lowest value so that legend order is correct
    stats_df = stats_df.sort_values(stat, ascending=False)
    # Add vertical lines for each player
    for _, row in stats_df.iterrows():
        x = row[stat]
        fig.add_scatter(
            x=[x, x],
            y=[0, ymax],
            mode="lines",
            name=players_dict[row["person_id"]],
            line=dict(color=row["color"], width=1),
            hoverinfo="x+name",
        )

    fig.update_layout(
        hovermode="x",
        xaxis=dict(title=stat + " per game (career average)", hoverformat=".1f"),
        legend=dict(orientation="h", y=1.03, yanchor="bottom", x=0.5, xanchor="center"),
    )

    # Convert Figure to FigureWidget so we can add click events
    fig = go.FigureWidget(fig.data, fig.layout)
    fig.data[1].on_click(on_rug_click)

    return fig


def position_vs_value_plot(
    per_base_df: pd.DataFrame,
    data_series: list,
    selected_range_low: int,
    selected_range_high: int,
    plot_range_low: int,
    plot_range_high: int,
    last_selected_series: str,
    column_colors_dict: dict,
    column_names_dict: dict,
    show_means: bool,
):
    if per_base_df.empty:
        return plt.figure()

    pos_plot = sns.scatterplot()

    # Plot vertical lines at selected range
    pos_plot.axvline(selected_range_low, color="black", linestyle="--")
    pos_plot.axvline(selected_range_high, color="black", linestyle="--")

    # Plot data series
    for series in data_series:
        sns.scatterplot(
            data=per_base_df,
            x="pos",
            y=series,
            color=column_colors_dict[series],
            edgecolors="black",
        )

    # Set plot limits and labels
    pos_plot.set_xlim(plot_range_low, plot_range_high)
    pos_plot.set_xlabel("Position")
    pos_plot.set_ylabel(column_names_dict[last_selected_series])

    # If show_means is selected, draw horizontal lines at means of last
    # selected series for the full range and selected range

    if show_means:
        if data_series:
            full_mean = per_base_df[last_selected_series].mean(skipna=True)
            selected_mean = per_base_df.loc[
                per_base_df["pos"].between(
                    selected_range_low,
                    selected_range_high,
                ),
                last_selected_series,
            ].mean(skipna=True)

            pos_plot.axhline(
                float(full_mean),
                color=column_colors_dict[last_selected_series],
                linestyle="solid",
            )
            pos_plot.axhline(
                float(selected_mean),
                color=column_colors_dict[last_selected_series],
                linestyle="dashed",
            )

            # Add legend with correct labels

            legend_handles = [
                Line2D(
                    [0],
                    [0],
                    color=column_colors_dict[last_selected_series],
                    linestyle="solid",
                    label="Full range",
                ),
                Line2D(
                    [0],
                    [0],
                    color=column_colors_dict[last_selected_series],
                    linestyle="dashed",
                    label="Selected range",
                ),
            ]
            pos_plot.legend(
                handles=legend_handles,
                loc="upper right",
            )

        return pos_plot


def violin_plot(
    per_base_df: pd.DataFrame, last_selected_series: str, column_names_dict: dict
):

    if per_base_df.empty:
        return plt.figure()
    if last_selected_series not in per_base_df:
        return plt.figure()

    # Now refactor the df to allow for easy plotting.
    # We want to plot the distributions of the selected series
    # for the full range of the data vs the selected range.
    # We can use the "is_selected" column to split the data.

    plotting_df = per_base_df.melt(
        id_vars=["is_selected"],
        value_vars=last_selected_series,
        var_name="series",
        value_name="value",
    )

    violin_plot = sns.violinplot(
        data=plotting_df,
        x="value",
        hue="is_selected",
        split=True,
        linewidth=1,
        linecolor="black",
        fill=True,
    )

    violin_plot.set_xlabel(column_names_dict[last_selected_series])
    violin_plot.set_ylabel("Density")
    violin_plot.set_title(f"Distribution of {column_names_dict[last_selected_series]}")
    violin_plot.legend(
        title="Selected range",
        loc="upper right",
    )

    return violin_plot
