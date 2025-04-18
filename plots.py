import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D


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
