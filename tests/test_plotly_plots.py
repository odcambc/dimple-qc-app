"""Tests for Plotly plotting helpers."""

import pandas as pd

from plotly_plots import base_position_vs_value_plot_plotly
from process_data import process_per_base_file, update_mean_values_per_base, update_per_base_df


class TestBasePositionVsValuePlot:
    def test_show_means_uses_flattened_mean_columns(
        self, minimal_per_base_df: pd.DataFrame
    ) -> None:
        processed = process_per_base_file(minimal_per_base_df, False)
        updated = update_per_base_df(processed, [(2, 4)])
        means = update_mean_values_per_base(updated)

        fig = base_position_vs_value_plot_plotly(
            updated,
            means,
            ["entropy"],
            [0, int(updated["pos"].max())],
            2,
            4,
            "entropy",
            True,
        )
        assert fig is not None
        assert len(fig.data) > 0

    def test_show_means_empty_mean_df_no_crash(
        self, minimal_per_base_df: pd.DataFrame
    ) -> None:
        processed = process_per_base_file(minimal_per_base_df, False)
        fig = base_position_vs_value_plot_plotly(
            processed,
            pd.DataFrame(),
            ["entropy"],
            [0, 10],
            0,
            10,
            "entropy",
            True,
        )
        assert fig is not None

    def test_feature_regions_adds_vrects(
        self, minimal_per_base_df: pd.DataFrame
    ) -> None:
        processed = process_per_base_file(minimal_per_base_df, False)
        regions = [
            {"start": 1, "end": 3, "label": "CDS1", "color": "rgba(100,149,237,0.15)"},
            {"start": 4, "end": 6, "label": "Promoter"},
        ]
        fig = base_position_vs_value_plot_plotly(
            processed,
            pd.DataFrame(),
            ["entropy"],
            [0, 6],
            0,
            6,
            "entropy",
            False,
            feature_regions=regions,
        )
        # vrects are stored in layout.shapes
        shapes = fig.layout.shapes
        vrects = [s for s in shapes if s.type == "rect"]
        assert len(vrects) == 2
        assert vrects[0].x0 == 1
        assert vrects[0].x1 == 3
        assert vrects[1].x0 == 4
        assert vrects[1].x1 == 6

    def test_no_feature_regions_no_extra_shapes(
        self, minimal_per_base_df: pd.DataFrame
    ) -> None:
        processed = process_per_base_file(minimal_per_base_df, False)
        fig = base_position_vs_value_plot_plotly(
            processed,
            pd.DataFrame(),
            ["entropy"],
            [0, 6],
            0,
            6,
            "entropy",
            False,
        )
        # Only shapes should be from vlines (selected range markers)
        shapes = fig.layout.shapes
        vrects = [s for s in shapes if s.type == "rect" and hasattr(s, "fillcolor") and s.fillcolor]
        assert len(vrects) == 0
