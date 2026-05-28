"""Regression tests that protect the app's session-isolation guarantee.

Multiple concurrent users share one Python process. Each user gets their own
Shiny session (and therefore their own reactive graph), but any module-level
mutable state would be shared across sessions and could leak user data.

These tests assert two complementary properties:

1. Structural: no module in the app introduces the known leak vectors
   (module-level reactive.Value, module-level Plotly Figure, lru_cache on
   data-handling functions).
2. Behavioral: the per-session data pipeline functions are pure -- running
   them on user A's data does not affect what they return for user B, and
   they do not mutate their inputs in place.

If these tests fail, the most likely cause is a recent change that
re-introduced shared mutable state. See CLAUDE.md for the project's
no-module-level-mutable-state rule.
"""

from __future__ import annotations

import importlib
from types import ModuleType

import pandas as pd
import plotly.graph_objects as go
import pytest
from shiny import reactive

from process_data import process_per_base_file

# Modules that participate in the per-session data pipeline. app.py is excluded
# because in Shiny Express the entire file is re-executed per session, so its
# module-scope reactive primitives are session-scoped by design.
APP_MODULES = [
    "evaluate_data",
    "per_base_io",
    "plotly_plots",
    "process_data",
    "process_reference",
    "shared",
    "validation",
]


@pytest.fixture(params=APP_MODULES)
def app_module(request: pytest.FixtureRequest) -> ModuleType:
    return importlib.import_module(request.param)


class TestNoSharedMutableState:
    """Structural assertions. Each catches one previously-considered leak vector."""

    def test_no_module_level_reactive_value(self, app_module: ModuleType) -> None:
        for name, obj in vars(app_module).items():
            assert not isinstance(obj, reactive.Value), (
                f"{app_module.__name__}.{name} is a module-level reactive.Value, "
                "which is shared across sessions. Use @reactive.calc inside a "
                "session-scoped function instead."
            )

    def test_no_module_level_plotly_figures(self, app_module: ModuleType) -> None:
        for name, obj in vars(app_module).items():
            if name.startswith("__"):
                continue
            assert not isinstance(obj, go.Figure), (
                f"{app_module.__name__}.{name} is a module-level go.Figure. "
                "Plotly figures are mutated in place by methods like add_trace "
                "and update_layout; return them from a factory function instead."
            )

    def test_no_lru_cache_at_module_scope(self, app_module: ModuleType) -> None:
        for name, obj in vars(app_module).items():
            if not callable(obj) or not hasattr(obj, "cache_info"):
                continue
            pytest.fail(
                f"{app_module.__name__}.{name} appears to be functools.lru_cache "
                "or @cache decorated. Caching across sessions on functions that "
                "take user-uploaded data leaks results between users."
            )


class TestPipelineIsPure:
    """The pipeline functions must be pure -- this is the second line of defense
    against cross-session leaks if shared state is ever reintroduced upstream."""

    def test_process_per_base_file_is_deterministic(
        self,
        minimal_per_base_df: pd.DataFrame,
        variant_region_per_base_df: pd.DataFrame,
    ) -> None:
        # Simulate user A, then user B, then user A again. A's two runs must match.
        a1 = process_per_base_file(minimal_per_base_df.copy(), False)
        _ = process_per_base_file(variant_region_per_base_df.copy(), False)
        a2 = process_per_base_file(minimal_per_base_df.copy(), False)

        pd.testing.assert_frame_equal(a1, a2)

    def test_process_per_base_file_does_not_mutate_input(
        self, minimal_per_base_df: pd.DataFrame
    ) -> None:
        before = minimal_per_base_df.copy()
        _ = process_per_base_file(minimal_per_base_df, False)
        pd.testing.assert_frame_equal(before, minimal_per_base_df)
