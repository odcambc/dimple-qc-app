import pandas as pd
import pytest
from unittest.mock import patch

from validation import validate_per_base_file, expected_columns


@pytest.fixture
def valid_df():
    """DataFrame with all required columns."""
    return pd.DataFrame(
        {col: [1] for col in expected_columns}
    )


class TestValidatePerBaseFile:
    def test_valid_file(self, valid_df):
        # Mock ui.notification_show since we don't have a Shiny session
        with patch("validation.ui.notification_show"):
            assert validate_per_base_file(valid_df) is True

    def test_empty_dataframe(self):
        with patch("validation.ui.notification_show"):
            assert validate_per_base_file(pd.DataFrame()) is False

    def test_missing_columns(self):
        df = pd.DataFrame({"pos": [1], "ref": ["A"]})
        with patch("validation.ui.notification_show") as mock_notify:
            assert validate_per_base_file(df) is False
            mock_notify.assert_called_once()

    def test_extra_columns_ok(self, valid_df):
        valid_df["extra_col"] = [42]
        with patch("validation.ui.notification_show"):
            assert validate_per_base_file(valid_df) is True
