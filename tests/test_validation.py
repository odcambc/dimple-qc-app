import pandas as pd

from validation import validate_per_base_file


class TestValidatePerBaseFile:
    def test_valid_file(self, valid_df, patch_validation_ui):
        assert validate_per_base_file(valid_df) is True

    def test_empty_dataframe(self, patch_validation_ui):
        assert validate_per_base_file(pd.DataFrame()) is False

    def test_missing_columns(self, patch_validation_ui):
        df = pd.DataFrame({"pos": [1], "ref": ["A"]})
        assert validate_per_base_file(df) is False
        patch_validation_ui.assert_called_once()

    def test_extra_columns_ok(self, valid_df, patch_validation_ui):
        valid_df["extra_col"] = [42]
        assert validate_per_base_file(valid_df) is True
