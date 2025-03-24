# validation.py

import pandas as pd
import scipy.stats as stats

from shared import test_cols


def test_per_base_file(
    processed_data: pd.DataFrame,
    selected_means_df: pd.DataFrame,
    full_mean_df: pd.DataFrame,
) -> pd.DataFrame:
    """Function to test the metrics and return a dataframe with summary pass/fail results"""

    # If the dataframes are empty, return an empty dataframe
    if selected_means_df.empty or full_mean_df.empty:
        return pd.DataFrame()

    # Create a dataframe to store the results
    test_results = pd.DataFrame(
        columns=pd.Index(["Metric", "Test", "Result", "p_value", "t_stat"])
    )
    test_results.set_index("Metric", inplace=True)

    # Test the mean values
    test_results = test_mean_values(
        processed_data, selected_means_df, full_mean_df, test_results
    )

    # Test the standard deviation values
    # test_results = test_std_values(selected_means_df, full_mean_df, test_results)

    # Test the p-values
    # test_results = test_p_values(selected_means_df, full_mean_df, test_results)

    return test_results


def test_mean_values(
    processed_data: pd.DataFrame,
    selected_means_df: pd.DataFrame,
    full_mean_df: pd.DataFrame,
    test_results: pd.DataFrame,
) -> pd.DataFrame:
    """Function to test the mean values and return a dataframe with summary pass/fail results"""

    # Create a dataframe to store the results
    test_results = pd.DataFrame(
        columns=pd.Index(["Metric", "Test", "Result", "p_value", "t_stat"])
    )
    test_results.set_index("Metric", inplace=True)

    # Test the mean values
    for column in test_cols:
        # Get the mean values
        selected_values = processed_data[processed_data["is_selected"]][column]
        unselected_values = processed_data[~processed_data["is_selected"]][column]

        # Perform a t-test
        t_test_result = stats.ttest_ind(
            selected_values,
            unselected_values,
            equal_var=False,
        )
        t_stat, p_value = t_test_result.statistic, t_test_result.pvalue  # type: ignore

        # Perform the test
        if p_value < 0.05:
            test_results.loc[column, "Metric"] = column
            test_results.loc[column, "Test"] = "Compare means"
            test_results.loc[column, "Result"] = "Pass"
            test_results.loc[column, "p_value"] = p_value
            test_results.loc[column, "t_stat"] = t_stat
        else:
            test_results.loc[column, "Metric"] = column
            test_results.loc[column, "Test"] = "Compare means"
            test_results.loc[column, "Result"] = "Fail"
            test_results.loc[column, "p_value"] = p_value
            test_results.loc[column, "t_stat"] = t_stat

    print(test_results.columns)

    return test_results


def test_std_values(
    selected_means_df: pd.DataFrame,
    full_mean_df: pd.DataFrame,
    test_results: pd.DataFrame,
) -> pd.DataFrame:
    """Function to test the standard deviation values and return a dataframe with summary pass/fail results"""

    # Test the standard deviation values
    for column in selected_means_df.columns:
        # Get the standard deviation values
        selected_std = selected_means_df[column].std()
        full_std = full_mean_df[column].std()

        # Perform the test
        if selected_std == full_std:
            test_results.loc[column, "Test"] = "Standard Deviation Test"
            test_results.loc[column, "Result"] = "Pass"
        else:
            test_results.loc[column, "Test"] = "Standard Deviation Test"
            test_results.loc[column, "Result"] = "Fail"

    return test_results


def test_p_values(
    selected_means_df: pd.DataFrame,
    full_mean_df: pd.DataFrame,
    test_results: pd.DataFrame,
) -> pd.DataFrame:
    """Function to test the p-values and return a dataframe with summary pass/fail results"""

    # Test the p-values
    for column in selected_means_df.columns:
        # Get the p-values
        selected_p = selected_means_df[column].p_value
        full_p = full_mean_df[column].p_value

        # Perform the test
        if selected_p == full_p:
            test_results.loc[column, "Test"] = "P-Value Test"
            test_results.loc[column, "Result"] = "Pass"
        else:
            test_results.loc[column, "Test"] = "P-Value Test"
            test_results.loc[column, "Result"] = "Fail"

    return test_results
